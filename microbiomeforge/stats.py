"""İstatistik motoru — çeşitlilik, ordinasyon, grup karşılaştırması, diferansiyel bolluk.

Rapor kuralı: her bulgu istatistikle desteklenir (test + etki büyüklüğü + p/q). Tüm
hesaplar numpy+scipy ile native yapılır (harici R'a bağımlı değil).

Kapsam:
  * Alfa çeşitlilik: Shannon, Simpson, gözlenen zenginlik (richness)
  * Alfa grup karşılaştırması: Mann-Whitney U (2 grup) / Kruskal-Wallis (>2)
  * Beta çeşitlilik: Bray-Curtis mesafe matrisi
  * Ordinasyon: PCoA (klasik MDS)
  * Grup ayrımı: PERMANOVA (permütasyon testi, pseudo-F)
  * Diferansiyel bolluk: CLR dönüşümü + taxon-başı test + BH-FDR (q-değeri)

Girdi: M4'ten taxon×örnek bağıl-bolluk matrisi + örnek→grup eşlemesi.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats as sps


# --------------------------------------------------------------------------- #
# Alfa çeşitlilik
# --------------------------------------------------------------------------- #
def alpha_diversity(matrix: pd.DataFrame) -> pd.DataFrame:
    """taxon×örnek matrisinden örnek başına alfa çeşitlilik metrikleri."""
    samples = [c for c in matrix.columns if c != "taxon"]
    rows = []
    for s in samples:
        p = matrix[s].to_numpy(dtype=float)
        p = p[p > 0]
        total = p.sum()
        if total <= 0:
            rows.append({"sample": s, "shannon": 0.0, "simpson": 0.0, "richness": 0})
            continue
        p = p / total
        shannon = float(-np.sum(p * np.log(p)))
        simpson = float(1.0 - np.sum(p ** 2))
        rows.append({
            "sample": s,
            "shannon": round(shannon, 4),
            "simpson": round(simpson, 4),
            "richness": int((matrix[s] > 0).sum()),
        })
    return pd.DataFrame(rows)


@dataclass
class GroupTest:
    metric: str
    test: str
    statistic: float
    p_value: float
    n_groups: int
    effect: float  # 2 grup: rank-biserial; >2: epsilon^2 yaklaşık


def compare_alpha(alpha_df: pd.DataFrame, groups: dict[str, str], metric: str = "shannon") -> GroupTest:
    """Alfa çeşitliliğini gruplar arasında karşılaştırır (parametrik olmayan)."""
    df = alpha_df.copy()
    df["group"] = df["sample"].map(groups)
    df = df.dropna(subset=["group"])
    grouped = [g[metric].to_numpy(dtype=float) for _, g in df.groupby("group")]
    n_groups = len(grouped)
    if n_groups < 2 or any(len(x) == 0 for x in grouped):
        return GroupTest(metric, "none", float("nan"), float("nan"), n_groups, float("nan"))
    if n_groups == 2:
        u, p = sps.mannwhitneyu(grouped[0], grouped[1], alternative="two-sided")
        n1, n2 = len(grouped[0]), len(grouped[1])
        effect = 1.0 - (2.0 * u) / (n1 * n2)  # rank-biserial korelasyon
        return GroupTest(metric, "Mann-Whitney U", float(u), float(p), 2, round(float(effect), 4))
    h, p = sps.kruskal(*grouped)
    n = sum(len(x) for x in grouped)
    eps2 = (h - n_groups + 1) / (n - n_groups) if n > n_groups else float("nan")
    return GroupTest(metric, "Kruskal-Wallis", float(h), float(p), n_groups, round(float(eps2), 4))


# --------------------------------------------------------------------------- #
# Beta çeşitlilik + ordinasyon + PERMANOVA
# --------------------------------------------------------------------------- #
def bray_curtis_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    """Örnekler arası Bray-Curtis mesafe matrisi."""
    samples = [c for c in matrix.columns if c != "taxon"]
    X = matrix[samples].to_numpy(dtype=float).T  # örnek × taxon
    n = len(samples)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.abs(X[i] - X[j]).sum()
            den = (X[i] + X[j]).sum()
            d = num / den if den > 0 else 0.0
            D[i, j] = D[j, i] = d
    return pd.DataFrame(D, index=samples, columns=samples)


def pcoa(distance: pd.DataFrame, n_axes: int = 2) -> pd.DataFrame:
    """Klasik MDS (PCoA) koordinatları + eksen başına açıklanan varyans."""
    D = distance.to_numpy(dtype=float)
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    pos = eigvals > 0
    coords = eigvecs[:, :n_axes] * np.sqrt(np.abs(eigvals[:n_axes]))
    explained = eigvals[:n_axes] / eigvals[pos].sum() if pos.any() else np.zeros(n_axes)
    out = pd.DataFrame(coords, columns=[f"PCo{i+1}" for i in range(n_axes)], index=distance.index)
    out.attrs["explained"] = [round(float(e), 4) for e in explained]
    return out.reset_index().rename(columns={"index": "sample"})


@dataclass
class Permanova:
    pseudo_f: float
    p_value: float
    permutations: int
    n_groups: int


def permanova(distance: pd.DataFrame, groups: dict[str, str], permutations: int = 999,
              seed: int = 42) -> Permanova:
    """PERMANOVA: grupların beta-çeşitlilikte ayrılıp ayrılmadığını test eder."""
    samples = list(distance.index)
    labels = np.array([groups.get(s) for s in samples])
    mask = labels != None  # noqa: E711
    samples = list(np.array(samples)[mask])
    labels = labels[mask]
    D = distance.loc[samples, samples].to_numpy(dtype=float)
    uniq = np.unique(labels)
    n_groups = len(uniq)
    if n_groups < 2:
        return Permanova(float("nan"), float("nan"), permutations, n_groups)

    def pseudo_f(lab):
        n = len(lab)
        total_ss = (D ** 2).sum() / (2 * n)
        within_ss = 0.0
        for g in np.unique(lab):
            idx = np.where(lab == g)[0]
            ng = len(idx)
            if ng > 1:
                sub = D[np.ix_(idx, idx)]
                within_ss += (sub ** 2).sum() / (2 * ng)
        among_ss = total_ss - within_ss
        a, b = n_groups - 1, n - n_groups
        return (among_ss / a) / (within_ss / b) if within_ss > 0 and b > 0 else float("nan")

    observed = pseudo_f(labels)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(permutations):
        perm = rng.permutation(labels)
        if pseudo_f(perm) >= observed:
            count += 1
    p = (count + 1) / (permutations + 1)
    return Permanova(round(float(observed), 4), round(float(p), 4), permutations, n_groups)


# --------------------------------------------------------------------------- #
# Diferansiyel bolluk (CLR + taxon-başı test + BH-FDR)
# --------------------------------------------------------------------------- #
def _clr(matrix_vals: np.ndarray, pseudocount: float = 1e-6) -> np.ndarray:
    """Örnek × taxon → centered log-ratio."""
    x = matrix_vals + pseudocount
    logx = np.log(x)
    return logx - logx.mean(axis=1, keepdims=True)


def benjamini_hochberg(pvals: np.ndarray) -> np.ndarray:
    """BH-FDR düzeltilmiş q-değerleri."""
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


def differential_abundance(matrix: pd.DataFrame, groups: dict[str, str]) -> pd.DataFrame:
    """CLR dönüşümü sonrası taxon-başı grup karşılaştırması + BH-FDR.

    2 grup: Mann-Whitney U; >2 grup: Kruskal-Wallis. Etki: gruplar arası ortalama
    CLR farkı (2 grup) / grup-içi varyansa göre H (>2). Dönüş q-değerine göre sıralı.
    """
    samples = [c for c in matrix.columns if c != "taxon"]
    labels = np.array([groups.get(s) for s in samples])
    keep = labels != None  # noqa: E711
    samples = list(np.array(samples)[keep])
    labels = labels[keep]
    taxa = matrix["taxon"].tolist()

    X = matrix[samples].to_numpy(dtype=float).T  # örnek × taxon
    clr = _clr(X)
    uniq = list(dict.fromkeys(labels))
    two = len(uniq) == 2

    pvals, effects = [], []
    for k in range(clr.shape[1]):
        col = clr[:, k]
        gvals = [col[labels == g] for g in uniq]
        try:
            if two:
                _, p = sps.mannwhitneyu(gvals[0], gvals[1], alternative="two-sided")
                effect = float(np.mean(gvals[1]) - np.mean(gvals[0]))  # CLR log-fark
            else:
                _, p = sps.kruskal(*gvals)
                effect = float(np.var([np.mean(g) for g in gvals]))
        except ValueError:
            p, effect = 1.0, 0.0
        pvals.append(p)
        effects.append(effect)

    q = benjamini_hochberg(np.array(pvals))
    out = pd.DataFrame({
        "taxon": taxa,
        "effect_clr": np.round(effects, 4),
        "p_value": np.round(pvals, 5),
        "q_value": np.round(q, 5),
    })
    out["significant"] = out["q_value"] < 0.05
    return out.sort_values("q_value").reset_index(drop=True)
