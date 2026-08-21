"""Şekil üretimi — yayın kalitesinde, ZORUNLU numaralı + isimli.

Rapor kuralı: her şekil "Şekil N. <başlık>" biçiminde numaralı ve isimlidir; metinde
bu numarayla atıfta bulunulur. Şekiller ≥300 DPI PNG olarak (rapora gömmek için) ve
ayrıca tek bir çok-sayfalı `06_figures.pdf` olarak üretilir.

`FigureManifest` her şeklin (numara, başlık, png yolu) kaydını tutar; M7 rapor
birleştirici bunu kullanarak şekilleri gömer ve numarayla atıf yapar.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # başsız (GUI yok)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

# Yayın stili.
plt.rcParams.update({
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.titlesize": 12,
    "axes.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.autolayout": True,
})


@dataclass
class FigureEntry:
    number: int
    title: str
    png_path: str


@dataclass
class FigureManifest:
    entries: list[FigureEntry] = field(default_factory=list)

    def ref(self, number: int) -> str:
        return f"Şekil {number}"


class Figures:
    """Numaralı şekil üreticisi. Her `add_*` çağrısı numara atar ve PNG yazar."""

    def __init__(self, figdir: str | Path):
        self.figdir = Path(figdir)
        self.figdir.mkdir(parents=True, exist_ok=True)
        self._n = 0
        self.manifest = FigureManifest()

    def _save(self, fig, title: str) -> int:
        self._n += 1
        num = self._n
        fig.suptitle(f"Şekil {num}. {title}", fontsize=12, y=1.02)
        png = self.figdir / f"figure_{num:02d}.png"
        fig.savefig(png, bbox_inches="tight")
        plt.close(fig)
        self.manifest.entries.append(FigureEntry(num, title, str(png)))
        return num

    # ---- Şekil türleri ---- #
    def raw_qc(self, qc_df: pd.DataFrame, title: str = "Ham okuma kalitesi (örnek başına)") -> int:
        raw = qc_df[qc_df["stage"] == "raw"] if "stage" in qc_df else qc_df
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 4))
        ax1.bar(raw["sample"], raw["n_reads"], color="#4C72B0")
        ax1.set_ylabel("Okuma sayısı")
        ax1.set_title("Okuma derinliği")
        ax1.tick_params(axis="x", rotation=45)
        ax2.bar(raw["sample"], raw["mean_quality"], color="#55A868")
        ax2.set_ylabel("Ortalama Phred Q")
        ax2.set_title("Kalite")
        ax2.tick_params(axis="x", rotation=45)
        return self._save(fig, title)

    def abundance_bar(self, matrix: pd.DataFrame, top_n: int = 10,
                      title: str = "Taksonomik bileşim (en bol taksonlar)") -> int:
        samples = [c for c in matrix.columns if c != "taxon"]
        m = matrix.copy()
        m["_mean"] = m[samples].mean(axis=1)
        m = m.sort_values("_mean", ascending=False)
        top = m.head(top_n)
        other = m[samples].iloc[top_n:].sum()
        fig, ax = plt.subplots(figsize=(8, 5))
        bottom = np.zeros(len(samples))
        cmap = plt.get_cmap("tab20")
        for i, (_, row) in enumerate(top.iterrows()):
            vals = row[samples].to_numpy(dtype=float)
            ax.bar(samples, vals, bottom=bottom, label=str(row["taxon"]), color=cmap(i % 20))
            bottom += vals
        other = other.to_numpy(dtype=float)
        if (other > 0).any():
            ax.bar(samples, other, bottom=bottom, label="Diğer", color="#BBBBBB")
        ax.set_ylabel("Bağıl bolluk")
        ax.tick_params(axis="x", rotation=45)
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        return self._save(fig, title)

    def alpha_box(self, alpha_df: pd.DataFrame, groups: dict[str, str], metric: str = "shannon",
                  title: str | None = None) -> int:
        df = alpha_df.copy()
        df["group"] = df["sample"].map(groups)
        df = df.dropna(subset=["group"])
        order = list(dict.fromkeys(df["group"]))
        data = [df[df["group"] == g][metric].to_numpy() for g in order]
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.boxplot(data, labels=order, showmeans=True)
        ax.set_ylabel(metric.capitalize())
        ax.set_xlabel("Grup")
        return self._save(fig, title or f"Alfa çeşitlilik ({metric}) grup karşılaştırması")

    def pcoa_scatter(self, coords: pd.DataFrame, groups: dict[str, str],
                     title: str = "Beta çeşitlilik ordinasyonu (PCoA, Bray-Curtis)") -> int:
        df = coords.copy()
        df["group"] = df["sample"].map(groups)
        explained = coords.attrs.get("explained", [0, 0])
        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        cmap = plt.get_cmap("Set1")
        for i, (g, sub) in enumerate(df.groupby("group")):
            ax.scatter(sub["PCo1"], sub["PCo2"], label=str(g), color=cmap(i), s=60)
        ax.set_xlabel(f"PCo1 ({explained[0]*100:.1f}%)")
        ax.set_ylabel(f"PCo2 ({explained[1]*100:.1f}%)")
        ax.legend(title="Grup")
        return self._save(fig, title)

    def mag_quality(self, mag_df: pd.DataFrame,
                    title: str = "MAG kalitesi (tamlık - kontaminasyon)") -> int:
        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        colors = {"high": "#2CA02C", "medium": "#FF7F0E", "low": "#D62728"}
        for tier, sub in mag_df.groupby("quality"):
            ax.scatter(sub["completeness"], sub["contamination"],
                       label=tier, color=colors.get(tier, "#333333"), s=50)
        ax.axvline(90, ls="--", c="grey", lw=0.8)
        ax.axhline(5, ls="--", c="grey", lw=0.8)
        ax.set_xlabel("Tamlık (%)")
        ax.set_ylabel("Kontaminasyon (%)")
        ax.legend(title="Kalite")
        return self._save(fig, title)

    def volcano(self, diff_df: pd.DataFrame,
                title: str = "Diferansiyel bolluk (volcano)") -> int:
        fig, ax = plt.subplots(figsize=(5.5, 4.5))
        q = diff_df["q_value"].clip(lower=1e-10)
        y = -np.log10(q)
        sig = diff_df["significant"]
        ax.scatter(diff_df["effect_clr"][~sig], y[~sig], color="#BBBBBB", s=25, label="ns")
        ax.scatter(diff_df["effect_clr"][sig], y[sig], color="#D62728", s=35, label="q<0.05")
        ax.axhline(-np.log10(0.05), ls="--", c="grey", lw=0.8)
        ax.set_xlabel("Etki (CLR log-fark)")
        ax.set_ylabel("-log10(q)")
        ax.legend()
        return self._save(fig, title)

    def compile_pdf(self, pdf_path: str | Path) -> Path:
        """Tüm şekilleri tek çok-sayfalı PDF'te toplar (06_figures.pdf)."""
        pdf_path = Path(pdf_path)
        with PdfPages(pdf_path) as pdf:
            for entry in self.manifest.entries:
                img = plt.imread(entry.png_path)
                fig, ax = plt.subplots(figsize=(8.27, 5.5))  # ~A4 genişlik
                ax.imshow(img)
                ax.axis("off")
                fig.savefig(pdf, format="pdf", bbox_inches="tight")
                plt.close(fig)
        return pdf_path
