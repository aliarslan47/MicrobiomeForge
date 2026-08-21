"""nf-core çıktı ayrıştırıcıları → normalize tablolar.

nf-core pipeline'larının ürettiği standart dosyaları MicrobiomeForge'un downstream
istatistik/rapor katmanının beklediği normalize tablolara çevirir:

  * taxprofiler (Bracken / sylph / MetaPhlAn) → taxon × örnek bağıl-bolluk matrisi
    → 02_taxonomic_profile.tsv
  * CheckM2 quality_report → MAG kalite tablosu → 04_mag_quality.tsv
  * funcscan (hAMRonization/ARG) → AMR/BGC tablosu → 05_amr_bgc_results.tsv

Ayrıştırıcılar toleranslıdır: eksik sütunlara makul varsayılanlar uygular, ama
sütun adlarını uydurmaz — tanınmayan biçimde açık hata verir.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def parse_bracken(path: str | Path, sample: str) -> pd.DataFrame:
    """Bir Bracken çıktı tablosunu (name, fraction_total_reads) okur.

    Dönüş: [taxon, abundance] (bağıl bolluk, örnek içinde toplamı ~1).
    """
    df = pd.read_csv(path, sep="\t")
    name_col = _first_present(df, ["name", "taxon", "clade_name"])
    frac_col = _first_present(df, ["fraction_total_reads", "new_est_reads", "abundance", "relative_abundance"])
    if name_col is None or frac_col is None:
        raise ValueError(f"{path}: Bracken benzeri sütunlar bulunamadı (name/fraction).")
    out = df[[name_col, frac_col]].copy()
    out.columns = ["taxon", "abundance"]
    total = out["abundance"].sum()
    if total > 0:
        out["abundance"] = out["abundance"] / total  # bağıl bolluğa normalize et
    out["sample"] = sample
    return out.reset_index(drop=True)


def build_abundance_matrix(per_sample: dict[str, pd.DataFrame]) -> pd.DataFrame:
    """{örnek: [taxon, abundance]} → taxon (satır) × örnek (sütun) matrisi.

    Eksik taxon'lar 0 ile doldurulur.
    """
    frames = []
    for sample, df in per_sample.items():
        s = df.set_index("taxon")["abundance"].rename(sample)
        frames.append(s)
    if not frames:
        return pd.DataFrame()
    matrix = pd.concat(frames, axis=1).fillna(0.0)
    matrix.index.name = "taxon"
    return matrix.reset_index()


def parse_checkm2(path: str | Path) -> pd.DataFrame:
    """CheckM2 quality_report.tsv → [bin, completeness, contamination, quality]."""
    df = pd.read_csv(path, sep="\t")
    name_col = _first_present(df, ["Name", "bin", "Bin Id", "genome"])
    comp_col = _first_present(df, ["Completeness", "completeness"])
    cont_col = _first_present(df, ["Contamination", "contamination"])
    if not all([name_col, comp_col, cont_col]):
        raise ValueError(f"{path}: CheckM2 sütunları bulunamadı (Name/Completeness/Contamination).")
    out = df[[name_col, comp_col, cont_col]].copy()
    out.columns = ["bin", "completeness", "contamination"]
    out["quality"] = out.apply(_mag_quality_tier, axis=1)
    return out.sort_values("completeness", ascending=False).reset_index(drop=True)


def _mag_quality_tier(row) -> str:
    """MIMAG benzeri kalite katmanı (yüksek/orta/düşük)."""
    comp, cont = row["completeness"], row["contamination"]
    if comp >= 90 and cont <= 5:
        return "high"
    if comp >= 50 and cont <= 10:
        return "medium"
    return "low"


def parse_amr(path: str | Path) -> pd.DataFrame:
    """funcscan/hAMRonization benzeri AMR tablosu → normalize sütunlar."""
    df = pd.read_csv(path, sep="\t")
    gene_col = _first_present(df, ["gene_symbol", "gene_name", "Best_Hit_ARO", "gene"])
    sample_col = _first_present(df, ["input_file_name", "sample", "Contig"])
    if gene_col is None:
        raise ValueError(f"{path}: AMR gen sütunu bulunamadı.")
    cols = {gene_col: "gene"}
    if sample_col:
        cols[sample_col] = "sample"
    out = df.rename(columns=cols)
    keep = [c for c in ["sample", "gene"] if c in out.columns]
    return out[keep].reset_index(drop=True)


def _first_present(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None
