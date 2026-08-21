"""Çıktı dosyalama kuralı (değişmez).

- Dosyalar 1'den numaralı, sıralı; RAPOR en son (en yüksek numara).
- Düz yapı: çıktı dizininde klasör-içinde-klasör YOK.
- Yalnız SONUÇLAR: ara/geçici/log dosyaları çıktı dizininde tutulmaz (ayrı work-dir).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# Kanonik numaralı çıktı adları. Rapor daima en son sırada.
OUTPUT_FILES: dict[str, str] = {
    "raw_qc": "01_rawdata_qc_stats.tsv",
    "taxonomic_profile": "02_taxonomic_profile.tsv",
    "diversity_stats": "03_diversity_stats.tsv",
    "mag_quality": "04_mag_quality.tsv",
    "amr_bgc": "05_amr_bgc_results.tsv",
    "figures": "06_figures.pdf",
    "report": "07_report.pdf",
}


def output_path(outdir: str | Path, key: str) -> Path:
    if key not in OUTPUT_FILES:
        raise KeyError(f"Bilinmeyen çıktı anahtarı: {key!r}. Geçerli: {list(OUTPUT_FILES)}")
    return Path(outdir) / OUTPUT_FILES[key]


def write_table(df: pd.DataFrame, outdir: str | Path, key: str) -> Path:
    """Bir sonuç tablosunu kanonik numaralı ada TSV olarak yazar (düz dizin)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = output_path(outdir, key)
    df.to_csv(path, sep="\t", index=False)
    return path
