"""Ham veri QC istatistiği — istatistik zincirinin BAŞI.

Rapor kuralı: tüm istatistik ham okumadan başlar. Bu modül her örnek için okuma
sayısı, taban sayısı, uzunluk dağılımı (min/medyan/ortalama/N50/maks), ortalama
Phred kalitesi ve GC yüzdesini FASTQ'tan doğrudan hesaplar ve trim öncesi→sonrası
karşılaştırma tablosunu üretir (`01_rawdata_qc_stats.tsv`).

Hesaplama MicrobiomeForge içinde yapılır (harici araca bağımlı değil). Üretimde
nf-core fastp/NanoPlot raporları mevcutsa `parsers` modülü onları ayrıca içe aktarır
ve o araçları kaynakçaya (references) işler; burada sahte atıf yapılmaz.
"""

from __future__ import annotations

import gzip
import statistics
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class RawReadStats:
    sample: str
    stage: str          # "raw" | "trimmed"
    category: str       # "short" | "long"
    n_reads: int
    total_bases: int
    min_length: int
    median_length: float
    mean_length: float
    max_length: int
    n50: int
    mean_quality: float
    gc_percent: float


def _open(path: Path):
    path = Path(path)
    return gzip.open(path, "rt") if path.suffix == ".gz" else open(path, "rt")


def _n50(lengths: list[int]) -> int:
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    half = sum(s) / 2
    acc = 0
    for x in s:
        acc += x
        if acc >= half:
            return x
    return s[-1]


def compute_raw_stats(
    fastq_path: str | Path,
    sample: str,
    category: str,
    stage: str = "raw",
    max_reads: int | None = None,
) -> RawReadStats:
    """Bir FASTQ'tan ham okuma istatistiğini hesaplar (Phred+33)."""
    lengths: list[int] = []
    total_q = 0.0
    total_bases = 0
    gc = 0
    with _open(fastq_path) as fh:
        n = 0
        while True:
            if max_reads is not None and n >= max_reads:
                break
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()
            qual = fh.readline().rstrip("\n")
            if not qual:
                break
            L = len(seq)
            lengths.append(L)
            total_bases += L
            gc += sum(1 for c in seq if c in "GCgc")
            if qual:
                total_q += sum(ord(c) - 33 for c in qual)
            n += 1

    if not lengths:
        return RawReadStats(sample, stage, category, 0, 0, 0, 0, 0, 0, 0, 0.0, 0.0)

    return RawReadStats(
        sample=sample,
        stage=stage,
        category=category,
        n_reads=len(lengths),
        total_bases=total_bases,
        min_length=min(lengths),
        median_length=round(statistics.median(lengths), 1),
        mean_length=round(total_bases / len(lengths), 1),
        max_length=max(lengths),
        n50=_n50(lengths),
        mean_quality=round(total_q / total_bases, 2) if total_bases else 0.0,
        gc_percent=round(100 * gc / total_bases, 2) if total_bases else 0.0,
    )


def qc_table(stats: list[RawReadStats]) -> pd.DataFrame:
    """RawReadStats listesini tabloya çevirir (örnek × aşama)."""
    df = pd.DataFrame(asdict(s) for s in stats)
    if not df.empty:
        df = df.sort_values(["sample", "stage"]).reset_index(drop=True)
    return df


def qc_delta(df: pd.DataFrame) -> pd.DataFrame:
    """Örnek başına trim öncesi→sonrası değişim (okuma/taban tutulma oranı)."""
    rows = []
    for sample, g in df.groupby("sample"):
        raw = g[g["stage"] == "raw"]
        trimmed = g[g["stage"] == "trimmed"]
        if raw.empty or trimmed.empty:
            continue
        r, t = raw.iloc[0], trimmed.iloc[0]
        rows.append({
            "sample": sample,
            "reads_retained_pct": round(100 * t["n_reads"] / r["n_reads"], 2) if r["n_reads"] else 0.0,
            "bases_retained_pct": round(100 * t["total_bases"] / r["total_bases"], 2) if r["total_bases"] else 0.0,
            "mean_q_raw": r["mean_quality"],
            "mean_q_trimmed": t["mean_quality"],
        })
    return pd.DataFrame(rows)
