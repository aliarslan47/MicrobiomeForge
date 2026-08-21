"""Platform otomatik tespiti (kısa vs uzun okuma) — ham FASTQ'tan.

Sistem, örneğin hangi sekanslama teknolojisiyle üretildiğini ham okumadan anlar ve
downstream'i (kısa→Illumina yolu, uzun→ONT/PacBio yolu) buna göre yönlendirir.

Sezgi (platform karakteristikleri):
  * Illumina (kısa): okumalar ~50–300 bp, dar uzunluk dağılımı, yüksek kalite (Q30+).
  * ONT (uzun):      okumalar çok değişken, genelde medyan >1 kb, orta kalite (~Q10–20).
  * PacBio HiFi:     okumalar ~10–20 kb, yüksek kalite (Q20+).

Karar yalnız ham okumanın uzunluk dağılımı + ortalama Phred kalitesine dayanır;
harici bir metadata gerektirmez.
"""

from __future__ import annotations

import gzip
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator, Literal

Category = Literal["short", "long"]
Platform = Literal["illumina", "ont", "pacbio_hifi", "unknown"]

# Eşikler (platform karakteristiklerinden; değiştirilebilir sabitler).
SHORT_MAX_MEDIAN_LEN = 400      # bu medyanın altı + az uzun-okuma → kısa
LONG_FRAC_OVER_1KB = 0.05       # okumaların >%5'i 1 kb üstündeyse uzun eğilimi
HIFI_MIN_MEAN_Q = 20.0          # uzun okumada bu ve üstü ortalama Q → PacBio HiFi


@dataclass(frozen=True)
class PlatformCall:
    category: Category
    platform: Platform
    median_length: float
    p90_length: float
    frac_over_1kb: float
    mean_quality: float
    n_reads_sampled: int
    confidence: float  # 0..1

    def as_row(self) -> dict:
        return {
            "category": self.category,
            "platform": self.platform,
            "median_length": round(self.median_length, 1),
            "p90_length": round(self.p90_length, 1),
            "frac_over_1kb": round(self.frac_over_1kb, 4),
            "mean_quality": round(self.mean_quality, 2),
            "n_reads_sampled": self.n_reads_sampled,
            "confidence": round(self.confidence, 3),
        }


def _open(path: Path):
    """Düz veya gzip FASTQ'yu şeffaf açar."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")


def _iter_fastq(path: Path, max_reads: int) -> Iterator[tuple[int, float]]:
    """(okuma_uzunluğu, ortalama_Phred_kalite) çiftlerini üretir. Phred+33 varsayar."""
    with _open(path) as fh:
        count = 0
        while count < max_reads:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline().rstrip("\n")
            fh.readline()  # '+'
            qual = fh.readline().rstrip("\n")
            if not qual:
                break
            length = len(seq)
            if qual:
                mean_q = sum(ord(c) - 33 for c in qual) / len(qual)
            else:
                mean_q = 0.0
            yield length, mean_q
            count += 1


def _percentile(values: list[int], pct: float) -> float:
    if not values:
        return 0.0
    s = sorted(values)
    k = (len(s) - 1) * pct
    lo = int(k)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def detect_platform(fastq_path: str | Path, max_reads: int = 5000) -> PlatformCall:
    """Bir FASTQ dosyasından platform kategorisini ve teknolojisini kestirir."""
    lengths: list[int] = []
    quals: list[float] = []
    for length, mean_q in _iter_fastq(Path(fastq_path), max_reads):
        lengths.append(length)
        quals.append(mean_q)

    if not lengths:
        return PlatformCall("short", "unknown", 0, 0, 0, 0, 0, confidence=0.0)

    median_len = statistics.median(lengths)
    p90 = _percentile(lengths, 0.90)
    frac_over_1kb = sum(1 for x in lengths if x > 1000) / len(lengths)
    mean_q = statistics.mean(quals)

    is_short = median_len <= SHORT_MAX_MEDIAN_LEN and frac_over_1kb < LONG_FRAC_OVER_1KB

    if is_short:
        category: Category = "short"
        platform: Platform = "illumina"
        # Medyan ne kadar kısa ve dağılım ne kadar darsa güven o kadar yüksek.
        confidence = min(1.0, 0.6 + (SHORT_MAX_MEDIAN_LEN - median_len) / SHORT_MAX_MEDIAN_LEN * 0.4)
    else:
        category = "long"
        if mean_q >= HIFI_MIN_MEAN_Q:
            platform = "pacbio_hifi"
        else:
            platform = "ont"
        confidence = min(1.0, 0.5 + frac_over_1kb * 0.5)

    return PlatformCall(
        category=category,
        platform=platform,
        median_length=median_len,
        p90_length=p90,
        frac_over_1kb=frac_over_1kb,
        mean_quality=mean_q,
        n_reads_sampled=len(lengths),
        confidence=confidence,
    )


def detect_sample(
    fastq_1: str | Path,
    fastq_2: str | Path | None = None,
    long_reads: str | Path | None = None,
    override: Platform | None = None,
    max_reads: int = 5000,
) -> PlatformCall:
    """Örnek-düzeyi tespit.

    - `long_reads` verilmişse onu; yoksa `fastq_1`'i inceler.
    - `override` verilmişse tespit atlanır (kullanıcı platformu zorlayabilir).
    """
    if override is not None:
        cat: Category = "short" if override == "illumina" else "long"
        return PlatformCall(cat, override, 0, 0, 0, 0, 0, confidence=1.0)
    target = long_reads if long_reads is not None else fastq_1
    return detect_platform(target, max_reads=max_reads)
