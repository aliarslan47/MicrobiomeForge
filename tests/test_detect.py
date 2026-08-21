"""M1 birim testleri: platform otomatik tespiti (sentetik FASTQ ile)."""

import gzip
import random
from pathlib import Path

from microbiomeforge.detect import detect_platform, detect_sample


def _write_fastq(path: Path, reads: list[tuple[str, str]], gz: bool = False):
    opener = gzip.open if gz else open
    with opener(path, "wt") as fh:
        for i, (seq, qual) in enumerate(reads):
            fh.write(f"@read{i}\n{seq}\n+\n{qual}\n")


def _illumina_reads(n=500, length=150, q=38):
    reads = []
    for _ in range(n):
        seq = "".join(random.choice("ACGT") for _ in range(length))
        qual = chr(q + 33) * length
        reads.append((seq, qual))
    return reads


def _long_reads(n=300, mean_len=8000, q=12):
    reads = []
    for _ in range(n):
        length = max(500, int(random.gauss(mean_len, mean_len * 0.3)))
        seq = "".join(random.choice("ACGT") for _ in range(length))
        qual = chr(q + 33) * length
        reads.append((seq, qual))
    return reads


def test_detects_illumina_short(tmp_path):
    fq = tmp_path / "sample_R1.fastq"
    _write_fastq(fq, _illumina_reads())
    call = detect_platform(fq)
    assert call.category == "short"
    assert call.platform == "illumina"
    assert call.median_length <= 400
    assert call.confidence > 0.6


def test_detects_ont_long(tmp_path):
    fq = tmp_path / "sample_ont.fastq"
    _write_fastq(fq, _long_reads(q=12))
    call = detect_platform(fq)
    assert call.category == "long"
    assert call.platform == "ont"
    assert call.frac_over_1kb > 0.5


def test_detects_pacbio_hifi_by_quality(tmp_path):
    fq = tmp_path / "sample_hifi.fastq"
    _write_fastq(fq, _long_reads(mean_len=12000, q=30))
    call = detect_platform(fq)
    assert call.category == "long"
    assert call.platform == "pacbio_hifi"


def test_gzip_supported(tmp_path):
    fq = tmp_path / "sample_R1.fastq.gz"
    _write_fastq(fq, _illumina_reads(), gz=True)
    call = detect_platform(fq)
    assert call.category == "short"


def test_override_skips_detection(tmp_path):
    call = detect_sample(fastq_1="nonexistent.fastq", override="ont")
    assert call.category == "long"
    assert call.platform == "ont"
    assert call.confidence == 1.0


def test_long_reads_field_takes_priority(tmp_path):
    short_fq = tmp_path / "short.fastq"
    long_fq = tmp_path / "long.fastq"
    _write_fastq(short_fq, _illumina_reads())
    _write_fastq(long_fq, _long_reads())
    call = detect_sample(fastq_1=short_fq, long_reads=long_fq)
    assert call.category == "long"
