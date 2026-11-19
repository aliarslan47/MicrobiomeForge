"""Yerel taksonomik sınıflandırma motoru — FASTQ → Kraken2 → Bracken.

nf-core/Nextflow/konteyner olmadan çalışır: Kraken2 ve Bracken'ı doğrudan çağırıp
her örnek için `<örnek>.bracken` profil dosyası üretir. Bu dosyalar nf-core/taxprofiler
Bracken çıktısıyla aynı biçimdedir; rapor katmanı ikisini ayırt etmeden okur.

Bellek ihtiyacı çıkan tür sayısına değil Kraken2 veritabanının boyutuna bağlıdır.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from .pipelines import Sample


@dataclass
class ClassifyResult:
    sample: str
    bracken: Path
    kraken_report: Path
    read_len: int
    classified_pct: Optional[float]


def tool_version(cmd: str) -> Optional[str]:
    """`kraken2 --version` / `bracken -v` çıktısından sürüm numarasını okur."""
    flag = "-v" if cmd == "bracken" else "--version"
    try:
        out = subprocess.run([cmd, flag], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return None
    m = re.search(r"(\d+\.\d+(?:\.\d+)?)", out.stdout + out.stderr)
    return m.group(1) if m else None


def check_tools() -> None:
    missing = [t for t in ("kraken2", "bracken") if shutil.which(t) is None]
    if missing:
        raise RuntimeError(f"Gerekli araç(lar) PATH'te yok: {missing}")


def available_read_lengths(db: Path) -> list[int]:
    """Veritabanındaki Bracken k-mer dağılımı dosyalarından okuma uzunluklarını bulur."""
    lens = []
    for p in db.glob("database*mers.kmer_distrib"):
        m = re.match(r"database(\d+)mers\.kmer_distrib", p.name)
        if m:
            lens.append(int(m.group(1)))
    return sorted(lens)


def pick_read_len(mean_len: float, available: list[int]) -> int:
    """Ortalama okuma uzunluğuna en uygun Bracken dağılımı: altında kalan en büyük.

    Okumalar tüm dağılımlardan kısaysa en küçüğü; uzun okumada (ONT) en büyüğü seçilir.
    """
    if not available:
        raise ValueError("Veritabanında Bracken k-mer dağılımı (database*mers.kmer_distrib) yok.")
    below = [L for L in available if L <= mean_len]
    return max(below) if below else min(available)


def _classified_pct(report: Path) -> Optional[float]:
    """Kraken2 raporundan sınıflanan okuma yüzdesi (100 - unclassified)."""
    try:
        for line in report.read_text().splitlines():
            parts = line.split("\t")
            if len(parts) >= 6 and parts[3].strip() == "U":
                return round(100.0 - float(parts[0]), 2)
    except OSError:
        return None
    return 100.0


def classify_sample(
    s: Sample,
    db: Path,
    outdir: Path,
    mean_len: float,
    threads: int = 8,
    level: str = "S",
    bracken_threshold: int = 10,
    memory_mapping: bool = False,
) -> ClassifyResult:
    """Tek örnek: Kraken2 (rapor) → Bracken (bolluk tahmini)."""
    outdir.mkdir(parents=True, exist_ok=True)
    report = outdir / f"{s.sample}.k2report"
    bracken = outdir / f"{s.sample}.bracken"
    done = outdir / f"{s.sample}.done"
    read_len = pick_read_len(mean_len, available_read_lengths(db))
    if done.exists() and report.exists() and bracken.exists():
        # Devam: bu örnek aynı çalışma dizininde daha önce tamamlanmış.
        print(f"[MicrobiomeForge] {s.sample}: önceki Kraken2/Bracken sonucu kullanılıyor", file=sys.stderr)
        return ClassifyResult(s.sample, bracken, report, read_len, _classified_pct(report))

    cmd = ["kraken2", "--db", str(db), "--threads", str(threads),
           "--report", str(report), "--output", "-"]
    if memory_mapping:
        cmd.append("--memory-mapping")
    if s.long_reads:
        cmd.append(s.long_reads)
    elif s.fastq_2:
        cmd += ["--paired", s.fastq_1, s.fastq_2]
    else:
        cmd.append(s.fastq_1)
    print(f"[MicrobiomeForge] Kraken2: {s.sample}", file=sys.stderr)
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)

    bcmd = ["bracken", "-d", str(db), "-i", str(report), "-o", str(bracken),
            "-w", str(outdir / f"{s.sample}.bracken_report"),
            "-r", str(read_len), "-l", level, "-t", str(bracken_threshold)]
    print(f"[MicrobiomeForge] Bracken: {s.sample} (r={read_len}, l={level}, t={bracken_threshold})",
          file=sys.stderr)
    subprocess.run(bcmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    done.write_text(f"db={db}\nread_len={read_len}\n")
    return ClassifyResult(s.sample, bracken, report, read_len, _classified_pct(report))
