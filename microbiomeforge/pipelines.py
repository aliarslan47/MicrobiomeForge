"""Örnek sayfası şeması + nf-core koşucu sarmalayıcıları (platform-farkında yönlendirme).

Bu modül örnekleri okur, platformu (gerekirse otomatik) çözer ve her örnek grubunu
uygun nf-core pipeline'ına yönlendirir:

  * kısa (Illumina) → taxprofiler (Kraken2+Bracken+sylph) + mag (metaSPAdes+SemiBin2)
  * uzun (ONT/PacBio) → taxprofiler (uzun profil: sylph/Kraken2) + mag (metaFlye/metaMDBG+SemiBin2)

Komutlar `NextflowPlan` olarak İNŞA edilir; fiili çalıştırma CLI/çalışma zamanında
yapılır (bu modül komutu üretmekle sorumludur, test edilebilir kalır). Bir pipeline
planlandığında kullandığı araçlar `ToolRegistry`'ye işlenir (kaynakça için).
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from .detect import Platform, detect_sample
from .references import ToolRegistry

REQUIRED_COLUMNS = {"sample", "group"}
VALID_SAMPLE_TYPES = {"microbiome", "microbiota", "environment"}


@dataclass
class Sample:
    sample: str
    group: str
    fastq_1: Optional[str] = None
    fastq_2: Optional[str] = None
    long_reads: Optional[str] = None
    sample_type: str = "microbiome"
    platform: str = "auto"          # auto | illumina | ont | pacbio_hifi
    category: Optional[str] = None  # çözüldükten sonra: short | long

    def has_reads(self) -> bool:
        return bool(self.fastq_1 or self.long_reads)


def load_samplesheet(path: str | Path) -> list[Sample]:
    """CSV örnek sayfasını okur ve doğrular."""
    path = Path(path)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        cols = set(reader.fieldnames or [])
        missing = REQUIRED_COLUMNS - cols
        if missing:
            raise ValueError(f"Örnek sayfasında eksik sütun(lar): {sorted(missing)}")
        samples: list[Sample] = []
        seen = set()
        for row in reader:
            name = (row.get("sample") or "").strip()
            if not name:
                raise ValueError("Boş 'sample' adı var.")
            if name in seen:
                raise ValueError(f"Yinelenen örnek adı: {name}")
            seen.add(name)
            stype = (row.get("sample_type") or "microbiome").strip() or "microbiome"
            if stype not in VALID_SAMPLE_TYPES:
                raise ValueError(f"{name}: geçersiz sample_type {stype!r} (geçerli: {sorted(VALID_SAMPLE_TYPES)})")
            s = Sample(
                sample=name,
                group=(row.get("group") or "").strip(),
                fastq_1=(row.get("fastq_1") or "").strip() or None,
                fastq_2=(row.get("fastq_2") or "").strip() or None,
                long_reads=(row.get("long_reads") or "").strip() or None,
                sample_type=stype,
                platform=(row.get("platform") or "auto").strip() or "auto",
            )
            if not s.group:
                raise ValueError(f"{name}: 'group' boş olamaz (karşılaştırmalı tasarım gerekli).")
            if not s.has_reads():
                raise ValueError(f"{name}: en az fastq_1 veya long_reads verilmeli.")
            samples.append(s)
    if len({s.group for s in samples}) < 2:
        # Karşılaştırmalı tasarım çekirdekte; tek grup uyarısı.
        pass  # tek gruplu betimsel koşuya izin ver, ama istatistik modülü uyaracak.
    return samples


def resolve_platforms(samples: list[Sample], detect: bool = True, max_reads: int = 5000) -> list[Sample]:
    """platform='auto' olan örnekleri ham okumadan tespit ederek doldurur."""
    for s in samples:
        if s.platform != "auto":
            s.category = "short" if s.platform == "illumina" else "long"
            continue
        if not detect:
            continue
        call = detect_sample(
            fastq_1=s.fastq_1, fastq_2=s.fastq_2, long_reads=s.long_reads, max_reads=max_reads
        )
        s.platform = call.platform
        s.category = call.category
    return samples


@dataclass
class NextflowPlan:
    pipeline: str                 # "nf-core/taxprofiler" | "nf-core/mag" | "nf-core/funcscan"
    input_samples: list[str]
    category: str                 # short | long | hybrid
    params: dict = field(default_factory=dict)
    tools: list[str] = field(default_factory=list)  # references.CATALOG anahtarları

    def command(self, outdir: str, container: str = "podman") -> list[str]:
        cmd = [
            "nextflow", "run", self.pipeline,
            "-profile", container,
            "--outdir", outdir,
        ]
        for k, v in self.params.items():
            cmd += [f"--{k}", str(v)]
        return cmd


def plan_pipelines(samples: list[Sample], registry: ToolRegistry) -> list[NextflowPlan]:
    """Örnekleri kategoriye göre gruplayıp nf-core planları üretir + araçları işler."""
    shorts = [s for s in samples if s.category == "short"]
    longs = [s for s in samples if s.category == "long"]
    plans: list[NextflowPlan] = []

    if shorts:
        plans.append(NextflowPlan(
            pipeline="nf-core/taxprofiler",
            input_samples=[s.sample for s in shorts],
            category="short",
            params={"run_kraken2": "true", "run_bracken": "true", "run_sylph": "true"},
            tools=["nfcore_taxprofiler", "kraken2", "bracken", "sylph", "fastp"],
        ))
        plans.append(NextflowPlan(
            pipeline="nf-core/mag",
            input_samples=[s.sample for s in shorts],
            category="short",
            params={"assembler": "spades", "binner": "semibin2"},
            tools=["nfcore_mag", "semibin2", "checkm2"],
        ))

    if longs:
        plans.append(NextflowPlan(
            pipeline="nf-core/taxprofiler",
            input_samples=[s.sample for s in longs],
            category="long",
            params={"run_sylph": "true", "run_kraken2": "true"},
            tools=["nfcore_taxprofiler", "sylph", "kraken2", "nanoplot"],
        ))
        plans.append(NextflowPlan(
            pipeline="nf-core/mag",
            input_samples=[s.sample for s in longs],
            category="long",
            params={"assembler": "metamdbg", "binner": "semibin2"},
            tools=["nfcore_mag", "metamdbg", "semibin2", "checkm2"],
        ))

    # AMR/BGC her kategori için (funcscan) — MAG'lar üzerinden.
    if shorts or longs:
        plans.append(NextflowPlan(
            pipeline="nf-core/funcscan",
            input_samples=[s.sample for s in samples],
            category="hybrid",
            params={"run_amp_screening": "true", "run_arg_screening": "true"},
            tools=[],  # funcscan alt-araç DOI'leri kullanıldıklarında ayrıca işlenecek
        ))

    # Planlanan tüm araçları kaynakçaya işle.
    for p in plans:
        for t in p.tools:
            registry.mark_used(t)
    return plans
