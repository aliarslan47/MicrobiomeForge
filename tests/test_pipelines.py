"""M3 birim testleri: örnek sayfası + platform-farkında nf-core planlama."""

import gzip
from pathlib import Path

import pytest

from microbiomeforge.pipelines import (
    Sample,
    load_samplesheet,
    plan_pipelines,
    resolve_platforms,
)
from microbiomeforge.references import ToolRegistry


def _write_fastq(path, reads, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "wt") as fh:
        for i, (s, q) in enumerate(reads):
            fh.write(f"@r{i}\n{s}\n+\n{q}\n")


def test_load_samplesheet_example():
    samples = load_samplesheet(Path("config/samplesheet.example.csv"))
    assert len(samples) == 6
    assert {s.group for s in samples} == {"pre", "post", "environment"}
    assert any(s.long_reads for s in samples)


def test_missing_group_column_rejected(tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("sample,fastq_1\ns1,a.fastq\n")
    with pytest.raises(ValueError, match="eksik sütun"):
        load_samplesheet(p)


def test_duplicate_sample_rejected(tmp_path):
    p = tmp_path / "dup.csv"
    p.write_text("sample,group,fastq_1\ns1,a,x.fq\ns1,b,y.fq\n")
    with pytest.raises(ValueError, match="Yinelenen"):
        load_samplesheet(p)


def test_resolve_platform_override():
    samples = [
        Sample("s1", "pre", fastq_1="x.fq", platform="illumina"),
        Sample("s2", "post", long_reads="y.fq", platform="ont"),
    ]
    resolve_platforms(samples, detect=False)
    assert samples[0].category == "short"
    assert samples[1].category == "long"


def test_resolve_platform_auto_detects(tmp_path):
    short_fq = tmp_path / "short.fastq"
    long_fq = tmp_path / "long.fastq"
    _write_fastq(short_fq, [("A" * 150, "I" * 150)] * 200)
    _write_fastq(long_fq, [("A" * 6000, "5" * 6000)] * 200)  # '5'=Q20-ish low
    samples = [
        Sample("s1", "pre", fastq_1=str(short_fq)),
        Sample("s2", "post", long_reads=str(long_fq)),
    ]
    resolve_platforms(samples)
    assert samples[0].category == "short"
    assert samples[1].category == "long"


def test_plan_routes_and_registers_tools():
    samples = [
        Sample("s1", "pre", fastq_1="x.fq", platform="illumina", category="short"),
        Sample("s2", "post", long_reads="y.fq", platform="ont", category="long"),
    ]
    reg = ToolRegistry()
    plans = plan_pipelines(samples, reg)
    pipelines = {p.pipeline for p in plans}
    assert "nf-core/taxprofiler" in pipelines
    assert "nf-core/mag" in pipelines
    assert "nf-core/funcscan" in pipelines
    # kısa yol metaSPAdes, uzun yol metaMDBG kullanmalı
    assemblers = {p.params.get("assembler") for p in plans if p.pipeline == "nf-core/mag"}
    assert {"spades", "metamdbg"} <= assemblers
    # kullanılan araçlar kaynakçaya işlendi
    keys = {r.key for r in reg.bibliography()}
    assert {"kraken2", "sylph", "semibin2", "metamdbg"} <= keys


def test_command_uses_podman_profile():
    samples = [Sample("s1", "pre", fastq_1="x.fq", platform="illumina", category="short")]
    reg = ToolRegistry()
    plans = plan_pipelines(samples, reg)
    cmd = plans[0].command(outdir="results", container="podman")
    assert "nextflow" in cmd and "-profile" in cmd and "podman" in cmd
