"""M8 entegrasyon testi: sentetik veriyle uçtan uca koşu (CLI dahil).

nf-core fiilen çalıştırılmaz; onun çıktıları (Bracken profilleri, CheckM2, AMR) sentetik
üretilir. Ham FASTQ'lar gerçekten okunup QC istatistiği hesaplanır (zincirin başı).
"""

import gzip
import random
from pathlib import Path

import pandas as pd

from microbiomeforge.cli import main
from microbiomeforge.outputs import OUTPUT_FILES


def _write_fastq(path: Path, n, length, q, gz=False):
    opener = gzip.open if gz else open
    with opener(path, "wt") as fh:
        for i in range(n):
            seq = "".join(random.choice("ACGT") for _ in range(length))
            fh.write(f"@r{i}\n{seq}\n+\n{chr(q+33)*length}\n")


def _bracken(path: Path, fractions: dict):
    pd.DataFrame({
        "name": list(fractions),
        "fraction_total_reads": list(fractions.values()),
    }).to_csv(path, sep="\t", index=False)


def _setup(tmp_path):
    data = tmp_path / "data"; data.mkdir()
    prof = tmp_path / "profiles"; prof.mkdir()
    # 4 örnek: pre (taxonA baskın) vs post (taxonB baskın), kısa okuma.
    samples = {
        "s_pre1": ("pre", {"TaxA": 80, "TaxB": 10, "TaxC": 10}),
        "s_pre2": ("pre", {"TaxA": 75, "TaxB": 15, "TaxC": 10}),
        "s_post1": ("post", {"TaxA": 10, "TaxB": 80, "TaxC": 10}),
        "s_post2": ("post", {"TaxA": 15, "TaxB": 75, "TaxC": 10}),
    }
    rows = []
    for name, (group, fracs) in samples.items():
        fq = data / f"{name}_R1.fastq.gz"
        _write_fastq(fq, n=300, length=150, q=36, gz=True)
        _bracken(prof / f"{name}.bracken", fracs)
        rows.append({"sample": name, "fastq_1": str(fq), "fastq_2": "",
                     "long_reads": "", "group": group,
                     "sample_type": "microbiome", "platform": "auto"})
    ss = tmp_path / "samplesheet.csv"
    pd.DataFrame(rows).to_csv(ss, index=False)

    # CheckM2 + AMR (sentetik nf-core mag/funcscan çıktıları)
    checkm2 = tmp_path / "quality_report.tsv"
    pd.DataFrame({"Name": ["bin.1", "bin.2"], "Completeness": [96.0, 62.0],
                  "Contamination": [1.5, 7.0]}).to_csv(checkm2, sep="\t", index=False)
    amr = tmp_path / "amr.tsv"
    pd.DataFrame({"input_file_name": ["s_pre1", "s_post1"],
                  "gene_symbol": ["blaTEM", "tetM"]}).to_csv(amr, sep="\t", index=False)
    return ss, prof, checkm2, amr


def test_end_to_end_cli(tmp_path):
    ss, prof, checkm2, amr = _setup(tmp_path)
    outdir = tmp_path / "results"
    rc = main([
        "run",
        "--samplesheet", str(ss),
        "--profiles-dir", str(prof),
        "--outdir", str(outdir),
        "--workdir", str(tmp_path / "work"),
        "--checkm2", str(checkm2),
        "--amr", str(amr),
        "--project", "IntegrationTest",
    ])
    assert rc == 0

    # Numaralı çıktılar 01..06 + rapor mevcut.
    for key in ["raw_qc", "taxonomic_profile", "diversity_stats", "mag_quality", "amr_bgc", "figures"]:
        assert (outdir / OUTPUT_FILES[key]).exists(), f"eksik: {OUTPUT_FILES[key]}"
    # Rapor: PDF (xelatex varsa) ya da MD fallback — her hâlde 07_ ile.
    report_files = list(outdir.glob("07_report.*"))
    assert report_files, "rapor üretilmedi"

    # Düz yapı: çıktı dizininde alt-klasör yok.
    assert not any(p.is_dir() for p in outdir.iterdir())

    # Ham QC istatistiği gerçekten hesaplandı (zincir başı).
    qc = pd.read_csv(outdir / OUTPUT_FILES["raw_qc"], sep="\t")
    assert set(qc["sample"]) == {"s_pre1", "s_pre2", "s_post1", "s_post2"}
    assert (qc["n_reads"] == 300).all()
    assert (qc["category"] == "short").all()  # otomatik platform tespiti

    # Profil matrisi taxon×örnek.
    prof_df = pd.read_csv(outdir / OUTPUT_FILES["taxonomic_profile"], sep="\t")
    assert set(prof_df["taxon"]) == {"TaxA", "TaxB", "TaxC"}


def test_run_analysis_reports_comparative_design_from_metadata(tmp_path):
    """İki grup → run_analysis karşılaştırmalı tasarımı bildirmeli (PERMANOVA koşar)."""
    from microbiomeforge.orchestrate import run_analysis
    from microbiomeforge.pipelines import load_samplesheet, resolve_platforms

    ss, prof, checkm2, amr = _setup(tmp_path)
    samples = resolve_platforms(load_samplesheet(ss))
    profiles = {s.sample: str(prof / f"{s.sample}.bracken") for s in samples}
    res = run_analysis(samples, profiles, tmp_path / "out", tmp_path / "wk")
    assert res["design"] == "comparative"
    assert res["permanova_p"] is not None


def test_single_design_skips_comparative_stats(tmp_path):
    """Sisteme 'tekli' denince karşılaştırmalı testler atlanmalı ve not verilmeli."""
    from microbiomeforge.orchestrate import run_analysis
    from microbiomeforge.pipelines import load_samplesheet, resolve_platforms

    ss, prof, checkm2, amr = _setup(tmp_path)
    samples = resolve_platforms(load_samplesheet(ss))
    profiles = {s.sample: str(prof / f"{s.sample}.bracken") for s in samples}
    res = run_analysis(samples, profiles, tmp_path / "out", tmp_path / "wk", design="single")
    assert res["design"] == "single"
    assert res["permanova_p"] is None          # karşılaştırma atlandı
    assert res["n_significant_taxa"] == 0
    assert "atlan" in res["design_note"].lower()


def test_cli_declared_comparative_on_single_group_fails(tmp_path):
    """CLI --design comparative ama tek grup → gürültülü hata (rc!=0)."""
    ss, prof, checkm2, amr = _setup(tmp_path)
    # Tüm örnekleri tek gruba indir.
    df = pd.read_csv(ss); df["group"] = "cohort"; df.to_csv(ss, index=False)
    rc = main([
        "run", "--samplesheet", str(ss), "--profiles-dir", str(prof),
        "--outdir", str(tmp_path / "results"), "--workdir", str(tmp_path / "work"),
        "--design", "comparative",
    ])
    assert rc != 0
