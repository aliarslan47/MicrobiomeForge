"""M2 birim testleri: ham veri QC istatistiği + numaralı çıktı kuralı."""

from pathlib import Path

from microbiomeforge.outputs import OUTPUT_FILES, output_path, write_table
from microbiomeforge.qc_stats import compute_raw_stats, qc_delta, qc_table


def _write_fastq(path: Path, seqs_quals):
    with open(path, "wt") as fh:
        for i, (s, q) in enumerate(seqs_quals):
            fh.write(f"@r{i}\n{s}\n+\n{q}\n")


def test_compute_raw_stats_basic(tmp_path):
    fq = tmp_path / "s.fastq"
    reads = [("ACGTACGTGG", "IIIIIIIIII"), ("ACGT", "IIII")]  # Q40 (I)
    _write_fastq(fq, reads)
    st = compute_raw_stats(fq, sample="s", category="short")
    assert st.n_reads == 2
    assert st.total_bases == 14
    assert st.min_length == 4 and st.max_length == 10
    assert st.mean_quality == 40.0
    # GC: ilk okuma 6 (C,G,C,G,G,G? -> "ACGTACGTGG": G,C,G,C,G,G =6) + ikinci "ACGT": G,C=2 => 8/14
    assert abs(st.gc_percent - round(100 * 8 / 14, 2)) < 0.01


def test_n50(tmp_path):
    fq = tmp_path / "s.fastq"
    reads = [("A" * L, "I" * L) for L in (100, 100, 100, 1000)]
    _write_fastq(fq, reads)
    st = compute_raw_stats(fq, "s", "long")
    assert st.n50 == 1000  # 1300 toplam, yarısı 650; en uzun 1000 tek başına aşar


def test_qc_delta_before_after(tmp_path):
    fq_raw = tmp_path / "raw.fastq"
    fq_trim = tmp_path / "trim.fastq"
    _write_fastq(fq_raw, [("A" * 100, "I" * 100)] * 10)
    _write_fastq(fq_trim, [("A" * 100, "I" * 100)] * 8)
    stats = [
        compute_raw_stats(fq_raw, "s1", "short", stage="raw"),
        compute_raw_stats(fq_trim, "s1", "short", stage="trimmed"),
    ]
    df = qc_table(stats)
    delta = qc_delta(df)
    assert delta.iloc[0]["reads_retained_pct"] == 80.0


def test_output_naming_rule():
    # Rapor en yüksek numara olmalı.
    keys = list(OUTPUT_FILES.values())
    assert keys[0].startswith("01_")
    assert keys[-1] == "07_report.pdf"
    # Numaralar sıralı ve artan.
    prefixes = [int(v.split("_")[0]) for v in OUTPUT_FILES.values()]
    assert prefixes == sorted(prefixes)


def test_write_table_flat(tmp_path):
    import pandas as pd
    df = pd.DataFrame({"sample": ["s1"], "n_reads": [10]})
    p = write_table(df, tmp_path, "raw_qc")
    assert p.name == "01_rawdata_qc_stats.tsv"
    assert p.parent == Path(tmp_path)  # düz: doğrudan outdir altında
    assert p.exists()
