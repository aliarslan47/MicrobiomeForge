import pandas as pd
import pytest

from microbiomeforge import classify, parsers


def test_pick_read_len():
    avail = [50, 75, 100, 150, 200, 250, 300]
    assert classify.pick_read_len(151, avail) == 150
    assert classify.pick_read_len(149.5, avail) == 100
    assert classify.pick_read_len(30, avail) == 50
    assert classify.pick_read_len(4500, avail) == 300
    with pytest.raises(ValueError):
        classify.pick_read_len(150, [])


def test_available_read_lengths(tmp_path):
    for n in (100, 150, 75):
        (tmp_path / f"database{n}mers.kmer_distrib").write_text("")
    (tmp_path / "hash.k2d").write_text("")
    assert classify.available_read_lengths(tmp_path) == [75, 100, 150]


def test_classified_pct(tmp_path):
    rep = tmp_path / "x.k2report"
    rep.write_text("12.50\t125\t125\tU\t0\tunclassified\n87.50\t875\t0\tR\t1\troot\n")
    assert classify._classified_pct(rep) == 87.5


def test_filter_low_abundance_drops_and_renormalizes():
    m = pd.DataFrame({"taxon": ["A", "B", "noise"],
                      "s1": [0.6, 0.39995, 0.00005],
                      "s2": [0.5, 0.49995, 0.00005]})
    out = parsers.filter_low_abundance(m, 1e-4)
    assert list(out["taxon"]) == ["A", "B"]
    assert out[["s1", "s2"]].sum().round(9).tolist() == [1.0, 1.0]


def test_filter_keeps_taxon_abundant_in_one_sample():
    m = pd.DataFrame({"taxon": ["A", "rare_here"], "s1": [1.0, 0.0], "s2": [0.9, 0.1]})
    out = parsers.filter_low_abundance(m, 1e-4)
    assert "rare_here" in set(out["taxon"])


def test_filter_off():
    m = pd.DataFrame({"taxon": ["A"], "s1": [1.0]})
    assert parsers.filter_low_abundance(m, 0).equals(m)


def test_single_sample_run(tmp_path):
    """Tek örnek: PCoA atlanır, rapor yine üretilir."""
    from microbiomeforge.orchestrate import run_analysis
    from microbiomeforge.pipelines import Sample
    prof = tmp_path / "s1.bracken"
    prof.write_text("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\t"
                    "new_est_reads\tfraction_total_reads\n"
                    "A\t1\tS\t60\t0\t60\t0.6\nB\t2\tS\t40\t0\t40\t0.4\n")
    fq = tmp_path / "s1.fastq"
    fq.write_text("".join(f"@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in range(5)))
    res = run_analysis([Sample("s1", "gaita", fastq_1=str(fq), platform="illumina",
                               category="short")],
                       {"s1": str(prof)}, tmp_path / "out", tmp_path / "work")
    assert res["design"] == "single" and res["n_taxa"] == 2


def test_host_removed_from_profile(tmp_path):
    from microbiomeforge.orchestrate import run_analysis
    from microbiomeforge.pipelines import Sample
    prof = tmp_path / "s1.bracken"
    prof.write_text("name\tfraction_total_reads\nA\t0.5\nB\t0.4\nHomo sapiens\t0.1\n")
    fq = tmp_path / "s1.fastq"
    fq.write_text("".join(f"@r{i}\nACGTACGTAC\n+\nIIIIIIIIII\n" for i in range(5)))
    res = run_analysis([Sample("s1", "gaita", fastq_1=str(fq), platform="illumina",
                               category="short")],
                       {"s1": str(prof)}, tmp_path / "out", tmp_path / "work")
    assert res["host_fraction"] == {"s1": 0.1}
    prof_out = pd.read_csv(tmp_path / "out" / "02_taxonomic_profile.tsv", sep="\t")
    assert "Homo sapiens" not in set(prof_out["taxon"])
    assert abs(prof_out["s1"].sum() - 1.0) < 1e-9


def test_filter_per_group_keeps_small_groups():
    m = pd.DataFrame({"taxon": ["B1", "B2", "V1", "V2"],
                      "s": [0.7, 0.29994, 0.00004, 0.00002]})
    groups = {"B1": "Bakteri", "B2": "Bakteri", "V1": "Virüs", "V2": "Virüs"}
    assert len(parsers.filter_low_abundance(m, 1e-4)) == 2
    out = parsers.filter_low_abundance(m, 1e-4, groups=groups)
    assert set(out["taxon"]) == {"B1", "B2", "V1", "V2"}
