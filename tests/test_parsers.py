"""M4 birim testleri: nf-core çıktı ayrıştırıcıları."""

import pandas as pd
import pytest

from microbiomeforge.parsers import (
    build_abundance_matrix,
    parse_amr,
    parse_bracken,
    parse_checkm2,
)


def test_parse_bracken_normalizes(tmp_path):
    p = tmp_path / "s1.bracken"
    pd.DataFrame({
        "name": ["E.coli", "B.fragilis"],
        "fraction_total_reads": [30.0, 10.0],  # normalize edilmemiş
    }).to_csv(p, sep="\t", index=False)
    df = parse_bracken(p, sample="s1")
    assert abs(df["abundance"].sum() - 1.0) < 1e-9
    assert set(df["taxon"]) == {"E.coli", "B.fragilis"}
    assert (df["sample"] == "s1").all()


def test_build_abundance_matrix_fills_zero():
    a = pd.DataFrame({"taxon": ["x", "y"], "abundance": [0.5, 0.5], "sample": "s1"})
    b = pd.DataFrame({"taxon": ["y", "z"], "abundance": [0.7, 0.3], "sample": "s2"})
    m = build_abundance_matrix({"s1": a, "s2": b})
    m = m.set_index("taxon")
    assert m.loc["x", "s2"] == 0.0
    assert m.loc["z", "s1"] == 0.0
    assert set(m.columns) == {"s1", "s2"}
    assert set(m.index) == {"x", "y", "z"}


def test_parse_checkm2_quality_tiers(tmp_path):
    p = tmp_path / "quality_report.tsv"
    pd.DataFrame({
        "Name": ["bin.1", "bin.2", "bin.3"],
        "Completeness": [95.0, 60.0, 20.0],
        "Contamination": [2.0, 8.0, 15.0],
    }).to_csv(p, sep="\t", index=False)
    df = parse_checkm2(p)
    tiers = dict(zip(df["bin"], df["quality"]))
    assert tiers["bin.1"] == "high"
    assert tiers["bin.2"] == "medium"
    assert tiers["bin.3"] == "low"


def test_parse_amr(tmp_path):
    p = tmp_path / "amr.tsv"
    pd.DataFrame({
        "input_file_name": ["s1", "s1"],
        "gene_symbol": ["blaTEM", "tetM"],
    }).to_csv(p, sep="\t", index=False)
    df = parse_amr(p)
    assert list(df.columns) == ["sample", "gene"]
    assert set(df["gene"]) == {"blaTEM", "tetM"}


def test_unknown_format_raises(tmp_path):
    p = tmp_path / "weird.tsv"
    pd.DataFrame({"foo": [1], "bar": [2]}).to_csv(p, sep="\t", index=False)
    with pytest.raises(ValueError):
        parse_bracken(p, sample="s1")
