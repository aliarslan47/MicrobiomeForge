import pandas as pd

from microbiomeforge import domains

K2 = """\
 10.00\t100\t100\tU\t0\tunclassified
 90.00\t900\t5\tR\t1\troot
 80.00\t800\t0\tR1\t131567\t  cellular organisms
 70.00\t700\t0\tD\t2\t    Bacteria
 70.00\t700\t700\tS\t111\t      Bacterium alpha
  6.00\t60\t0\tD\t2157\t    Archaea
  6.00\t60\t60\tS\t222\t      Archaeon beta
  4.00\t40\t0\tD\t2759\t    Eukaryota
  3.00\t30\t0\tK\t4751\t      Fungi
  3.00\t30\t30\tS\t333\t        Fungus gamma
  1.00\t10\t10\tS\t9606\t      Homo sapiens
  9.50\t95\t0\tR1\t10239\t  Viruses
  9.50\t95\t95\tS\t444\t    Phage delta
"""


def test_parse_k2report_groups(tmp_path):
    p = tmp_path / "s.k2report"
    p.write_text(K2)
    s = domains.parse_k2report(p, "s")
    assert s.total_reads == 1000 and s.unclassified == 100
    assert s.clade_reads["Bakteri"] == 700 and s.clade_reads["Virüs"] == 95
    assert s.clade_reads["Mantar"] == 30 and s.clade_reads["Konak (insan)"] == 10
    assert s.clade_reads["Diğer ökaryot"] == 0
    assert s.taxid_category[111] == "Bakteri" and s.taxid_category[222] == "Arke"
    assert s.taxid_category[333] == "Mantar" and s.taxid_category[9606] == "Konak (insan)"
    assert s.taxid_category[444] == "Virüs"
    row = s.percent_row()
    assert row["Bakteri"] == 70.0 and row["Sınıflanamayan"] == 10.0


def test_within_category_renormalizes():
    m = pd.DataFrame({"taxon": ["a", "b", "v"], "s": [0.6, 0.2, 0.2]})
    cats = {"a": "Bakteri", "b": "Bakteri", "v": "Virüs"}
    sub = domains.within_category(m, cats, "Bakteri")
    assert list(sub["taxon"]) == ["a", "b"]
    assert abs(sub["s"].sum() - 1.0) < 1e-9 and abs(sub["s"].iloc[0] - 0.75) < 1e-9


def test_fungi_note_when_db_lacks_fungi():
    note = domains.category_note("Mantar", pd.DataFrame(), fungi_in_db=False)
    assert "içermiyor" in note
