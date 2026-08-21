"""M6 birim testleri: numaralı/isimli şekil üretimi."""

from pathlib import Path

import pandas as pd

from microbiomeforge.figures import Figures
from microbiomeforge.stats import alpha_diversity, bray_curtis_matrix, differential_abundance, pcoa

MATRIX = pd.DataFrame({
    "taxon": ["A", "B", "C"],
    "pre1": [0.8, 0.1, 0.1], "pre2": [0.75, 0.15, 0.1],
    "post1": [0.1, 0.8, 0.1], "post2": [0.15, 0.75, 0.1],
})
GROUPS = {"pre1": "pre", "pre2": "pre", "post1": "post", "post2": "post"}
QC = pd.DataFrame({
    "sample": ["pre1", "post1"], "stage": ["raw", "raw"],
    "n_reads": [1000, 1200], "mean_quality": [35.0, 34.0],
})
MAG = pd.DataFrame({
    "bin": ["b1", "b2"], "completeness": [95.0, 60.0],
    "contamination": [2.0, 8.0], "quality": ["high", "medium"],
})


def test_numbered_sequential_manifest(tmp_path):
    figs = Figures(tmp_path / "figs")
    n1 = figs.raw_qc(QC)
    n2 = figs.abundance_bar(MATRIX)
    assert n1 == 1 and n2 == 2
    assert [e.number for e in figs.manifest.entries] == [1, 2]
    assert all(Path(e.png_path).exists() for e in figs.manifest.entries)


def test_titles_prefixed_and_referenced(tmp_path):
    figs = Figures(tmp_path / "figs")
    figs.raw_qc(QC, title="Ham kalite")
    e = figs.manifest.entries[0]
    assert e.title == "Ham kalite"
    assert figs.manifest.ref(1) == "Şekil 1"


def test_all_figure_types_and_pdf(tmp_path):
    figs = Figures(tmp_path / "figs")
    figs.raw_qc(QC)
    figs.abundance_bar(MATRIX)
    a = alpha_diversity(MATRIX)
    figs.alpha_box(a, GROUPS)
    coords = pcoa(bray_curtis_matrix(MATRIX))
    figs.pcoa_scatter(coords, GROUPS)
    figs.mag_quality(MAG)
    figs.volcano(differential_abundance(MATRIX, GROUPS))
    assert len(figs.manifest.entries) == 6
    pdf = figs.compile_pdf(tmp_path / "06_figures.pdf")
    assert pdf.exists() and pdf.stat().st_size > 0
