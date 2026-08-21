"""M7 birim testleri: rapor birleştirici (Markdown + PDF)."""

import shutil

import pandas as pd
import pytest

from microbiomeforge.figures import Figures
from microbiomeforge.references import ToolRegistry
from microbiomeforge.report import ReportContext, df_to_md, generate, render_markdown
from microbiomeforge.stats import (
    alpha_diversity,
    bray_curtis_matrix,
    differential_abundance,
    pcoa,
)

MATRIX = pd.DataFrame({
    "taxon": ["A", "B", "C"],
    "pre1": [0.8, 0.1, 0.1], "pre2": [0.75, 0.15, 0.1],
    "post1": [0.1, 0.8, 0.1], "post2": [0.15, 0.75, 0.1],
})
GROUPS = {"pre1": "pre", "pre2": "pre", "post1": "post", "post2": "post"}
QC = pd.DataFrame({
    "sample": ["pre1", "post1"], "stage": ["raw", "raw"], "category": ["short", "short"],
    "n_reads": [1000, 1200], "total_bases": [150000, 180000],
    "min_length": [150, 150], "median_length": [150.0, 150.0], "mean_length": [150.0, 150.0],
    "max_length": [150, 150], "n50": [150, 150], "mean_quality": [35.0, 34.0], "gc_percent": [50.0, 51.0],
})


def _build_ctx(tmp_path):
    reg = ToolRegistry()
    reg.mark_used("kraken2", "2.1.3")
    reg.mark_used("sylph", "0.6.1")
    figs = Figures(tmp_path / "figs")
    fi = {}
    fi["raw_qc"] = figs.raw_qc(QC)
    fi["abundance"] = figs.abundance_bar(MATRIX)
    a = alpha_diversity(MATRIX)
    fi["alpha"] = figs.alpha_box(a, GROUPS)
    fi["pcoa"] = figs.pcoa_scatter(pcoa(bray_curtis_matrix(MATRIX)), GROUPS)
    diff = differential_abundance(MATRIX, GROUPS)
    fi["volcano"] = figs.volcano(diff)
    return ReportContext(
        project="TestRun", groups=["pre", "post"], platforms=["illumina"], n_samples=4,
        qc_df=QC, top_taxa=MATRIX, alpha_df=a, diff_df=diff,
        registry=reg, manifest=figs.manifest,
        alpha_test_text="Gruplar arası Shannon farkı test edildi.",
        beta_test_text="PERMANOVA ile gruplar ayrıştı.",
        diff_text="CLR-tabanlı diferansiyel bolluk uygulandı.",
        fig_index=fi,
    )


def test_df_to_md_basic():
    md = df_to_md(pd.DataFrame({"a": [1], "b": [2.5]}))
    assert "| a | b |" in md
    assert "| --- | --- |" in md


def test_render_markdown_has_numbered_tables_and_refs(tmp_path):
    ctx = _build_ctx(tmp_path)
    md = render_markdown(ctx)
    assert "**Tablo 1.**" in md
    assert "**Tablo 2.**" in md          # birden çok numaralı tablo
    assert "Şekil 1." in md              # gömülü şekil atfı
    assert "## 8. Kaynakça" in md
    assert "doi:10.1038/s41587-024-02412-y" in md  # sylph DOI kaynakçada
    assert "## 2. Ham veri kalitesi" in md         # istatistik zinciri ham okumadan


def test_references_only_used_tools(tmp_path):
    ctx = _build_ctx(tmp_path)
    md = render_markdown(ctx)
    assert "Kraken 2" in md
    assert "MetaPhlAn" not in md  # kullanılmayan araç kaynakçada yok


@pytest.mark.skipif(shutil.which("pandoc") is None, reason="pandoc yok")
def test_generate_pdf(tmp_path):
    ctx = _build_ctx(tmp_path)
    out = generate(ctx, outdir=tmp_path / "results", workdir=tmp_path / "work")
    # xelatex mevcutsa PDF, değilse markdown fallback — ikisi de 07_ ile başlar.
    produced = out["pdf"] or out["markdown"]
    assert produced is not None
    assert "07_report" in produced
