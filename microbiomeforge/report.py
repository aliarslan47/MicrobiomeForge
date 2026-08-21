"""Rapor birleştirici — Jinja2 → Markdown → PDF (pandoc/xelatex).

Rapor kuralları burada uygulanır:
  * Her tablo "Tablo N." ile numaralı+isimli, metinde atıflı.
  * Her şekil "Şekil N." (M6 manifest'inden) gömülü + atıflı.
  * İstatistik ham okumadan (bölüm 2) başlar, yukarı akar.
  * Kaynakça yalnız FİİLEN kullanılan araçları versiyon + DOI ile listeler (ToolRegistry).

Çıktı: 07_report.pdf (en son numara). Markdown ara dosyası work-dir'de tutulur.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Optional

import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

from .figures import FigureEntry, FigureManifest
from .outputs import output_path
from .references import ToolRegistry

_TEMPLATE_DIR = Path(__file__).parent / "templates"


# --------------------------------------------------------------------------- #
# Markdown yardımcıları
# --------------------------------------------------------------------------- #
def df_to_md(df: pd.DataFrame, max_rows: Optional[int] = None) -> str:
    """DataFrame → GitHub-flavored markdown tablosu (tabulate'siz)."""
    if df is None or df.empty:
        return "_(veri yok)_"
    if max_rows is not None:
        df = df.head(max_rows)
    cols = list(df.columns)
    header = "| " + " | ".join(str(c) for c in cols) + " |"
    sep = "| " + " | ".join("---" for _ in cols) + " |"
    lines = [header, sep]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(_fmt(row[c]) for c in cols) + " |")
    return "\n".join(lines)


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.4g}"
    return str(v)


class TableNumberer:
    """Ardışık 'Tablo N.' başlıkları üretir."""

    def __init__(self) -> None:
        self._n = 0

    def caption(self, title: str, df: pd.DataFrame, max_rows: Optional[int] = None) -> str:
        self._n += 1
        return f"**Tablo {self._n}.** {title}\n\n{df_to_md(df, max_rows)}"


# --------------------------------------------------------------------------- #
# Bağlam
# --------------------------------------------------------------------------- #
@dataclass
class ReportContext:
    project: str
    groups: list[str]
    platforms: list[str]
    n_samples: int
    qc_df: pd.DataFrame
    top_taxa: pd.DataFrame
    alpha_df: pd.DataFrame
    diff_df: pd.DataFrame
    registry: ToolRegistry
    manifest: FigureManifest
    qc_delta_df: Optional[pd.DataFrame] = None
    mag_df: Optional[pd.DataFrame] = None
    amr_df: Optional[pd.DataFrame] = None
    alpha_test_text: str = ""
    beta_test_text: str = ""
    diff_text: str = ""
    summary: str = ""
    methods: str = ""
    fig_index: dict[str, int] = field(default_factory=dict)  # rol → şekil no


def _fig(manifest: FigureManifest, number: Optional[int]) -> Optional[FigureEntry]:
    if number is None:
        return None
    for e in manifest.entries:
        if e.number == number:
            return e
    return None


def render_markdown(ctx: ReportContext) -> str:
    env = Environment(
        loader=FileSystemLoader(str(_TEMPLATE_DIR)),
        autoescape=select_autoescape(enabled_extensions=()),
        trim_blocks=True, lstrip_blocks=True,
    )
    tpl = env.get_template("report.md.j2")
    tn = TableNumberer()
    fi = ctx.fig_index

    data = {
        "project": ctx.project,
        "date": date.today().isoformat(),
        "n_samples": ctx.n_samples,
        "groups": ctx.groups,
        "platforms": ctx.platforms,
        "summary": ctx.summary or "Karşılaştırmalı mikrobiyom analizi sonuçları aşağıda özetlenmiştir.",
        "table_raw_qc": tn.caption("Ham okuma istatistiği (örnek başına)", ctx.qc_df),
        "table_qc_delta": tn.caption("Trim öncesi→sonrası tutulma oranı", ctx.qc_delta_df)
            if ctx.qc_delta_df is not None and not ctx.qc_delta_df.empty else "",
        "table_top_taxa": tn.caption("En bol taksonlar (bağıl bolluk)", ctx.top_taxa, max_rows=15),
        "alpha_narrative": ctx.alpha_test_text,
        "table_alpha": tn.caption("Örnek başına alfa çeşitlilik", ctx.alpha_df),
        "beta_narrative": ctx.beta_test_text,
        "diff_narrative": ctx.diff_text,
        "table_diff": tn.caption("Diferansiyel bol taksonlar (q<0.05)", ctx.diff_df[ctx.diff_df["significant"]])
            if "significant" in ctx.diff_df.columns and ctx.diff_df["significant"].any() else "",
        "table_mag": tn.caption("MAG kalite özeti", ctx.mag_df)
            if ctx.mag_df is not None and not ctx.mag_df.empty else "",
        "table_amr": tn.caption("Saptanan AMR genleri", ctx.amr_df)
            if ctx.amr_df is not None and not ctx.amr_df.empty else "",
        "methods": ctx.methods or _default_methods(ctx),
        "references": ctx.registry.formatted_bibliography(),
        "fig_raw_qc": _fig(ctx.manifest, fi.get("raw_qc")),
        "fig_abundance": _fig(ctx.manifest, fi.get("abundance")),
        "fig_alpha": _fig(ctx.manifest, fi.get("alpha")),
        "fig_pcoa": _fig(ctx.manifest, fi.get("pcoa")),
        "fig_volcano": _fig(ctx.manifest, fi.get("volcano")),
        "fig_mag": _fig(ctx.manifest, fi.get("mag")),
    }
    return tpl.render(**data)


def _default_methods(ctx: ReportContext) -> str:
    tools = ", ".join(r.name.split(":")[0].split(" ")[0] for r in ctx.registry.bibliography())
    return (
        "Ham okumalar platform (kısa/uzun) açısından otomatik sınıflandırıldı; kalite "
        "istatistikleri okuma düzeyinden hesaplandı. Taksonomik profilleme ve derleme "
        f"nf-core pipeline'larıyla yürütüldü. Kullanılan araçlar: {tools}. Alfa/beta "
        "çeşitlilik, PERMANOVA ve CLR-tabanlı diferansiyel bolluk (BH-FDR) analizleri "
        "uygulandı."
    )


def markdown_to_pdf(md_path: Path, pdf_path: Path) -> bool:
    """pandoc + xelatex ile PDF üretir. pandoc yoksa False döner (MD kalır)."""
    if shutil.which("pandoc") is None:
        return False
    cmd = [
        "pandoc", str(md_path), "-o", str(pdf_path),
        "--pdf-engine=xelatex", "-V", "geometry:margin=2cm",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, timeout=300)
        return pdf_path.exists()
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False


def generate(ctx: ReportContext, outdir: str | Path, workdir: str | Path) -> dict:
    """Raporu üretir: work-dir'de MD, çıktı dizininde 07_report.pdf (yoksa .md)."""
    outdir = Path(outdir)
    workdir = Path(workdir)
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    md = render_markdown(ctx)
    md_path = workdir / "07_report.md"
    md_path.write_text(md, encoding="utf-8")

    pdf_path = output_path(outdir, "report")
    ok = markdown_to_pdf(md_path, pdf_path)
    if not ok:
        # PDF motoru yoksa markdown'ı çıktıya kopyala (yine numaralı: 07_report.md).
        fallback = outdir / "07_report.md"
        fallback.write_text(md, encoding="utf-8")
        return {"pdf": None, "markdown": str(fallback), "md_source": str(md_path)}
    return {"pdf": str(pdf_path), "markdown": str(md_path)}
