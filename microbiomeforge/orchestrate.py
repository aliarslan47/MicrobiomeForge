"""Uçtan uca orkestrasyon — ham okuma + nf-core çıktıları → numaralı rapor.

Bu modül tüm katmanları birleştirir:
  QC (M2) → profil ayrıştırma (M4) → istatistik (M5) → şekiller (M6) → rapor (M7)

nf-core'un fiili çalıştırılması ortama bağlıdır (Nextflow + podman). Bu orkestrasyon,
nf-core çıktıları (örnek-başına Bracken profilleri, CheckM2, AMR) hazır olduğunda onları
alıp raporu üretir; böylece rapor katmanı nf-core'dan bağımsız test edilebilir kalır.

Çıktılar (düz, numaralı, rapor son):
  01_rawdata_qc_stats.tsv, 02_taxonomic_profile.tsv, 03_diversity_stats.tsv,
  04_mag_quality.tsv, 05_amr_bgc_results.tsv, 06_figures.pdf, 07_report.pdf
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from . import figures as figmod
from . import parsers, report, stats
from .outputs import output_path, write_table
from .pipelines import Sample, plan_pipelines
from .qc_stats import compute_raw_stats, qc_delta, qc_table
from .references import ToolRegistry


def run_analysis(
    samples: list[Sample],
    profile_files: dict[str, str],       # örnek → Bracken benzeri profil dosyası
    outdir: str | Path,
    workdir: str | Path,
    checkm2_file: Optional[str] = None,
    amr_file: Optional[str] = None,
    project: str = "MicrobiomeForge",
) -> dict:
    outdir, workdir = Path(outdir), Path(workdir)
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    registry = ToolRegistry()
    groups = {s.sample: s.group for s in samples}

    # nf-core plan → hangi araçlar kullanıldı (kaynakça).
    plan_pipelines(samples, registry)

    # --- 01: Ham QC (istatistik zincirinin başı) ---
    qc_stats_list = []
    for s in samples:
        target = s.long_reads or s.fastq_1
        if target and Path(target).exists():
            registry.mark_used("seqkit")
            cat = s.category or "short"
            registry.mark_used("nanoplot" if cat == "long" else "fastp")
            qc_stats_list.append(compute_raw_stats(target, s.sample, cat, stage="raw"))
    qc_df = qc_table(qc_stats_list)
    write_table(qc_df, outdir, "raw_qc")
    delta_df = qc_delta(qc_df) if "stage" in qc_df and (qc_df["stage"] == "trimmed").any() else None

    # --- 02: Taksonomik profil ---
    per_sample = {name: parsers.parse_bracken(path, name) for name, path in profile_files.items()}
    matrix = parsers.build_abundance_matrix(per_sample)
    write_table(matrix, outdir, "taxonomic_profile")

    # --- 03: Çeşitlilik + karşılaştırma ---
    alpha = stats.alpha_diversity(matrix)
    write_table(alpha, outdir, "diversity_stats")
    at = stats.compare_alpha(alpha, groups, metric="shannon")
    D = stats.bray_curtis_matrix(matrix)
    perm = stats.permanova(D, groups)
    coords = stats.pcoa(D)
    diff = stats.differential_abundance(matrix, groups)

    # --- 04: MAG kalite ---
    mag_df = None
    if checkm2_file and Path(checkm2_file).exists():
        mag_df = parsers.parse_checkm2(checkm2_file)
        write_table(mag_df, outdir, "mag_quality")

    # --- 05: AMR ---
    amr_df = None
    if amr_file and Path(amr_file).exists():
        amr_df = parsers.parse_amr(amr_file)
        write_table(amr_df, outdir, "amr_bgc")

    # --- 06: Şekiller (numaralı) ---
    figs = figmod.Figures(workdir / "figs")
    fi = {}
    fi["raw_qc"] = figs.raw_qc(qc_df)
    fi["abundance"] = figs.abundance_bar(matrix)
    fi["alpha"] = figs.alpha_box(alpha, groups)
    fi["pcoa"] = figs.pcoa_scatter(coords, groups)
    fi["volcano"] = figs.volcano(diff)
    if mag_df is not None and not mag_df.empty:
        fi["mag"] = figs.mag_quality(mag_df)
    figs.compile_pdf(output_path(outdir, "figures"))

    # --- Anlatı metinleri (istatistikten) ---
    alpha_text = (
        f"Shannon çeşitliliği gruplar arasında {at.test} ile karşılaştırıldı "
        f"(istatistik={at.statistic:.3g}, p={at.p_value:.3g}, etki={at.effect:.3g})."
    )
    beta_text = (
        f"Bray-Curtis mesafesi üzerinde PERMANOVA gruplar arası ayrımı test etti "
        f"(pseudo-F={perm.pseudo_f:.3g}, p={perm.p_value:.3g}, {perm.permutations} permütasyon)."
    )
    n_sig = int(diff["significant"].sum()) if "significant" in diff else 0
    diff_text = (
        f"CLR-tabanlı taxon-başı test + BH-FDR ile {n_sig} takson anlamlı (q<0.05) "
        f"diferansiyel bolluk gösterdi."
    )

    # --- 07: Rapor ---
    ctx = report.ReportContext(
        project=project,
        groups=sorted({s.group for s in samples}),
        platforms=sorted({s.platform for s in samples}),
        n_samples=len(samples),
        qc_df=qc_df, top_taxa=matrix, alpha_df=alpha, diff_df=diff,
        registry=registry, manifest=figs.manifest,
        qc_delta_df=delta_df, mag_df=mag_df, amr_df=amr_df,
        alpha_test_text=alpha_text, beta_test_text=beta_text, diff_text=diff_text,
        summary=(
            f"{len(samples)} örnek ({', '.join(sorted({s.group for s in samples}))} grupları) "
            f"karşılaştırmalı olarak analiz edildi."
        ),
        fig_index=fi,
    )
    report_out = report.generate(ctx, outdir, workdir)

    return {
        "outdir": str(outdir),
        "n_samples": len(samples),
        "permanova_p": perm.p_value,
        "n_significant_taxa": n_sig,
        "report": report_out,
        "n_references": len(registry.bibliography()),
    }
