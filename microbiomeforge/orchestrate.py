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

import sys
from pathlib import Path
from typing import Optional

import pandas as pd

from . import figures as figmod
from . import domains, parsers, report, stats
from .outputs import output_path, write_table
from .pipelines import DesignKind, Sample, plan_pipelines, resolve_design
from .qc_stats import compute_raw_stats, qc_delta, qc_table
from .references import ToolRegistry


def _find_k2report(profile: Path, sample: str) -> Optional[Path]:
    """Bracken profilinin yanındaki Kraken2 raporunu bulur (yerel motor veya nf-core)."""
    for ext in (".k2report", ".kraken2.report.txt", ".kreport", ".report"):
        p = profile.parent / f"{sample}{ext}"
        if p.exists():
            return p
    return None


def _db_text(kraken_db, classify_rows, host_fraction, min_abundance, n_raw, n_kept,
             min_group_abundance=0.005, low_groups=()) -> str:
    parts = []
    if kraken_db:
        parts.append(f"Kraken2 veritabanı: `{Path(kraken_db).name}`.")
    for r in classify_rows:
        parts.append(f"{r['sample']}: okumaların %{r['classified_pct']} kadarı sınıflandırıldı.")
    if host_fraction and any(v > 0 for v in host_fraction.values()):
        hf = ", ".join(f"{k} %{v * 100:.3f}" for k, v in host_fraction.items())
        parts.append(f"Konak (insan) okumaları mikrobiyom profilinden çıkarıldı (Bracken payı: {hf}).")
    if min_abundance > 0:
        parts.append(f"Tüm profilde %{min_abundance * 100:g} veya kendi grubu içinde "
                     f"%{min_group_abundance * 100:g} bağıl bolluğa hiçbir örnekte ulaşmayan türler "
                     f"gürültü olarak elendi ({n_raw} → {n_kept} tür).")
    if low_groups:
        parts.append(f"Toplam okumanın %{domains.LOW_SIGNAL_PCT:g}'inden azını oluşturan gruplar "
                     f"({', '.join(low_groups)}) düşük sinyal olarak işaretlendi ve çeşitlilik "
                     f"hesaplarına katılmadı.")
    return " ".join(parts)


def run_analysis(
    samples: list[Sample],
    profile_files: dict[str, str],       # örnek → Bracken benzeri profil dosyası
    outdir: str | Path,
    workdir: str | Path,
    checkm2_file: Optional[str] = None,
    amr_file: Optional[str] = None,
    project: str = "MicrobiomeForge",
    design: Optional[DesignKind] = None,  # None=metadata'dan çıkar; "single"|"comparative"=açık bildirim
    kraken_db: Optional[str] = None,      # verilirse yerel motor: FASTQ → Kraken2 → Bracken
    threads: int = 8,
    min_abundance: float = 0.0,           # bu bağıl bolluğun altındaki taksonlar elenir (ör. 1e-4)
    bracken_threshold: int = 10,          # Bracken -t: taksonun en az okuma sayısı
    memory_mapping: bool = False,
    host_taxa: tuple[str, ...] = ("Homo sapiens",),
    min_group_abundance: float = 0.005,   # grup içi eşik (küçük gruplar için)  # mikrobiyom profilinden çıkarılır, ayrıca raporlanır
) -> dict:
    outdir, workdir = Path(outdir), Path(workdir)
    outdir.mkdir(parents=True, exist_ok=True)
    workdir.mkdir(parents=True, exist_ok=True)

    registry = ToolRegistry()
    groups = {s.sample: s.group for s in samples}

    # Tasarım çözümü: metadata (group sütunu) ya da açık bildirim (design=...).
    dspec = resolve_design(samples, declared=design)
    # Kararı gürültülü şekilde bildir (sessiz varsayım yok).
    print(f"[MicrobiomeForge] Tasarım: {dspec.kind} — {dspec.note}", file=sys.stderr)

    # Kaynakçaya YALNIZ fiilen çalışan araçlar girer.
    # Yerel motorda nf-core planı işlenmez; nf-core modunda plan (mevcut davranış).
    if kraken_db is None:
        plan_pipelines(samples, registry)

    # --- 01: Ham QC (istatistik zincirinin başı) ---
    qc_stats_list = []
    for s in samples:
        target = s.long_reads or s.fastq_1
        if target and Path(target).exists():
            cat = s.category or "short"
            # QC istatistikleri Python içinde hesaplanır; seqkit/fastp/NanoPlot çalışmaz.
            # nf-core modunda bu araçlar plan üzerinden zaten kaynakçaya girer.
            qc_stats_list.append(compute_raw_stats(target, s.sample, cat, stage="raw"))
    qc_df = qc_table(qc_stats_list)
    write_table(qc_df, outdir, "raw_qc")
    delta_df = qc_delta(qc_df) if "stage" in qc_df and (qc_df["stage"] == "trimmed").any() else None

    # --- 02a: Yerel sınıflandırma (FASTQ → Kraken2 → Bracken) ---
    classify_rows = []
    if kraken_db is not None:
        from . import classify
        classify.check_tools()
        db = Path(kraken_db)
        # Bracken dağılımı için ortanca uzunluk (kırpılmış birkaç kısa okuma ortalamayı düşürür).
        mean_len = {st.sample: st.median_length for st in qc_stats_list}
        profile_files = dict(profile_files or {})
        for s in samples:
            if s.sample in profile_files or not s.has_reads():
                continue
            r = classify.classify_sample(
                s, db, workdir / "kraken2", mean_len.get(s.sample, 150.0), threads=threads,
                bracken_threshold=bracken_threshold, memory_mapping=memory_mapping,
            )
            profile_files[s.sample] = str(r.bracken)
            classify_rows.append({"sample": r.sample, "classified_pct": r.classified_pct,
                                  "bracken_read_len": r.read_len})
        registry.mark_used("kraken2", classify.tool_version("kraken2"))
        registry.mark_used("bracken", classify.tool_version("bracken"))
        for row in classify_rows:
            print(f"[MicrobiomeForge] {row['sample']}: sınıflanan okuma %{row['classified_pct']}",
                  file=sys.stderr)

    # --- 02: Taksonomik profil ---
    per_sample = {name: parsers.parse_bracken(path, name) for name, path in profile_files.items()}
    matrix = parsers.build_abundance_matrix(per_sample)

    # Biyolojik gruplar (Bakteri/Arke/Virüs/Mantar/...) — Kraken2 raporu varsa.
    domain_summaries, taxon_cat = [], {}
    for name, path in profile_files.items():
        rep_path = _find_k2report(Path(path), name)
        if rep_path is None:
            continue
        ds = domains.parse_k2report(rep_path, name)
        domain_summaries.append(ds)
        for t, c in domains.taxon_categories(path, ds).items():
            taxon_cat.setdefault(t, c)
    domain_df = domains.summary_table(domain_summaries) if domain_summaries else None
    if domain_df is not None:
        write_table(domain_df, outdir, "taxonomic_groups")

    host_fraction = {}
    if host_taxa and not matrix.empty:
        is_host = matrix["taxon"].isin(host_taxa) | (
            matrix["taxon"].map(taxon_cat.get) == "Konak (insan)")
        num_cols = [c for c in matrix.columns if c != "taxon"]
        host_fraction = {c: round(float(matrix.loc[is_host, c].sum()), 6) for c in num_cols}
        if is_host.any():
            matrix = matrix[~is_host].reset_index(drop=True)
            totals = matrix[num_cols].sum(axis=0).replace(0, 1.0)
            matrix[num_cols] = matrix[num_cols] / totals
            print(f"[MicrobiomeForge] Konak ({', '.join(host_taxa)}) profilden çıkarıldı; "
                  f"Bracken payı: {host_fraction}", file=sys.stderr)
    n_taxa_raw = len(matrix)
    if min_abundance > 0:
        matrix = parsers.filter_low_abundance(matrix, min_abundance, groups=taxon_cat or None,
                                              min_group_abundance=min_group_abundance)
        print(f"[MicrobiomeForge] Düşük bolluk filtresi (<{min_abundance:g}): "
              f"{n_taxa_raw} → {len(matrix)} takson", file=sys.stderr)
    if not matrix.empty:  # en bol takson üstte
        _num = matrix.select_dtypes("number")
        matrix = matrix.loc[_num.mean(axis=1).sort_values(ascending=False).index].reset_index(drop=True)
    if taxon_cat:
        _out = matrix.copy()
        _out.insert(1, "group", _out["taxon"].map(lambda t: taxon_cat.get(t, "Bilinmiyor")))
        write_table(_out, outdir, "taxonomic_profile")
    else:
        write_table(matrix, outdir, "taxonomic_profile")

    # Düşük sinyalli gruplar (toplam okumanın %0,01'inden azı) ayrı bölümde uyarıyla
    # gösterilir ama çeşitlilik/istatistiğe katılmaz (gürültü zenginliği şişirmesin).
    low_groups = []
    if domain_df is not None:
        low_groups = [c for c in domains.MICROBIAL
                      if c in domain_df and 0 < domain_df[c].mean() < domains.LOW_SIGNAL_PCT]
    full_matrix = matrix
    if low_groups and taxon_cat:
        keep = ~matrix["taxon"].map(lambda t: taxon_cat.get(t) in low_groups)
        matrix = matrix[keep].reset_index(drop=True)
        _nc = [c for c in matrix.columns if c != "taxon"]
        matrix[_nc] = matrix[_nc] / matrix[_nc].sum(axis=0).replace(0, 1.0)
        print(f"[MicrobiomeForge] Düşük sinyalli gruplar çeşitliliğe katılmadı: {low_groups}",
              file=sys.stderr)

    # --- 03: Çeşitlilik (+ tasarım karşılaştırmalıysa grup testleri) ---
    alpha = stats.alpha_diversity(matrix)
    write_table(alpha, outdir, "diversity_stats")
    n_profiled = len([c for c in matrix.columns if c != "taxon"])
    D = stats.bray_curtis_matrix(matrix)
    # Beta çeşitlilik ordinasyonu en az 3 örnekle anlamlı; 1-2 örnekte atlanır.
    coords = stats.pcoa(D) if n_profiled >= 3 else None
    if dspec.kind == "comparative":
        at = stats.compare_alpha(alpha, groups, metric="shannon")
        perm = stats.permanova(D, groups)
        diff = stats.differential_abundance(matrix, groups)
    else:
        # Tekli/betimsel: gruplar arası karşılaştırma yapılmaz.
        at = perm = None
        diff = pd.DataFrame()

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
    category_sections = []
    if taxon_cat:
        fungi_in_db = domains.db_has_fungi(kraken_db) if kraken_db else None
        for cat in domains.MICROBIAL:
            sub = domains.within_category(full_matrix, taxon_cat, cat)
            fig_no = None
            if not sub.empty:
                fig_no = figs.abundance_bar(sub, title=f"{cat}: en bol türler (grup içi bağıl bolluk)")
            category_sections.append({
                "category": cat, "df": sub, "fig": fig_no,
                "note": domains.category_note(
                    cat, sub, fungi_in_db,
                    float(domain_df[cat].mean()) if domain_df is not None and cat in domain_df else None),
            })
    else:
        fi["abundance"] = figs.abundance_bar(matrix)
    fi["alpha"] = figs.alpha_box(alpha, groups)
    if coords is not None:
        fi["pcoa"] = figs.pcoa_scatter(coords, groups)
    if dspec.kind == "comparative" and not diff.empty:
        fi["volcano"] = figs.volcano(diff)
    if mag_df is not None and not mag_df.empty:
        fi["mag"] = figs.mag_quality(mag_df)
    figs.compile_pdf(output_path(outdir, "figures"))

    # --- Anlatı metinleri (istatistikten) ---
    if dspec.kind == "comparative":
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
    else:
        n_sig = 0
        alpha_text = ("Tekli/betimsel tasarım: alfa çeşitliliği örnekler için betimsel "
                      "olarak raporlandı; gruplar arası karşılaştırma yapılmadı.")
        beta_text = ("Betimsel beta-çeşitlilik (Bray-Curtis + PCoA) sunuldu; "
                     "PERMANOVA gibi grup-ayrım testleri uygulanmadı." if coords is not None else
                     f"Beta-çeşitlilik ordinasyonu {n_profiled} örnekle anlamlı olmadığından "
                     "uygulanmadı (en az 3 örnek gerekir).")
        diff_text = "Diferansiyel bolluk testi tekli tasarımda uygulanmadı."

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
            f"{len(samples)} örnek ({', '.join(dspec.groups)} grubu/grupları) "
            f"{'karşılaştırmalı' if dspec.kind == 'comparative' else 'tekli/betimsel'} "
            f"olarak analiz edildi. {dspec.note}"
        ),
        fig_index=fi,
        domain_df=domain_df,
        category_sections=category_sections,
        db_text=_db_text(kraken_db, classify_rows, host_fraction, min_abundance, n_taxa_raw,
                         len(full_matrix), min_group_abundance, low_groups),
    )
    report_out = report.generate(ctx, outdir, workdir)

    return {
        "outdir": str(outdir),
        "n_samples": len(samples),
        "design": dspec.kind,
        "design_note": dspec.note,
        "permanova_p": perm.p_value if perm is not None else None,
        "n_significant_taxa": n_sig,
        "report": report_out,
        "n_references": len(registry.bibliography()),
        "n_taxa_raw": n_taxa_raw,
        "n_taxa": len(full_matrix),
        "classification": classify_rows,
        "host_fraction": host_fraction,
    }
