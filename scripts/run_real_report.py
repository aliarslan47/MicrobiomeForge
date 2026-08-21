#!/usr/bin/env python3
"""Gerçek veri rapor sürücüsü — PRJNA827663 (gut, öncesi/sonrası, kısa+uzun).

nf-core/taxprofiler'ın GERÇEK Bracken çıktılarını (Kraken2 2.1.5 + Bracken 3.0.1,
standard-8 DB 2025-04-02) MicrobiomeForge rapor katmanına bağlar ve numaralı çıktılar
+ 07_report.pdf üretir.

Referanslar YALNIZ fiilen çalışan araçlarla işaretlenir (uydurma yok):
nf-core/taxprofiler 1.2.4, Kraken2 2.1.5, Bracken 3.0.1, fastp (Illumina QC).
"""

from pathlib import Path

from microbiomeforge import figures as figmod
from microbiomeforge import parsers, report, stats
from microbiomeforge.outputs import output_path, write_table
from microbiomeforge.pipelines import Sample, resolve_platforms
from microbiomeforge.qc_stats import compute_raw_stats, qc_table
from microbiomeforge.references import ToolRegistry

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
BRACKEN = ROOT / "results" / "bracken" / "k2std08"
OUT = ROOT / "real_report"
WORK = Path("/tmp/mf_real_report_work")

# Örnek → (grup, ham fastq, bracken profil)
SAMPLES = [
    Sample("gut_pre_illumina", "pre",
           fastq_1=str(DATA / "pre_Day1_illumina_R1.fastq.gz"),
           fastq_2=str(DATA / "pre_Day1_illumina_R2.fastq.gz"), platform="illumina"),
    Sample("gut_post_illumina", "post",
           fastq_1=str(DATA / "post_Day30_illumina_R1.fastq.gz"),
           fastq_2=str(DATA / "post_Day30_illumina_R2.fastq.gz"), platform="illumina"),
    Sample("gut_pre_nanopore", "pre",
           long_reads=str(DATA / "pre_Day1_nanopore.fastq.gz"), platform="ont"),
    Sample("gut_post_nanopore", "post",
           long_reads=str(DATA / "post_Day30_nanopore.fastq.gz"), platform="ont"),
]
PROFILES = {
    "gut_pre_illumina": BRACKEN / "gut_pre_illumina_day1_k2std08.bracken.tsv",
    "gut_post_illumina": BRACKEN / "gut_post_illumina_day30_k2std08.bracken.tsv",
    "gut_pre_nanopore": BRACKEN / "gut_pre_nanopore_day1_k2std08.bracken.tsv",
    "gut_post_nanopore": BRACKEN / "gut_post_nanopore_day30_k2std08.bracken.tsv",
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    resolve_platforms(SAMPLES, detect=False)
    groups = {s.sample: s.group for s in SAMPLES}

    # Referanslar: yalnız GERÇEKTEN çalışan araçlar + doğru sürümler.
    reg = ToolRegistry()
    reg.mark_used("nfcore_taxprofiler", "1.2.4")
    reg.mark_used("kraken2", "2.1.5")
    reg.mark_used("bracken", "3.0.1")
    reg.mark_used("fastp")  # Illumina ön-işleme (taxprofiler)

    # 01 — Ham QC (istatistik zincirinin başı)
    qc = [compute_raw_stats(s.long_reads or s.fastq_1, s.sample, s.category or "short") for s in SAMPLES]
    qc_df = qc_table(qc)
    write_table(qc_df, OUT, "raw_qc")

    # 02 — Taksonomik profil (gerçek Bracken)
    per_sample = {n: parsers.parse_bracken(p, n) for n, p in PROFILES.items()}
    matrix = parsers.build_abundance_matrix(per_sample)
    write_table(matrix, OUT, "taxonomic_profile")

    # 03 — Çeşitlilik + karşılaştırma
    alpha = stats.alpha_diversity(matrix)
    write_table(alpha, OUT, "diversity_stats")
    at = stats.compare_alpha(alpha, groups, "shannon")
    D = stats.bray_curtis_matrix(matrix)
    perm = stats.permanova(D, groups, permutations=999)
    coords = stats.pcoa(D)
    diff = stats.differential_abundance(matrix, groups)

    # 06 — Şekiller
    figs = figmod.Figures(WORK / "figs")
    fi = {"raw_qc": figs.raw_qc(qc_df), "abundance": figs.abundance_bar(matrix, top_n=12),
          "alpha": figs.alpha_box(alpha, groups), "pcoa": figs.pcoa_scatter(coords, groups),
          "volcano": figs.volcano(diff)}
    figs.compile_pdf(output_path(OUT, "figures"))

    # Anlatı (gerçek istatistikten)
    n_sig = int(diff["significant"].sum())
    alpha_text = (f"Shannon çeşitliliği zaman noktaları (öncesi/sonrası) arasında {at.test} ile "
                  f"karşılaştırıldı (istatistik={at.statistic:.3g}, p={at.p_value:.3g}, etki={at.effect:.3g}).")
    beta_text = (f"Bray-Curtis mesafesinde PERMANOVA öncesi/sonrası ayrımını test etti "
                 f"(pseudo-F={perm.pseudo_f:.3g}, p={perm.p_value:.3g}, {perm.permutations} permütasyon). "
                 f"NOT: tek birey, zaman noktası başına biyolojik tekrar yok; ayrıca platform "
                 f"(Illumina/Nanopore) baskın teknik eksendir. Sonuçlar betimseldir.")
    diff_text = (f"CLR-tabanlı taxon-başı test + BH-FDR ile {n_sig} takson q<0.05 gösterdi "
                 f"(tek-özne tasarımda güç sınırlı).")

    # 07 — Rapor
    ctx = report.ReportContext(
        project="MicrobiomeForge — gerçek veri (PRJNA827663 gut, öncesi/sonrası, kısa+uzun)",
        groups=["pre", "post"], platforms=sorted({s.platform for s in SAMPLES}), n_samples=len(SAMPLES),
        qc_df=qc_df, top_taxa=matrix, alpha_df=alpha, diff_df=diff, registry=reg, manifest=figs.manifest,
        alpha_test_text=alpha_text, beta_test_text=beta_text, diff_text=diff_text,
        summary=("Tek sağlıklı bireyin bağırsak mikrobiyomu Gün 1 (öncesi) ve Gün 30 (sonrası) "
                 "örneklerinde, aynı örnekler hem Illumina (kısa) hem Oxford Nanopore (uzun) ile "
                 "sekanslanarak Kraken2 2.1.5 + Bracken 3.0.1 (standard-8 DB, 2025-04-02) ile "
                 "profillendi ve karşılaştırıldı."),
        fig_index=fi,
    )
    out = report.generate(ctx, OUT, WORK)
    print("Rapor:", out.get("pdf") or out.get("markdown"))
    print("PERMANOVA p =", perm.p_value, "| anlamlı takson =", n_sig,
          "| kaynak =", len(reg.bibliography()))
    print("Çıktılar:", sorted(p.name for p in OUT.iterdir()))


if __name__ == "__main__":
    main()
