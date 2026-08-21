"""Araç referans kayıt defteri.

ÇEKIRDEK İLKE: Rapor kaynakçasındaki her atıf DOI-doğrulanmıştır (uydurma yok).
Aşağıdaki DOI'ler yayın kaynaklarından (Nature, Oxford Bioinformatics, PeerJ, PLOS,
Genome Biology, NAR GAB, bioRxiv) doğrulanarak girilmiştir; her biri
`https://doi.org/<DOI>` ile çözülebilir.

Kullanım deseni:
    reg = ToolRegistry()
    reg.mark_used("kraken2", version="2.1.3")
    reg.mark_used("sylph", version="0.6.1")
    ...
    for i, ref in enumerate(reg.bibliography(), start=1):
        print(f"[{i}] {ref.citation()}")

Yalnız `mark_used` ile işaretlenen (yani fiilen çalıştırılan) araçlar kaynakçaya girer.
Bir araç güncellenir / yeni bir araç eklenirse KATALOG'a doğrulanmış DOI'siyle eklenir.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class ToolReference:
    """Tek bir aracın kalıcı künyesi. `version` çalışma anında doldurulur."""

    key: str
    name: str
    doi: str
    year: int
    journal: str
    role: str  # pipeline içindeki görevi (ör. "taksonomik sınıflandırma")
    url: Optional[str] = None

    def doi_url(self) -> str:
        return self.url or f"https://doi.org/{self.doi}"

    def citation(self, version: Optional[str] = None) -> str:
        ver = f" (v{version})" if version else ""
        return f"{self.name}{ver}. {self.journal}, {self.year}. doi:{self.doi}"


# --------------------------------------------------------------------------- #
# KATALOG — doğrulanmış DOI'ler (web'den teyitli, 2026-08-21).
# Yeni araç / güncelleme buraya doğrulanmış DOI ile eklenir.
# --------------------------------------------------------------------------- #
CATALOG: dict[str, ToolReference] = {
    # --- Ham QC (istatistik zincirinin başı) ---
    "fastp": ToolReference(
        "fastp", "fastp: an ultra-fast all-in-one FASTQ preprocessor",
        "10.1093/bioinformatics/bty560", 2018, "Bioinformatics",
        "kısa-okuma QC/trimming",
    ),
    "seqkit": ToolReference(
        "seqkit", "SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation",
        "10.1371/journal.pone.0163962", 2016, "PLOS ONE",
        "okuma istatistiği (uzunluk/GC/sayı)",
    ),
    "nanoplot": ToolReference(
        "nanoplot", "NanoPack: visualizing and processing long-read sequencing data",
        "10.1093/bioinformatics/bty149", 2018, "Bioinformatics",
        "uzun-okuma QC/istatistik",
    ),
    "nanoplot2": ToolReference(
        "nanoplot2", "NanoPack2: population-scale evaluation of long-read sequencing data",
        "10.1093/bioinformatics/btad311", 2023, "Bioinformatics",
        "uzun-okuma QC/istatistik (v2)",
    ),
    # --- Taksonomik profilleme ---
    "kraken2": ToolReference(
        "kraken2", "Improved metagenomic analysis with Kraken 2",
        "10.1186/s13059-019-1891-0", 2019, "Genome Biology",
        "k-mer taksonomik sınıflandırma",
    ),
    "bracken": ToolReference(
        "bracken", "Bracken: estimating species abundance in metagenomics data",
        "10.7717/peerj-cs.104", 2017, "PeerJ Computer Science",
        "tür-düzeyi bolluk yeniden kestirimi",
    ),
    "sylph": ToolReference(
        "sylph", "Rapid species-level metagenome profiling and containment estimation with sylph",
        "10.1038/s41587-024-02412-y", 2024, "Nature Biotechnology",
        "hızlı tür-düzeyi profil (kısa+uzun)",
    ),
    "metaphlan4": ToolReference(
        "metaphlan4", "Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4",
        "10.1038/s41587-023-01688-w", 2023, "Nature Biotechnology",
        "marker-gen taksonomik profil",
    ),
    # --- Assembly + binning ---
    "metamdbg": ToolReference(
        "metamdbg", "High-quality metagenome assembly from long accurate reads with metaMDBG",
        "10.1038/s41587-023-01983-6", 2024, "Nature Biotechnology",
        "uzun-okuma (HiFi) metagenom assembly",
    ),
    "semibin2": ToolReference(
        "semibin2", "SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing",
        "10.1093/bioinformatics/btad209", 2023, "Bioinformatics",
        "MAG binning (kısa+uzun)",
    ),
    "checkm2": ToolReference(
        "checkm2", "CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning",
        "10.1038/s41592-023-01940-w", 2023, "Nature Methods",
        "MAG kalite (tamlık/kontaminasyon)",
    ),
    # --- İstatistik ---
    "ancombc": ToolReference(
        "ancombc", "Analysis of compositions of microbiomes with bias correction (ANCOM-BC)",
        "10.1038/s41467-020-17041-7", 2020, "Nature Communications",
        "diferansiyel bolluk testi",
    ),
    # --- Pipeline motorları ---
    "nfcore_taxprofiler": ToolReference(
        "nfcore_taxprofiler", "nf-core/taxprofiler: highly parallelised and flexible pipeline for metagenomic taxonomic classification and profiling",
        "10.1101/2023.10.20.563221", 2023, "bioRxiv",
        "taksonomik profil pipeline (kısa+uzun)",
    ),
    "nfcore_mag": ToolReference(
        "nfcore_mag", "nf-core/mag: a best-practice pipeline for metagenome hybrid assembly and binning",
        "10.1093/nargab/lqac007", 2022, "NAR Genomics and Bioinformatics",
        "assembly+binning pipeline (hibrit)",
    ),
}


class ToolRegistry:
    """Çalışma boyunca fiilen kullanılan araçları versiyonlarıyla toplar."""

    def __init__(self) -> None:
        self._used: dict[str, Optional[str]] = {}

    def mark_used(self, key: str, version: Optional[str] = None) -> None:
        if key not in CATALOG:
            raise KeyError(
                f"'{key}' KATALOG'da yok. Yeni araç ancak doğrulanmış DOI ile "
                f"references.CATALOG'a eklenerek kullanılabilir (uydurma yok)."
            )
        # Aynı araç birden çok işaretlenirse ilk bilinen versiyonu koru.
        if key not in self._used or (version and not self._used[key]):
            self._used[key] = version

    def bibliography(self) -> list[ToolReference]:
        """Kullanılan araçları yıl→ad sırasıyla döndürür (kaynakça için)."""
        refs = [CATALOG[k] for k in self._used]
        return sorted(refs, key=lambda r: (r.year, r.name))

    def version_of(self, key: str) -> Optional[str]:
        return self._used.get(key)

    def formatted_bibliography(self) -> list[str]:
        """[1] ... [2] ... numaralı, versiyonlu kaynakça satırları."""
        out = []
        for i, ref in enumerate(self.bibliography(), start=1):
            out.append(f"[{i}] {ref.citation(self.version_of(ref.key))}")
        return out
