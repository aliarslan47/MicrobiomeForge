# MicrobiomeForge

A layer that turns raw metagenomic reads (FASTQ) into a **comparative, statistical, literature-cited** microbiome report.

[![type](https://img.shields.io/badge/type-microbiome%20report-0d6b8f)](https://github.com/aliarslan47/MicrobiomeForge)
[![engine](https://img.shields.io/badge/engine-nf--core-2f8f5b)](https://github.com/aliarslan47/MicrobiomeForge)
[![reads](https://img.shields.io/badge/reads-short%20%C2%B7%20long-c07211)](https://github.com/aliarslan47/MicrobiomeForge)

[Türkçe](README.tr.md) · **English**

## What is it?

MicrobiomeForge is the microbiome member of the Forge family. It does **not** build a new aligner/classifier — tool installation and integration are already solved by the [nf-core](https://nf-co.re) pipelines (taxprofiler, mag, funcscan, ampliseq). MicrobiomeForge uses those as its **engine** and adds its own layer on top: turning nf-core's raw output into an interpreted report.

## What it does

It converts nf-core's raw output (profile table, MAG quality, AMR/BGC) into an interpreted report with **numbered-and-named figures, per-step statistics and a DOI-verified bibliography**. Invariant rules:

- Every figure/table is numbered + named (Figure 1., Table 1.), cited in text, publication-quality.
- Statistics are always present (value + appropriate test + effect size + p/q), and the chain starts from the raw reads (reads/base/length/Q/GC; before→after trim).
- Every tool used enters a numbered bibliography with version + verified DOI. No fabrication.
- Comparative design is core: before/after, different environments, microbiome/microbiota.
- Platform (short/long reads) is auto-detected from the raw data and the system routes accordingly.

## Installation

```bash
git clone https://github.com/aliarslan47/MicrobiomeForge.git
cd MicrobiomeForge
```

> In development. Requires an nf-core execution environment (Nextflow + a container engine); see `docs/` for the module plan and setup. Install instructions are finalized as the engine wiring lands.

## Usage

```bash
# target interface (in development)
microbiomeforge run --samplesheet config/samplesheet.csv --outdir results/
```

## Modules

Output is flat, numbered from 1, results only, report last:

| # | Output | Content |
|---|---|---|
| 01 | `01_rawdata_qc_stats.tsv` | raw-read QC statistics |
| 02 | `02_taxonomic_profile.tsv` | taxonomic profile |
| 03 | `03_diversity_stats.tsv` | alpha/beta diversity + stats |
| 04 | `04_mag_quality.tsv` | MAG quality |
| 05 | `05_amr_bgc_results.tsv` | AMR / BGC results |
| 06 | `06_figures.pdf` | numbered, named figures |
| 07 | `07_report.pdf` | interpreted final report |

The full module plan lives in `docs/` and the project memory.

---

Forge family: [RNAForge](https://github.com/aliarslan47/RNAForge) (bulk RNA-seq) · [BacForge](https://github.com/aliarslan47/BacForge) (bacteria) · [VirusForge](https://github.com/aliarslan47/VirusForge) (virus/phage) · **MicrobiomeForge** (microbiome) · [Vaxforge](https://github.com/aliarslan47/Vaxforge) (reverse vaccinology) · [ImmForge](https://github.com/aliarslan47/ImmForge) (immune simulation) · [PipelineForge](https://github.com/aliarslan47/PipelineForge) (DAG generator). Licensed under [MIT](LICENSE).
