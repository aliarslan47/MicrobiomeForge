# MicrobiomeForge

Ham metagenomik okumadan (FASTQ) **karşılaştırmalı, istatistiksel ve literatüre atıflı** bir mikrobiyom raporu üreten katman.

[![type](https://img.shields.io/badge/type-microbiome%20report-0d6b8f)](https://github.com/aliarslan47/MicrobiomeForge)
[![engine](https://img.shields.io/badge/engine-nf--core-2f8f5b)](https://github.com/aliarslan47/MicrobiomeForge)
[![reads](https://img.shields.io/badge/reads-short%20%C2%B7%20long-c07211)](https://github.com/aliarslan47/MicrobiomeForge)

**Türkçe** · [English](README.md)

## Nedir?

MicrobiomeForge, Forge ailesinin mikrobiyom üyesidir. Yeni bir hizalayıcı/sınıflandırıcı **geliştirmez** — araç kurulumu ve entegrasyonu problemi zaten [nf-core](https://nf-co.re) pipeline'larıyla (taxprofiler, mag, funcscan, ampliseq) çözülmüştür. MicrobiomeForge bunları **çalışma motoru** olarak kullanır ve üstüne kendi katmanını ekler: nf-core'un ham çıktısını yorumlanmış bir rapora çevirmek.

## Ne yapar?

nf-core'un ham çıktısını (profil tablosu, MAG kalite, AMR/BGC) **numaralı-isimli şekiller, her adımda istatistik ve DOI-doğrulanmış kaynakça** içeren yorumlanmış bir rapora çevirir. Değişmez kurallar:

- Her şekil/tablo numaralı + isimli (Şekil 1., Tablo 1.), metinde atıflı, yayın kalitesinde.
- İstatistik her zaman vardır (değer + uygun test + etki büyüklüğü + p/q) ve zincir ham okumadan başlar (okuma/taban/uzunluk/Q/GC; trim öncesi→sonrası).
- Kullanılan her araç versiyon + doğrulanmış DOI ile numaralı kaynakçaya girer. Uydurma yok.
- Karşılaştırmalı tasarım çekirdektedir: öncesi/sonrası, farklı çevreler, mikrobiyom/mikrobiyota.
- Platform (kısa/uzun okuma) ham datadan otomatik tespit edilir ve sistem ona göre yönlenir.

## Kurulum

```bash
git clone https://github.com/aliarslan47/MicrobiomeForge.git
cd MicrobiomeForge
```

> Geliştirme aşamasında. Bir nf-core çalışma ortamı gerektirir (Nextflow + bir konteyner motoru); modül planı ve kurulum için `docs/`'a bakın. Kurulum talimatları motor bağlanması ilerledikçe kesinleşir.

## Kullanım

```bash
# hedeflenen arayüz (geliştirme aşamasında)
microbiomeforge run --samplesheet config/samplesheet.csv --outdir results/
```

## Modüller

Çıktı düz, 1'den numaralı, yalnız sonuçlar, rapor en son:

| # | Çıktı | İçerik |
|---|---|---|
| 01 | `01_rawdata_qc_stats.tsv` | ham-okuma QC istatistikleri |
| 02 | `02_taxonomic_profile.tsv` | taksonomik profil |
| 03 | `03_diversity_stats.tsv` | alfa/beta çeşitlilik + istatistik |
| 04 | `04_mag_quality.tsv` | MAG kalitesi |
| 05 | `05_amr_bgc_results.tsv` | AMR / BGC sonuçları |
| 06 | `06_figures.pdf` | numaralı, isimli şekiller |
| 07 | `07_report.pdf` | yorumlanmış nihai rapor |

Tam modül planı `docs/`'ta ve proje hafızasındadır.

---

Forge ailesi: [RNAForge](https://github.com/aliarslan47/RNAForge) (bulk RNA-seq) · [BacForge](https://github.com/aliarslan47/BacForge) (bakteri) · [VirusForge](https://github.com/aliarslan47/VirusForge) (virüs/faj) · **MicrobiomeForge** (mikrobiyom) · [Vaxforge](https://github.com/aliarslan47/Vaxforge) (ters aşılama) · [ImmForge](https://github.com/aliarslan47/ImmForge) (bağışıklık simülasyonu) · [PipelineForge](https://github.com/aliarslan47/PipelineForge) (DAG üreticisi). [MIT](LICENSE) lisansı altında.
