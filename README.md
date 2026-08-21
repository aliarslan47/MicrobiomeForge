# MicrobiomeForge

Ham metagenomik okumadan (FASTQ) **karşılaştırmalı, istatistiksel ve literatüre atıflı**
bir mikrobiyom raporu üreten katman.

## Felsefe

MicrobiomeForge yeni bir hizalayıcı/sınıflandırıcı **geliştirmez**. Araç kurulumu ve
entegrasyonu problemi zaten [nf-core](https://nf-co.re) pipeline'larıyla (taxprofiler,
mag, funcscan, ampliseq) çözülmüştür. MicrobiomeForge bunları **çalışma motoru** olarak
kullanır ve özgün katkısını şu katmanda sunar:

> nf-core'un ham çıktısını (profil tablosu, MAG kalite, AMR/BGC) **numaralı-isimli
> şekiller, her adımda istatistik ve DOI-doğrulanmış kaynakça** içeren yorumlanmış bir
> **rapora** çevirmek.

## Değişmez kurallar

1. **Her şekil/tablo numaralı + isimli** (Şekil 1., Tablo 1.), metinde atıflı, yayın kalitesinde.
2. **İstatistik her zaman** vardır (sayısal + uygun test + etki büyüklüğü + p/q).
3. **İstatistik zinciri ham okumadan başlar** (okuma/taban/uzunluk/Q/GC; trim öncesi→sonrası).
4. **Kullanılan her araç** versiyon + doğrulanmış DOI ile numaralı kaynakçaya girer. Uydurma yok.
5. **Çıktı dosyaları 1'den numaralı, düz yapı, yalnız sonuçlar, rapor en son.**
6. **Karşılaştırmalı tasarım** çekirdekte: öncesi/sonrası, farklı çevreler, mikrobiyom/mikrobiyota.
7. **Platform (kısa/uzun okuma) ham datadan otomatik tespit edilir**, sistem ona göre yönlenir.

## Durum

Geliştirme aşamasında. Modül planı için `docs/` ve proje hafızasına bakınız.

## Kullanım (hedeflenen)

```bash
microbiomeforge run --samplesheet config/samplesheet.csv --outdir results/
```

Çıktı (düz, numaralı, rapor son):

```
results/
├── 01_rawdata_qc_stats.tsv
├── 02_taxonomic_profile.tsv
├── 03_diversity_stats.tsv
├── 04_mag_quality.tsv
├── 05_amr_bgc_results.tsv
├── 06_figures.pdf
└── 07_report.pdf
```
