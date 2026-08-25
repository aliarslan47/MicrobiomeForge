# DURUM — MicrobiomeForge

> Bu dosya "nerede kaldık" anlık görüntüsüdür. Tüm modül/karar detayı Claude belleğindedir
> (`microbiomeforge-project` memory). Claude bunu anlamlı her durakta ve `/clear` öncesi günceller.

**Konum:** `/home/ali/MicrobiomeForge/`
**GitHub:** `github.com/aliarslan47/MicrobiomeForge` (Forge ailesi mikrobiyom üyesi)
**Son güncelleme:** 2026-08-25 (tasarım çözümü eklendi)

## Nedir
Kendi hizalayıcı/sınıflandırıcı yazmaz. **nf-core motoru** (taxprofiler, mag, funcscan,
ampliseq) üstüne "yorumlanmış rapor" katmanı ekler: numaralı+isimli şekil, adım-adım
istatistik (değer+test+etki+p/q), DOI-doğrulanmış kaynakça, karşılaştırmalı tasarım,
platform (kısa/uzun okuma) otomatik tespiti.

## Çıktı (düz, 1'den numaralı, rapor en sonda)
`01_rawdata_qc_stats.tsv` · `02_taxonomic_profile.tsv` · `03_diversity_stats.tsv` ·
`04_mag_quality.tsv` · `05_amr_bgc_results.tsv` · `06_figures.pdf` · `07_report.pdf`

## Şu an nerede kaldık
- **M0–M8 TAMAM** (iskele+ref kayıt · platform tespit · QC-stat · samplesheet+nf-core planlama ·
  çıktı ayrıştırıcı · istatistik motoru [çeşitlilik+PERMANOVA+diferansiyel bolluk] · yayın-kalite
  şekil · rapor birleştirici Jinja2→MD→PDF/xelatex · CLI+uçtan-uca orkestrasyon+entegrasyon testi).
- **43/43 test geçiyor** (`python -m pytest -q`). Yalnız matplotlib deprecation uyarısı (hata değil).
- **GERÇEK KOŞU doğrulandı:** PRJNA827663 gut öncesi/sonrası (kısa+uzun), tam rapor →
  `real_report/` (01_qc, 02_taxonomic_profile, 03_diversity_stats + 06_figures.pdf + 07_report.pdf).
- README Forge standardına göre çift dilli (EN + `README.tr.md`).
- **TASARIM ÇÖZÜMÜ (2026-08-25):** tekli/betimsel vs karşılaştırmalı artık ya `group`
  metadata'sından çıkarılıyor ya da `--design single|comparative` ile açıkça bildiriliyor
  (`pipelines.resolve_design` → `DesignSpec`). Tutarsızlık gürültülü hata (comparative+tek grup
  → ValueError/rc!=0); tekli tasarımda grup testleri (compare_alpha/PERMANOVA/diferansiyel bolluk)
  atlanıp rapora+loga not düşülüyor. Platform tespiti (`detect.py`) veriden; tasarım ise
  metadata/bildirim — çünkü tasarım verinin fiziksel özelliği değil. **50/50 test geçiyor.**

## Sırada (Ali seçecek)
- nf-core motor bağlantısının (Nextflow + konteyner) canlı uçtan-uca kurulumu / `docs/` modül planı
- funcscan AMR/BGC (05) ve mag MAG-kalite (04) modüllerinin gerçek-veri doğrulaması
- Aile Pipeline DAG şeması (artifact) + PipelineForge spec
