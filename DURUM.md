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

## 2026-09-21 — AİLE PIPELINE DAG TAMAM + nf-core feasibility notu
- **✅ Pipeline DAG üretildi** (family standart, PipelineForge): `PipelineForge/specs/microbiomeforge.yml` →
  `docs/pipeline_architecture.html` (12 kart/15 kenar, self-validate ✓). Omurga girdi→QC→plan; plandan 3
  nf-core motoru (taxprofiler·mag·funcscan; **funcscan←mag MAG bağı**); yorum katmanı (P02 profil, 03
  çeşitlilik+PERMANOVA, 04 MAG-kalite, 05 AMR/BGC) → 06 şekiller → 07 rapor. 4 kategori (ortak/nf-core
  motoru/yorum/çıktı), çift-dilli TR/EN. **Görsel QA** chrome-headless-shell PNG ile yapıldı (temiz). README'lere
  DAG linki eklendi. (Pages: repo Settings→Pages→main /docs açılmalı; dosya docs/'ta hazır.)
- **nf-core feasibility (probe):** Java 11 + Docker daemon OK + 582G disk var → canlı nf-core **prensipte mümkün**.
  ⚠️ Java 11 son Nextflow için eski olabilir (17+ ister); tam koşu ağır (Nextflow kur + taxprofiler pull +
  DB/veri indir). Bounded proof-of-life henüz yapılmadı — büyük compute öncesi Ali kararı bekleniyor.

## Sırada (Ali seçecek)
- **nf-core motorunun CANLI uçtan-uca kurulumu** (Nextflow + konteyner): şu ana dek parser/stat/rapor katmanı
  canlı-OLMAYAN profil verisiyle doğrulandı; motor gerçek çalıştırılmadı. Ağır/uzun iş (bkz feasibility notu).
- funcscan AMR/BGC (05) ve mag MAG-kalite (04) modüllerinin gerçek-veri doğrulaması (04/05 real_report'ta yoktu).
