"""Taksonları biyolojik gruplara ayırma: Bakteri / Arke / Virüs / Mantar / Diğer ökaryot / Konak.

Kraken2 raporu (k2report) bir ağaçtır: her satırın adı girintiyle (2 boşluk = 1 düzey)
atasının altına yerleşir. Ata zinciri izlenerek her taxid bir gruba atanır ve grup
başına okuma sayısı (clade reads) toplanır. Bracken tür tablosundaki taxonomy_id
bu eşlemeyle gruplanır.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# Görüntüleme sırası
CATEGORIES = ["Bakteri", "Arke", "Virüs", "Mantar", "Diğer ökaryot", "Konak (insan)"]
MICROBIAL = ["Bakteri", "Arke", "Virüs", "Mantar", "Diğer ökaryot"]

# Öncelik: daha özgül ata önce kontrol edilir (insan ⊂ ökaryot, mantar ⊂ ökaryot).
_ANCHORS = [
    (9606, "Konak (insan)"),
    (4751, "Mantar"),
    (10239, "Virüs"),
    (2157, "Arke"),
    (2, "Bakteri"),
    (2759, "Diğer ökaryot"),
]


@dataclass
class DomainSummary:
    sample: str
    total_reads: int
    unclassified: int
    clade_reads: dict[str, int] = field(default_factory=dict)   # grup → okuma
    taxid_category: dict[int, str] = field(default_factory=dict)

    def percent_row(self) -> dict:
        row = {"sample": self.sample}
        t = self.total_reads or 1
        for c in CATEGORIES:
            row[c] = round(100.0 * self.clade_reads.get(c, 0) / t, 3)
        row["Sınıflanamayan"] = round(100.0 * self.unclassified / t, 3)
        known = sum(self.clade_reads.values()) + self.unclassified
        row["Diğer/kök düzeyi"] = round(max(0.0, 100.0 * (self.total_reads - known) / t), 3)
        return row


def parse_k2report(path: str | Path, sample: str) -> DomainSummary:
    anchors = dict(_ANCHORS)
    stack: list[int] = []          # derinliğe göre ata taxid'leri
    taxid_cat: dict[int, str] = {}
    clade: dict[str, int] = {}
    unclassified = 0
    root_reads = 0
    for line in Path(path).read_text().splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        reads, rank, taxid, name = int(parts[1]), parts[3].strip(), int(parts[4]), parts[5]
        if rank == "U":
            unclassified = reads
            continue
        depth = (len(name) - len(name.lstrip(" "))) // 2
        del stack[depth:]
        stack.append(taxid)
        if taxid == 1:
            root_reads = reads
        cat = None
        for anc in reversed(stack):  # en yakın ata önce → en özgül grup
            if anc in anchors:
                cat = anchors[anc]
                break
        if cat is not None:
            taxid_cat[taxid] = cat
        if taxid in anchors:
            clade[anchors[taxid]] = reads
    # Ökaryot clade'i mantar + insanı içerir → "Diğer ökaryot"tan düş.
    if "Diğer ökaryot" in clade:
        clade["Diğer ökaryot"] = max(
            0, clade["Diğer ökaryot"] - clade.get("Mantar", 0) - clade.get("Konak (insan)", 0))
    return DomainSummary(sample, root_reads + unclassified, unclassified, clade, taxid_cat)


def taxon_categories(bracken_path: str | Path, summary: DomainSummary) -> dict[str, str]:
    """Bracken tablosundaki her tür adı → grup (taxonomy_id üzerinden)."""
    df = pd.read_csv(bracken_path, sep="\t")
    if "taxonomy_id" not in df.columns or "name" not in df.columns:
        return {}
    return {str(n): summary.taxid_category.get(int(t), "Bilinmiyor")
            for n, t in zip(df["name"], df["taxonomy_id"])}


def summary_table(summaries: list[DomainSummary]) -> pd.DataFrame:
    return pd.DataFrame([s.percent_row() for s in summaries])


def within_category(matrix: pd.DataFrame, cats: dict[str, str], category: str) -> pd.DataFrame:
    """Profil matrisinden bir grubu seçip grup içinde yeniden 1'e normalize eder."""
    sub = matrix[matrix["taxon"].map(lambda t: cats.get(t) == category)].copy()
    if sub.empty:
        return sub
    num = [c for c in sub.columns if c != "taxon"]
    totals = sub[num].sum(axis=0).replace(0, 1.0)
    sub[num] = sub[num] / totals
    mean = sub[num].mean(axis=1)
    return sub.loc[mean.sort_values(ascending=False).index].reset_index(drop=True)


def db_has_fungi(db: str | Path, min_minimizers: int = 1_000_000) -> bool | None:
    """Kraken2 DB'nin inspect.txt'inden mantar (taxid 4751) kapsamını kontrol eder.
    inspect.txt yoksa None (bilinmiyor)."""
    p = Path(db) / "inspect.txt"
    if not p.exists():
        return None
    for line in p.read_text().splitlines():
        parts = line.split("\t")
        if len(parts) >= 6 and parts[4].strip() == "4751":
            return int(parts[1]) >= min_minimizers
    return False


LOW_SIGNAL_PCT = 0.01   # toplam okumanın %0,01'inden az → güvenilir değil


def category_note(category: str, sub: pd.DataFrame, fungi_in_db: bool | None,
                  group_pct: float | None = None) -> str:
    n = len(sub)
    if category in ("Mantar", "Diğer ökaryot") and fungi_in_db is False:
        return ("Kullanılan veritabanı mantar/protozoa genomlarını içermiyor; bu grup bu "
                "analizle değerlendirilemez (saptanmaması yokluk anlamına gelmez).")
    if n == 0:
        return "Bu grupta filtre eşiğini geçen tür saptanmadı."
    note = f"Bu grupta {n} tür saptandı; oranlar grup içi bağıl bolluktur (grup toplamı = %100)."
    if group_pct is not None and group_pct < LOW_SIGNAL_PCT:
        note = (f"**Düşük sinyal:** bu grup toplam okumanın yalnız %{group_pct:.4f}'i. Okumalar çok "
                f"sayıda ilgisiz türe birkaçar düzine dağılmış; bu desen büyük olasılıkla "
                f"yanlış-pozitif/gürültüdür ve grubun anlamlı varlığına kanıt sayılmamalıdır. " + note)
    return note
