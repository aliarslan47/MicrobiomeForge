"""MicrobiomeForge komut satırı arayüzü.

    microbiomeforge run \
        --samplesheet config/samplesheet.csv \
        --profiles-dir nf_out/profiles \
        --outdir results/ [--checkm2 quality_report.tsv] [--amr amr.tsv]

Yerel motor (nf-core/Docker gerekmez) — ham FASTQ'dan başlar:

    microbiomeforge run --samplesheet s.csv --kraken-db ~/databases/kraken2 \
        --outdir results/ [--threads 12] [--min-abundance 0.0001]

`--profiles-dir` verilirse her örnek için `<sample>.bracken` (veya .tsv) profil dosyası
okunur (ör. nf-core/taxprofiler çıktısı). `--kraken-db` verilirse profili olmayan
örnekler Kraken2 → Bracken ile yerelde sınıflandırılır. Platform (kısa/uzun) ham
okumadan otomatik tespit edilir (örnek sayfasında platform=auto ise).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .orchestrate import run_analysis
from .pipelines import load_samplesheet, resolve_platforms


def _find_profile(profiles_dir: Path, sample: str) -> str | None:
    for ext in (".bracken", ".tsv", ".txt"):
        p = profiles_dir / f"{sample}{ext}"
        if p.exists():
            return str(p)
    return None


def cmd_run(args: argparse.Namespace) -> int:
    samples = load_samplesheet(args.samplesheet)
    resolve_platforms(samples, detect=not args.no_detect)

    if not args.profiles_dir and not args.kraken_db:
        print("HATA: --profiles-dir veya --kraken-db verilmeli.", file=sys.stderr)
        return 2
    profile_files = {}
    if args.profiles_dir:
        profiles_dir = Path(args.profiles_dir)
        for s in samples:
            pf = _find_profile(profiles_dir, s.sample)
            if pf:
                profile_files[s.sample] = pf
    if args.kraken_db:
        missing_reads = [s.sample for s in samples
                         if s.sample not in profile_files and not s.has_reads()]
        if missing_reads:
            print(f"UYARI: okuma dosyası olmayan örnekler atlandı: {missing_reads}", file=sys.stderr)
        samples = [s for s in samples if s.sample in profile_files or s.has_reads()]
        missing_files = [p for s in samples if s.sample not in profile_files
                         for p in (s.fastq_1, s.fastq_2, s.long_reads) if p and not Path(p).exists()]
        if missing_files:
            print(f"HATA: okuma dosyaları bulunamadı: {missing_files}", file=sys.stderr)
            return 2
    else:
        missing = [s.sample for s in samples if s.sample not in profile_files]
        if missing:
            print(f"UYARI: profil bulunamayan örnekler atlandı: {missing}", file=sys.stderr)
        samples = [s for s in samples if s.sample in profile_files]
    if not samples:
        print("HATA: analiz edilecek örnek kalmadı.", file=sys.stderr)
        return 2

    design = None if args.design == "auto" else args.design
    try:
        result = run_analysis(
            samples=samples,
            profile_files=profile_files,
            outdir=args.outdir,
            workdir=args.workdir or (Path(args.outdir).parent / ".microbiomeforge_work"),
            checkm2_file=args.checkm2,
            amr_file=args.amr,
            project=args.project,
            design=design,
            kraken_db=args.kraken_db,
            threads=args.threads,
            min_abundance=args.min_abundance,
            min_group_abundance=args.min_group_abundance,
            bracken_threshold=args.bracken_threshold,
            memory_mapping=args.memory_mapping,
        )
    except ValueError as e:
        print(f"HATA: {e}", file=sys.stderr)
        return 2
    print(f"Tamamlandı: {result['outdir']}")
    print(f"  Örnek: {result['n_samples']} · tasarım={result['design']} · "
          f"PERMANOVA p={result['permanova_p']} · "
          f"anlamlı takson={result['n_significant_taxa']} · kaynak={result['n_references']}")
    rep = result["report"].get("pdf") or result["report"].get("markdown")
    print(f"  Takson: {result['n_taxa_raw']} → filtre sonrası {result['n_taxa']}")
    print(f"  Rapor: {rep}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="microbiomeforge",
                                description="Karşılaştırmalı, literatüre atıflı mikrobiyom raporu.")
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="Analiz + rapor üret")
    r.add_argument("--samplesheet", required=True)
    r.add_argument("--profiles-dir", default=None,
                   help="Hazır Bracken profilleri (<örnek>.bracken), ör. nf-core çıktısı")
    r.add_argument("--kraken-db", default=None,
                   help="Kraken2 veritabanı dizini: verilirse FASTQ'dan yerel Kraken2→Bracken çalışır")
    r.add_argument("--threads", type=int, default=8)
    r.add_argument("--min-abundance", type=float, default=0.0001,
                   help="Hiçbir örnekte bu bağıl bolluğa ulaşmayan taksonları ele (varsayılan 0.0001 = %%0,01; 0=kapalı)")
    r.add_argument("--min-group-abundance", type=float, default=0.005,
                   help="Grup içi eşik: kendi grubunda (arke, virüs...) bu payı geçen tür tutulur (varsayılan 0.005)")
    r.add_argument("--bracken-threshold", type=int, default=10,
                   help="Bracken -t: taksonun tahmine girmesi için en az okuma (varsayılan 10)")
    r.add_argument("--memory-mapping", action="store_true",
                   help="Kraken2 DB'yi RAM'e yüklemeden diskten oku (yavaş, az bellek)")
    r.add_argument("--outdir", required=True)
    r.add_argument("--workdir", default=None)
    r.add_argument("--checkm2", default=None, help="CheckM2 quality_report.tsv")
    r.add_argument("--amr", default=None, help="funcscan/AMR tablosu")
    r.add_argument("--project", default="MicrobiomeForge")
    r.add_argument("--no-detect", action="store_true", help="Platform otomatik tespitini kapat")
    r.add_argument("--design", choices=["auto", "single", "comparative"], default="auto",
                   help="Çalışma tasarımı: auto=group metadata'sından çıkar (varsayılan), "
                        "single=tekli/betimsel, comparative=karşılaştırmalı")
    r.set_defaults(func=cmd_run)
    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
