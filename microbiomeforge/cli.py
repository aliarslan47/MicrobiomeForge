"""MicrobiomeForge komut satırı arayüzü.

    microbiomeforge run \
        --samplesheet config/samplesheet.csv \
        --profiles-dir nf_out/profiles \
        --outdir results/ [--checkm2 quality_report.tsv] [--amr amr.tsv]

`--profiles-dir` içinde her örnek için `<sample>.bracken` (veya .tsv) profil dosyası
beklenir. Bu dosyalar nf-core/taxprofiler çıktısından gelir. Platform (kısa/uzun) ham
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

    profiles_dir = Path(args.profiles_dir)
    profile_files = {}
    for s in samples:
        pf = _find_profile(profiles_dir, s.sample)
        if pf:
            profile_files[s.sample] = pf
    missing = [s.sample for s in samples if s.sample not in profile_files]
    if missing:
        print(f"UYARI: profil bulunamayan örnekler atlandı: {missing}", file=sys.stderr)
    if not profile_files:
        print("HATA: hiçbir örnek için profil bulunamadı.", file=sys.stderr)
        return 2

    design = None if args.design == "auto" else args.design
    try:
        result = run_analysis(
            samples=[s for s in samples if s.sample in profile_files],
            profile_files=profile_files,
            outdir=args.outdir,
            workdir=args.workdir or (Path(args.outdir).parent / ".microbiomeforge_work"),
            checkm2_file=args.checkm2,
            amr_file=args.amr,
            project=args.project,
            design=design,
        )
    except ValueError as e:
        print(f"HATA: {e}", file=sys.stderr)
        return 2
    print(f"Tamamlandı: {result['outdir']}")
    print(f"  Örnek: {result['n_samples']} · tasarım={result['design']} · "
          f"PERMANOVA p={result['permanova_p']} · "
          f"anlamlı takson={result['n_significant_taxa']} · kaynak={result['n_references']}")
    rep = result["report"].get("pdf") or result["report"].get("markdown")
    print(f"  Rapor: {rep}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="microbiomeforge",
                                description="Karşılaştırmalı, literatüre atıflı mikrobiyom raporu.")
    sub = p.add_subparsers(dest="command", required=True)
    r = sub.add_parser("run", help="Analiz + rapor üret")
    r.add_argument("--samplesheet", required=True)
    r.add_argument("--profiles-dir", required=True)
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
