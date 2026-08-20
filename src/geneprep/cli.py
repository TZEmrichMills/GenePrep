"""Command-line interface for GenePrep."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .codons import gc_bounds, gc_content
from .io import read_input, write_fasta, write_optimised_xlsx
from .optimize import optimize, OptParams
from .qc import assess


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="geneprep",
        description="Codon-optimise protein sequences for E. coli expression and "
                    "GenScript synthesis. Reads FASTA or Excel (protein or DNA), "
                    "minimises GC while preserving codon diversity, breaks repeats, "
                    "and writes an optimised FASTA plus an Excel workbook with the "
                    "order-ready constructs and QC.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Excel input columns: name, protein_sequence, [5prime_flank], [3prime_flank]

Examples:
  geneprep designs.fasta
  geneprep designs.xlsx --flank-5 CTGAGC... --flank-3 TCCGGC...
  geneprep designs.xlsx --copy-flanks-from 1
  geneprep designs.fasta --exclude-sites BsaI=GGTCTC,EcoRI=GAATTC
  geneprep designs.fasta --candidates 100 --seed 1
""")
    p.add_argument("input", type=Path, help="Input .fasta or .xlsx")
    p.add_argument("-o", "--output", type=Path,
                   help="Output FASTA (default: <input>_optimised.fasta)")
    p.add_argument("--xlsx", type=Path,
                   help="Output Excel workbook (default: <input>_optimised.xlsx). "
                        "Contains an Optimised sheet and a QC sheet.")
    p.add_argument("--flank-5", default="", help="5' flanking DNA added to every gene")
    p.add_argument("--flank-3", default="", help="3' flanking DNA added to every gene")
    p.add_argument("--copy-flanks-from", type=int, metavar="N",
                   help="Use gene N's flanks (1-indexed) for all genes")
    p.add_argument("--exclude-sites", default="",
                   help="Restriction sites to remove, e.g. BsaI=GGTCTC,EcoRI=GAATTC "
                        "(default: none — opt in for your cloning strategy)")
    p.add_argument("--stop", choices=("auto", "always", "never"), default="auto",
                   help="How to handle the stop codon (default: auto). "
                        "auto = append TAA only if the AA sequence does not end "
                        "in '*' and the 3' flank has no in-frame stop in its "
                        "first codon; always = always append TAA; never = never "
                        "append, and FAIL any gene with no detectable stop.")
    p.add_argument("--candidates", type=int, default=40,
                   help="Number of weighted-random candidates to sample (default: 40)")
    p.add_argument("--gc-target", "--target-gc", type=float, default=0.52,
                   dest="gc_target", metavar="F",
                   help="GC fraction the sampler tilts toward (default: 0.52). "
                        "AA composition may push actual GC higher. Raising this "
                        "toward 0.55–0.58 can help expression of AT-biased "
                        "proteins; alanine-rich designs are already GC-rich, so "
                        "for those the default is usually right.")
    p.add_argument("--no-mfe", action="store_true",
                   help="Skip the 5' mRNA secondary-structure (ΔG) scoring. "
                        "The MFE factor uses seqfold; disable this to shave "
                        "~10-40 ms/gene at the cost of ignoring translation-"
                        "initiation folding.")
    p.add_argument("--seed", type=int, help="Random seed for reproducible output")
    p.add_argument("--yes", "-y", action="store_true",
                   help="Skip the column-mapping confirmation prompt for Excel input")
    p.add_argument("-v", "--verbose", action="store_true")
    return p.parse_args(argv)


STOP_CODONS = ("TAA", "TAG", "TGA")


def resolve_stop(protein: str, flank_3: str, mode: str) -> tuple[bool, str]:
    """Decide whether to append TAA and describe where the stop comes from.

    Returns (append_taa, source) where `source` is one of:
      'AA seq', '3\\' flank', 'auto-added', 'always-added', 'MISSING'.
    """
    aa_has_stop = protein.endswith("*")
    flank_has_stop = flank_3[:3].upper() in STOP_CODONS

    if mode == "always":
        return True, "always-added"
    if mode == "never":
        # 'never' governs whether we *auto-add* a stop; an explicit '*' in the
        # AA sequence is a direct user request and is still honoured.
        if aa_has_stop:
            return True, "AA seq"
        if flank_has_stop:
            return False, "3' flank"
        return False, "MISSING"
    # auto: the AA-terminal '*' is a request for a stop -- honour it by adding TAA
    if aa_has_stop:
        return True, "AA seq"
    if flank_has_stop:
        return False, "3' flank"
    return True, "auto-added"


def parse_sites(spec: str) -> dict[str, str]:
    sites = {}
    if not spec:
        return sites
    for pair in spec.split(","):
        name, sep, seq = pair.partition("=")
        if not sep or not seq.strip():
            raise SystemExit(f"Error: bad --exclude-sites entry {pair!r} (want Name=SEQUENCE)")
        sites[name.strip()] = seq.strip().upper()
    return sites


def main(argv=None) -> None:
    # Ensure Unicode (ΔG, ≤, etc.) prints cleanly on Windows terminals whose
    # legacy default encoding is cp1252. Safe no-op on macOS/Linux.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            try:
                stream.reconfigure(encoding="utf-8", errors="replace")
            except Exception:
                pass

    args = parse_args(argv)
    if not args.input.exists():
        raise SystemExit(f"Error: input not found: {args.input}")

    out_fasta = args.output or args.input.with_name(args.input.stem + "_optimised.fasta")
    out_xlsx = args.xlsx or args.input.with_name(args.input.stem + "_optimised.xlsx")
    sites = parse_sites(args.exclude_sites)

    print(f"Reading {args.input}...")
    entries = read_input(args.input, auto_yes=args.yes)
    if not entries:
        raise SystemExit("Error: no sequences found")
    print(f"  {len(entries)} sequence(s)")

    # Resolve flanks
    default_5, default_3 = args.flank_5.upper(), args.flank_3.upper()
    if args.copy_flanks_from is not None:
        i = args.copy_flanks_from - 1
        if not 0 <= i < len(entries):
            raise SystemExit(f"Error: --copy-flanks-from {args.copy_flanks_from} out of range (1-{len(entries)})")
        src = entries[i]
        default_5, default_3 = src.flank_5 or default_5, src.flank_3 or default_3
        print(f"  Using flanks from '{src.name}' for all genes")
    for e in entries:
        e.flank_5 = e.flank_5 or default_5
        e.flank_3 = e.flank_3 or default_3

    params_base = dict(n_candidates=args.candidates, target_gc=args.gc_target,
                       use_mfe=not args.no_mfe, seed=args.seed)

    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    out_seqs, qc_results, sheet_rows, errors = [], [], [], 0
    for idx, e in enumerate(entries, 1):
        origin = " (from DNA input)" if e.source_kind == "dna" else ""
        print(f"\n[{idx}/{len(entries)}] {e.name} ({len(e.protein_sequence.rstrip('*'))} aa){origin}")
        protein_body = e.protein_sequence.rstrip("*")
        bad = sorted(set(protein_body) - valid_aa)
        if bad:
            print(f"  SKIPPED: invalid residue(s) {bad} — expected the 20 standard amino acids")
            errors += 1
            continue

        append_taa, stop_source = resolve_stop(e.protein_sequence, e.flank_3, args.stop)
        result = optimize(protein_body, sites=sites, params=OptParams(**params_base))
        body = result.dna + ("TAA" if append_taa else "")
        full = e.flank_5 + body + e.flank_3
        out_seqs.append((e.name, full))

        floor, _ceil = gc_bounds(protein_body)
        qc = assess(e.name, full, floor, sites)
        qc.stop_source = stop_source
        if stop_source == "MISSING":
            qc.status = "FAIL"
            qc.notes.append("no stop codon (AA seq has no '*', 3' flank has no in-frame stop, --stop=never)")
        elif stop_source == "auto-added":
            print(f"  note: appended TAA (no stop in AA seq or 3' flank start)")
        qc_results.append(qc)
        sheet_rows.append({
            "name": e.name, "protein": e.protein_sequence,
            "flank_5": e.flank_5, "flank_3": e.flank_3,
            "optimised_dna": body, "full_dna": full,
            "stop_source": stop_source, "status": qc.status,
        })

        swaps = result.repeat_swaps + result.site_swaps
        mfe_txt = f"  5'ΔG={qc.five_prime_mfe}" if qc.five_prime_mfe is not None else ""
        print(f"  GC={qc.gc_pct}% (floor {qc.gc_floor_pct}%)  CAI={qc.cai}  "
              f"5'GC={qc.five_prime_gc}%{mfe_txt}  window≤{qc.window_gc_max}%  "
              f"longest-repeat={qc.longest_repeat}bp  [{qc.status}]")
        if swaps:
            print(f"  repair swaps: {result.site_swaps} site, {result.repeat_swaps} repeat/homopolymer")
        if qc.notes and (args.verbose or qc.status != "PASS"):
            for n in qc.notes:
                print(f"    note: {n}")

    write_fasta(out_seqs, out_fasta)
    write_optimised_xlsx(sheet_rows, qc_results, out_xlsx)
    print(f"\nWrote {out_fasta}")
    print(f"Wrote {out_xlsx}  (sheets: Optimised, QC)")

    n_pass = sum(1 for r in qc_results if r.status == "PASS")
    n_review = sum(1 for r in qc_results if r.status == "REVIEW")
    n_fail = sum(1 for r in qc_results if r.status == "FAIL")
    print(f"\nSummary: {n_pass} PASS, {n_review} REVIEW, {n_fail} FAIL"
          + (f", {errors} SKIPPED" if errors else "") + f" (of {len(entries)})")
    if n_fail:
        print("  FAIL = excluded restriction site or homopolymer >=6 present; review before ordering.")
    if n_review:
        print("  REVIEW = high local GC or a long repeat; usually still synthesizable, but check.")


if __name__ == "__main__":
    main()
