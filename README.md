# GenePrep

Codon-optimise protein (or DNA) sequences for **E. coli expression** and **GenScript synthesis**, locally and privately. Built for GC-rich de novo designs. Batch in FASTA or Excel; get back an order-ready FASTA and an Excel workbook.

## Quickstart

```bash
cd GenePrep
pip install -e .
geneprep designs.fasta
```

Two files come out: `designs_optimised.fasta` and `designs_optimised.xlsx`. The workbook has an **Optimised** sheet (name, full construct ready to paste into an order form, protein, flanks, gene body, stop-source, status) and a **QC** sheet (GC%, GC floor, CAI, GC window range, longest repeat, homopolymer, restriction sites, notes).

## Common flags

```bash
geneprep designs.fasta --flank-5 CTGAGC...ATGCC --flank-3 TCCGGC...GAGTTGG
geneprep designs.xlsx --copy-flanks-from 1            # reuse gene 1's flanks
geneprep designs.fasta --exclude-sites BsaI=GGTCTC,EcoRI=GAATTC
geneprep designs.fasta --seed 1 --candidates 100      # reproducible; deeper search
geneprep designs.xlsx --yes                           # skip column-mapping prompt
geneprep designs.fasta --gc-target 0.55               # tilt AT-biased proteins up
geneprep designs.fasta --no-mfe                       # skip seqfold 5' ΔG scoring
```

`geneprep --help` lists everything.

## Input flexibility

- **Format**: `.fasta` or `.xlsx`.
- **Type**: protein or DNA — DNA is translated (frame 1), then re-optimised.
- **Excel columns**: headers are autodetected (`Construct Name`, `AA Seq`, `5' Flank`, `dna_sequence`, etc. all work); the mapping is printed and confirmed interactively. A generic `sequence` column is classified by content.

## How it works

Weighted-random codon sampling (E. coli K-12 usage) tilted toward lower GC, drawn many times and scored — the best candidate is kept, then a targeted sweep clears excluded restriction sites, long repeats, and homopolymers. Preserving codon diversity is what breaks up repeats in alanine-rich sequences.

Factors considered:

- Overall and local (50 bp window) GC content.
- Codon adaptation (CAI) against E. coli K-12; classic problem codons (AGG, AGA, CGA, CTA, ATA) are excluded from the pool.
- Direct and reverse-complement repeats; homopolymers.
- 5′ translation-initiation region: low-GC ramp, N-terminal codon-2 preference (A/T-starting codons at position 2 for better expression), and — when `seqfold` is available — the minimum free energy (ΔG) of the first 40 nt (Kudla et al. showed this factor alone can span >250× in expression).
- Anti-Shine-Dalgarno motifs in the ORF (`AGGAGG` and near-variants that create unwanted internal ribosome binding).
- Tandem stalling-codon pairs (rare-codon pairs cause ribosome pauses beyond what individual rare codons do).
- User-excluded restriction sites.

The 5′ ΔG factor adds a hard dependency on `seqfold` (pure Python, pip-installed automatically). Disable it with `--no-mfe` to shave ~10–40 ms/gene.

Because every alanine codon is ≥67% GC, high GC is intrinsic to these proteins — GenePrep minimises and reports GC against the protein's theoretical floor rather than failing it.

## Status labels

- **PASS** — order it.
- **REVIEW** — a long (≥15 bp) repeat or a local GC window over 85%; usually still synthesisable, but check.
- **FAIL** — an excluded restriction site is present, or a homopolymer ≥6, or (with `--stop never`) no detectable stop.

## Stop codon handling

`--stop auto` (default) appends TAA only when the AA sequence doesn't end in `*` and the 3' flank doesn't begin with an in-frame stop. Alternatives: `--stop always`, `--stop never`. The workbook's "Source of stop codon" column tells you which route each gene took.

## Notes

- The protein is used as-is (no Met added); the start codon usually lives in your 5′ flank.
- Output FASTA is the full construct (flanks + gene body) — ready to paste into an order.
