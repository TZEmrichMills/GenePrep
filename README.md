# GenePrep

A quick tool that turns protein sequences into codon-optimised, order-ready gene sequences for **E. coli expression** and **GenScript synthesis**. Runs locally, keeps your designs private, and hands back an Excel file where every row has a full construct you can paste straight into the order form.

Built for GC-rich de novo designs (alanine-rich helical bundles), but works fine for any protein.

## Quickstart

You'll need [Python 3.9 or newer](https://www.python.org/downloads/). Works on macOS, Windows, and Linux.

```bash
git clone https://github.com/TZEmrichMills/GenePrep.git
cd GenePrep
pip install -e .
```

Then feed it a FASTA or Excel file:

```bash
geneprep designs.fasta
```

That's it. Two files come out next to the input:

- **`designs_optimised.fasta`** — the optimised DNA sequences.
- **`designs_optimised.xlsx`** — an Excel workbook with two sheets:
  - **Optimised**: name, then the full construct (paste this into GenScript), then protein, flanks, gene body, and a PASS/REVIEW/FAIL status.
  - **QC**: every quality metric — GC%, CAI, repeats, 5′ folding energy, and more.

### Adding your flanks

If every construct uses the same 5′ and 3′ flanks (a signal peptide, a linker, a tag):

```bash
geneprep designs.fasta --flank-5 CTGAGC...ATGCC --flank-3 TCCGGC...GAGTTGG
```

Or, if you keep per-gene flanks in your Excel sheet, one gene's flanks can be reused for all:

```bash
geneprep designs.xlsx --copy-flanks-from 1
```

## More flags

```bash
geneprep designs.fasta --exclude-sites BsaI=GGTCTC,EcoRI=GAATTC
geneprep designs.fasta --seed 1 --candidates 100      # reproducible; deeper search
geneprep designs.xlsx --yes                           # skip the column-mapping prompt
geneprep designs.fasta --gc-target 0.55               # nudge AT-biased proteins up
geneprep designs.fasta --no-mfe                       # skip the seqfold 5′ folding scoring
```

`geneprep --help` lists everything.

## Input flexibility

- **Format**: `.fasta` or `.xlsx`.
- **Type**: protein or DNA — DNA gets translated (frame 1) and re-optimised.
- **Excel headers**: autodetected. `Construct Name`, `AA Seq`, `5' Flank`, `dna_sequence` — all fine. The detected mapping is printed and confirmed before it runs. A generic `sequence` column is classified by content.

## How it works

GenePrep samples many codon-choice candidates using weighted-random selection (tilted toward E. coli's preferred codons *and* toward lower GC), scores them, and keeps the best. A short repair pass then clears any excluded restriction sites, long repeats, or homopolymers. Preserving **codon diversity** helps breaks up repeats in alanine-rich sequences.

What it scores against:

- Overall and local (50 bp window) GC content.
- Codon adaptation (CAI) against E. coli K-12; classic problem codons (AGG, AGA, CGA, CTA, ATA) are excluded from the pool.
- Direct and reverse-complement repeats; homopolymers.
- 5′ translation-initiation region: low-GC ramp, an A/T-starting codon preference at position 2, and the minimum free energy (ΔG) of the first 40 nt via `seqfold` (Kudla et al. showed this factor alone can span >250× in expression). Disable with `--no-mfe`.
- Anti-Shine-Dalgarno motifs in the ORF (`AGGAGG` and near-variants that create unwanted internal ribosome binding).
- Tandem stalling-codon pairs.
- Any restriction sites you asked to exclude.

Because every alanine codon is ≥67% GC, high GC is intrinsic to alanine-rich proteins — GenePrep minimises and *reports* GC against the protein's theoretical floor rather than failing it.

## Status labels

- **PASS** — order it.
- **REVIEW** — has a long (≥15 bp) repeat or a local GC window over 85%; usually still synthesisable, but worth a look.
- **FAIL** — an excluded restriction site slipped through, a homopolymer ≥6 is present, or (with `--stop never`) no stop was detected.

## Stop codon handling

`--stop auto` (default) appends TAA only when the AA sequence doesn't end in `*` and the 3′ flank doesn't begin with an in-frame stop. Alternatives: `--stop always`, `--stop never`. The workbook's "Source of stop codon" column tells you which route each gene took.

## Notes

- The protein is used as-is (no Met added); the start codon usually lives in your 5′ flank.
- Output FASTA is the full construct (flanks + gene body) — ready to paste into an order.
- If `pip install -e .` gives you permissions trouble, use `pip install --user -e .`, or set up a virtual environment first (`python -m venv .venv && source .venv/bin/activate` on macOS/Linux, `.venv\Scripts\activate` on Windows).
