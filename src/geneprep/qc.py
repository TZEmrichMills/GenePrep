"""Sequence QC detectors and per-gene QC result.

Thresholds are calibrated against real GenScript-accepted output for GC-rich
de novo designs (see README). High GC and short/RC repeats are *reported*, not
failed, because they are intrinsic to these proteins and GenScript accepts them.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .codons import gc_content, cai as _cai, STALLING_CODONS, CODON_TO_AA

try:
    from seqfold import dg as _seqfold_dg
    HAS_SEQFOLD = True
except Exception:  # pragma: no cover
    HAS_SEQFOLD = False

COMPLEMENT = str.maketrans("ATCG", "TAGC")

# --- Calibrated thresholds ------------------------------------------------
GC_WINDOW = 50            # bp, sliding window for local GC
GC_WINDOW_HIGH = 0.85     # local windows above this are flagged (GenScript hit ~0.86)
GC_WINDOW_LOW = 0.30      # local windows below this are flagged
HOMOPOLYMER_FLAG = 6      # runs >= this are a problem (GenScript tolerated 5)
REPEAT_FLAG = 15          # direct repeats >= this (bp) are flagged (GenScript max ~12)
RAMP_CODONS = 15          # length of the 5' translation-initiation ramp
MFE_WINDOW = 40           # nt of 5' region used for MFE calculation
MFE_WARN = -10.0          # kcal/mol; more negative than this is flagged (Kudla et al.)

# Anti-Shine-Dalgarno motifs the ribosome sees as internal binding sites.
# Perfect anti-SD is AGGAGG (complement of the 16S rRNA 3' tail CCUCCU).
# The variants below are the strongest matches; we penalise all of them.
ANTI_SD_MOTIFS: tuple[str, ...] = (
    "AGGAGG", "AAGGAG", "AGGAGA", "AGAGGA", "GGAGGA", "AAGGAGG", "AGGAGGA",
)


def revcomp(seq: str) -> str:
    return seq.upper().translate(COMPLEMENT)[::-1]


@dataclass
class QCResult:
    name: str
    length: int = 0
    gc_pct: float = 0.0
    gc_floor_pct: float = 0.0
    cai: float = 0.0
    window_gc_min: float = 0.0
    window_gc_max: float = 0.0
    windows_flagged: int = 0
    longest_repeat: int = 0
    repeats_flagged: int = 0
    max_homopolymer: int = 0
    homopolymers_flagged: int = 0
    restriction_sites: list[tuple[str, int]] = field(default_factory=list)
    five_prime_gc: float = 0.0
    five_prime_mfe: float | None = None
    stalling_pairs: int = 0
    anti_sd_hits: int = 0
    stop_source: str = ""
    status: str = "PASS"
    notes: list[str] = field(default_factory=list)

    @property
    def restriction_str(self) -> str:
        return "; ".join(f"{n}@{p}" for n, p in self.restriction_sites) or "None"

    @property
    def notes_str(self) -> str:
        return " | ".join(self.notes) or ""


def homopolymer_runs(seq: str, min_run: int) -> list[tuple[str, int, int]]:
    return [(m.group()[0], m.start(), len(m.group()))
            for m in re.finditer(r"(.)\1{" + str(min_run - 1) + r",}", seq.upper())]


def max_homopolymer(seq: str) -> int:
    runs = homopolymer_runs(seq, 2)
    return max((r[2] for r in runs), default=1)


def direct_repeats(seq: str, min_len: int) -> list[tuple[int, int, int]]:
    """Longest non-overlapping exact direct repeats >= min_len.

    Fast: index kmers at min_len, then extend each duplicated seed to its full
    length. Returns (start_of_first, start_of_second, length).
    """
    seq = seq.upper()
    n = len(seq)
    positions: dict[str, list[int]] = {}
    for i in range(n - min_len + 1):
        positions.setdefault(seq[i:i + min_len], []).append(i)

    found: list[tuple[int, int, int]] = []
    reported: set[int] = set()
    for kmer, pos in positions.items():
        if len(pos) < 2:
            continue
        a, b = pos[0], pos[1]
        if a in reported:
            continue
        length = min_len
        while (b + length < n and seq[a + length] == seq[b + length]
               and a + length < b):
            length += 1
        found.append((a, b, length))
        for p in range(a, a + length):
            reported.add(p)
    found.sort(key=lambda x: -x[2])
    return found


def longest_repeat(seq: str, probe: int = 10) -> int:
    reps = direct_repeats(seq, probe)
    return max((r[2] for r in reps), default=0)


def restriction_hits(seq: str, sites: dict[str, str]) -> list[tuple[str, int]]:
    seq = seq.upper()
    hits = []
    for name, site in sites.items():
        site = site.upper()
        for pattern in ({site, revcomp(site)}):
            start = 0
            while (pos := seq.find(pattern, start)) != -1:
                hits.append((name, pos))
                start = pos + 1
    hits.sort(key=lambda x: x[1])
    return hits


def count_stalling_pairs(dna: str) -> int:
    """Count adjacent codon pairs where both codons are known stallers."""
    n = 0
    prev_bad = False
    for i in range(0, len(dna) - 2, 3):
        c = dna[i:i + 3]
        bad = c in STALLING_CODONS
        if bad and prev_bad:
            n += 1
        prev_bad = bad
    return n


def count_anti_sd(dna: str) -> int:
    """Count anti-Shine-Dalgarno-like motifs anywhere in the sequence."""
    dna = dna.upper()
    reported: set[int] = set()
    for motif in ANTI_SD_MOTIFS:
        start = 0
        while (pos := dna.find(motif, start)) != -1:
            if pos not in reported:
                reported.add(pos)
            start = pos + 1
    return len(reported)


def five_prime_mfe(dna: str, window: int = MFE_WINDOW) -> float | None:
    """Minimum free energy (kcal/mol) of the first `window` nt using seqfold.
    Returns None if seqfold isn't installed. More positive = weaker structure
    = easier translation initiation."""
    if not HAS_SEQFOLD or len(dna) < 10:
        return None
    try:
        val = _seqfold_dg(dna[:window].replace("U", "T"))
        return float(val) if val is not None else None
    except Exception:
        return None


def window_gc_extremes(seq: str, window: int = GC_WINDOW,
                       high: float = GC_WINDOW_HIGH,
                       low: float = GC_WINDOW_LOW) -> tuple[float, float, int]:
    seq = seq.upper()
    if len(seq) < window:
        g = gc_content(seq)
        return g, g, 0
    lo = 1.0
    hi = 0.0
    flagged = 0
    for i in range(len(seq) - window + 1):
        g = gc_content(seq[i:i + window])
        lo = min(lo, g)
        hi = max(hi, g)
        if g > high or g < low:
            flagged += 1
    return lo, hi, flagged


def assess(name: str, dna: str, gc_floor: float,
           restriction_sites: dict[str, str]) -> QCResult:
    """Run all detectors on a coding sequence (already including flanks if any)."""
    r = QCResult(name=name)
    r.length = len(dna)
    r.gc_pct = round(gc_content(dna) * 100, 1)
    r.gc_floor_pct = round(gc_floor * 100, 1)
    r.cai = round(_cai(dna), 3)

    lo, hi, flagged = window_gc_extremes(dna)
    r.window_gc_min = round(lo * 100, 1)
    r.window_gc_max = round(hi * 100, 1)
    r.windows_flagged = flagged

    r.longest_repeat = longest_repeat(dna)
    r.repeats_flagged = sum(1 for rep in direct_repeats(dna, REPEAT_FLAG))
    r.max_homopolymer = max_homopolymer(dna)
    r.homopolymers_flagged = len(homopolymer_runs(dna, HOMOPOLYMER_FLAG))
    r.restriction_sites = restriction_hits(dna, restriction_sites)
    r.five_prime_gc = round(gc_content(dna[:RAMP_CODONS * 3]) * 100, 1)
    r.five_prime_mfe = five_prime_mfe(dna)
    if r.five_prime_mfe is not None:
        r.five_prime_mfe = round(r.five_prime_mfe, 2)
    r.stalling_pairs = count_stalling_pairs(dna)
    r.anti_sd_hits = count_anti_sd(dna)

    # --- Status + notes ---------------------------------------------------
    if r.restriction_sites:
        r.status = "FAIL"
        r.notes.append(f"excluded restriction site(s) present: {r.restriction_str}")
    if r.homopolymers_flagged:
        r.status = "FAIL"
        r.notes.append(f"homopolymer run >= {HOMOPOLYMER_FLAG}")
    if r.repeats_flagged:
        if r.status != "FAIL":
            r.status = "REVIEW"
        r.notes.append(f"{r.repeats_flagged} direct repeat(s) >= {REPEAT_FLAG} bp (longest {r.longest_repeat})")
    if r.window_gc_max > GC_WINDOW_HIGH * 100:
        if r.status == "PASS":
            r.status = "REVIEW"
        r.notes.append(f"local GC window peaks at {r.window_gc_max}%")
    if r.five_prime_mfe is not None and r.five_prime_mfe < MFE_WARN:
        if r.status == "PASS":
            r.status = "REVIEW"
        r.notes.append(f"strong 5' RNA structure (ΔG {r.five_prime_mfe} kcal/mol)")
    if r.stalling_pairs:
        r.notes.append(f"{r.stalling_pairs} tandem stalling-codon pair(s)")
    if r.anti_sd_hits:
        r.notes.append(f"{r.anti_sd_hits} anti-SD motif(s) in ORF")
    # GC is reported honestly relative to the intrinsic floor, never failed.
    if r.gc_pct - r.gc_floor_pct <= 3.0 and r.gc_pct > 60:
        r.notes.append(f"GC near intrinsic floor ({r.gc_floor_pct}%); residual GC is protein-inherent")

    return r
