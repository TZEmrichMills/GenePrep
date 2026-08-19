"""E. coli K-12 codon usage, synonymous maps, and GC-aware reverse translation.

Everything is derived from a single frequency table so there is no external
codon-table dependency. Frequencies are per-thousand (Kazusa / GenBank).
"""

from __future__ import annotations

import math
import random
from typing import Optional

ECOLI_CODON_FREQ: dict[str, dict[str, float]] = {
    "F": {"TTT": 22.0, "TTC": 16.2},
    "L": {"TTA": 13.7, "TTG": 13.3, "CTT": 11.0, "CTC": 10.9, "CTA": 3.9, "CTG": 52.6},
    "I": {"ATT": 30.3, "ATC": 24.8, "ATA": 4.3},
    "M": {"ATG": 27.8},
    "V": {"GTT": 18.3, "GTC": 15.2, "GTA": 10.8, "GTG": 26.0},
    "S": {"TCT": 8.5, "TCC": 8.5, "TCA": 7.2, "TCG": 8.9, "AGT": 8.8, "AGC": 16.0},
    "P": {"CCT": 7.0, "CCC": 5.5, "CCA": 8.4, "CCG": 23.2},
    "T": {"ACT": 9.0, "ACC": 23.4, "ACA": 7.1, "ACG": 14.4},
    "A": {"GCT": 15.3, "GCC": 25.5, "GCA": 20.0, "GCG": 33.7},
    "Y": {"TAT": 16.2, "TAC": 12.1},
    "H": {"CAT": 12.8, "CAC": 9.7},
    "Q": {"CAA": 15.2, "CAG": 28.8},
    "N": {"AAT": 17.7, "AAC": 21.5},
    "K": {"AAA": 33.6, "AAG": 10.3},
    "D": {"GAT": 32.2, "GAC": 19.1},
    "E": {"GAA": 39.6, "GAG": 17.8},
    "C": {"TGT": 5.2, "TGC": 6.5},
    "W": {"TGG": 15.2},
    "R": {"CGT": 20.9, "CGC": 22.0, "CGA": 3.6, "CGG": 5.4, "AGA": 2.1, "AGG": 1.2},
    "G": {"GGT": 24.7, "GGC": 29.4, "GGA": 8.0, "GGG": 11.1},
    "*": {"TAA": 2.0, "TAG": 0.2, "TGA": 1.0},
}

# Codons whose relative adaptiveness (freq / max-for-that-aa) is below this are
# excluded from seeding to avoid translationally poor codons. Kept mild so we do
# not over-constrain: excludes CTA, AGA, AGG (the classic E. coli rare codons).
RARE_REL_ADAPT = 0.10

# --- Derived lookup tables ------------------------------------------------
CODON_TO_AA: dict[str, str] = {}
SYNONYMS: dict[str, list[str]] = {}          # aa -> codons
REL_ADAPT: dict[str, float] = {}             # codon -> freq / max-for-aa
for _aa, _codons in ECOLI_CODON_FREQ.items():
    SYNONYMS[_aa] = list(_codons)
    _mx = max(_codons.values())
    for _c, _f in _codons.items():
        CODON_TO_AA[_c] = _aa
        REL_ADAPT[_c] = _f / _mx

# aa -> list of usable (non-rare) codons, GC of each
USABLE_CODONS: dict[str, list[str]] = {
    aa: [c for c in codons if REL_ADAPT[c] >= RARE_REL_ADAPT or len(codons) == 1]
    for aa, codons in SYNONYMS.items()
}

GC_BASES = frozenset("GC")


def codon_gc(codon: str) -> float:
    return sum(1 for b in codon if b in GC_BASES) / 3.0


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    s = seq.upper()
    return (s.count("G") + s.count("C")) / len(s)


def synonymous_codons(codon: str, include_rare: bool = False) -> list[str]:
    """All codons coding the same amino acid as `codon`."""
    aa = CODON_TO_AA.get(codon.upper())
    if aa is None:
        return [codon]
    return SYNONYMS[aa] if include_rare else USABLE_CODONS[aa]


def gc_bounds(protein: str) -> tuple[float, float]:
    """Theoretical min and max GC achievable for a protein, using the lowest-
    and highest-GC synonymous codon at each position. Lets us report how close
    a design is to its intrinsic GC floor (some GC is unavoidable)."""
    lo = hi = 0
    n = 0
    for aa in protein.upper():
        if aa in ("*", "X") or aa not in SYNONYMS:
            continue
        gcs = [codon_gc(c) for c in USABLE_CODONS[aa]]
        lo += min(gcs)
        hi += max(gcs)
        n += 1
    if n == 0:
        return 0.0, 0.0
    return lo / n, hi / n


def _weighted_choice(codons: list[str], weights: list[float],
                     rng: random.Random) -> str:
    return rng.choices(codons, weights=weights, k=1)[0]


def seed_translation(protein: str, target_gc: float = 0.55,
                     ramp_codons: int = 15, ramp_gc: float = 0.45,
                     seed: Optional[int] = None) -> str:
    """Produce an initial DNA sequence by GC-aware weighted sampling.

    Codons are drawn proportional to E. coli frequency, tilted toward AT-rich
    or GC-rich synonyms to steer the running GC toward `target_gc`. The first
    `ramp_codons` are pulled toward the lower `ramp_gc` to build a low-GC
    translation-initiation ramp (mirrors what GenScript does at the 5' end).
    """
    rng = random.Random(seed)
    out: list[str] = []
    running_gc_sum = 0.0
    running_len = 0

    for i, aa in enumerate(protein.upper()):
        if aa == "*":
            break
        if aa not in USABLE_CODONS:
            raise ValueError(f"Unknown amino acid {aa!r} at position {i + 1}")

        want = ramp_gc if i < ramp_codons else target_gc
        # How far current GC is from where we want it; steer next codon.
        current = running_gc_sum / running_len if running_len else want
        pull = want - current  # positive -> need more GC

        codons = USABLE_CODONS[aa]
        weights = []
        for c in codons:
            base = ECOLI_CODON_FREQ[aa][c]
            # tilt: codons on the needed side of 0.5 GC get boosted
            tilt = math.exp(4.0 * pull * (codon_gc(c) - 0.5))
            weights.append(base * tilt)
        choice = _weighted_choice(codons, weights, rng)
        out.append(choice)
        running_gc_sum += codon_gc(choice) * 3
        running_len += 3

    return "".join(out)


def cai(dna: str) -> float:
    """Codon Adaptation Index against E. coli K-12."""
    log_sum = 0.0
    n = 0
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        if codon in ("TAA", "TAG", "TGA"):
            continue
        w = REL_ADAPT.get(codon)
        if w and w > 0:
            log_sum += math.log(w)
            n += 1
    return math.exp(log_sum / n) if n else 0.0


def translate(dna: str) -> str:
    """Translate DNA to protein using the standard table (stop -> '*')."""
    out = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i + 3].upper()
        out.append(CODON_TO_AA.get(codon, "X"))
    return "".join(out)
