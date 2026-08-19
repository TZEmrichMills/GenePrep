"""Regression tests locking in GenScript-calibrated behaviour.

Run with `pytest`, or directly: `python tests/test_geneprep.py`.
"""

from geneprep.codons import translate, gc_content, cai
from geneprep.optimize import optimize, OptParams
from geneprep.qc import (restriction_hits, direct_repeats, max_homopolymer,
                         count_anti_sd, count_stalling_pairs, five_prime_mfe,
                         HAS_SEQFOLD)

DESIGNS = {
    "V1":  "HGYVEEGTVEQLAQAIAAVRAAHPDAAVLQVGRVFIVVAPTAAAHDAALAALEAEAAALGVKIVTLSAALAAADPALKAIWDAWLAATAALLAALAAAVAAGDAAAAAALAAQLAPALLATLRAVAAVRAAA",
    "S1":  "HGYVEAAEGDRSRVRVTAVGADGEETWTVEWDASREAAEAALAAAYAAAEAAAAALEAALGRPLTLAEAVALHRAALAAQ",
    "P6":  "HGYVEAKTATALTFRAATIDPATGRILTKTFTAASAAAAAAAAVDWHAAQTAANLAANAGKLAPAVVATIRTNLATKKAAVTAALAAALATAAVGDVIAINWGGEAMRAAVLALARAIAAALGVRADRITVIT",
    "HS1": "HGYVEERTETKVTAVVWQEVNGLRREVTVTADSIEALLAAAAAATRELLRAALAAAPDASLTAEQQAALTEGHALMVVGEILLQLGRHDEAIAAFRKALAIYEAALGPDHPAVAAALYLLGVALLAAGKKEEAAAAFREAVAIAPDAPWGAAARAALEEL",
}


def test_translation_roundtrip():
    for prot in DESIGNS.values():
        dna = optimize(prot, params=OptParams(seed=1)).dna
        assert translate(dna) == prot


def test_metrics_match_genscript_envelope():
    for name, prot in DESIGNS.items():
        dna = optimize(prot, params=OptParams(seed=1)).dna
        assert 0.45 <= gc_content(dna) <= 0.70, f"{name} GC {gc_content(dna):.2f}"
        assert cai(dna) >= 0.65, f"{name} CAI {cai(dna):.3f}"
        assert not [r for r in direct_repeats(dna, 15) if r[2] >= 15], f"{name} has long repeat"
        assert max_homopolymer(dna) < 6, f"{name} homopolymer {max_homopolymer(dna)}"


def test_restriction_site_exclusion():
    sites = {"BsaI": "GGTCTC", "EcoRI": "GAATTC", "BamHI": "GGATCC"}
    for prot in DESIGNS.values():
        dna = optimize(prot, sites=sites, params=OptParams(seed=1)).dna
        assert restriction_hits(dna, sites) == []


def test_no_anti_sd_or_stalling_pairs():
    for prot in DESIGNS.values():
        dna = optimize(prot, params=OptParams(seed=1)).dna
        assert count_anti_sd(dna) == 0
        assert count_stalling_pairs(dna) == 0


def test_five_prime_mfe_weak_when_available():
    if not HAS_SEQFOLD:
        return  # skip silently
    for prot in DESIGNS.values():
        dna = optimize(prot, params=OptParams(seed=1)).dna
        mfe = five_prime_mfe(dna)
        assert mfe is None or mfe > -12.0, f"5' ΔG too strong: {mfe}"


def test_determinism():
    # Same seed -> identical output (best-of-N over overlapping seed ranges may
    # legitimately pick the same optimum for different seeds, so we only assert
    # the reproducibility guarantee here).
    prot = DESIGNS["V1"]
    a = optimize(prot, params=OptParams(seed=5)).dna
    b = optimize(prot, params=OptParams(seed=5)).dna
    assert a == b


if __name__ == "__main__":
    for fn in [test_translation_roundtrip, test_metrics_match_genscript_envelope,
               test_restriction_site_exclusion, test_no_anti_sd_or_stalling_pairs,
               test_five_prime_mfe_weak_when_available, test_determinism]:
        fn()
        print(f"PASS {fn.__name__}")
    print("\nAll tests passed.")
