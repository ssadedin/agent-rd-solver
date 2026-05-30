"""
ACMG/AMP classification via the ClinGen Bayesian point system (Tavtigian et al. 2020).

The reasoning agent (Stage 6) decides WHICH criteria apply and at WHAT strength, with cited
rationale. This module does the deterministic arithmetic so the final class is reproducible and
cannot be fudged by prose. It also computes the maximum credible allele frequency
(Whiffin et al. 2017) used by the frequency filter.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Iterable

# Point values per strength, signed by direction. (ClinGen / Tavtigian 2020)
_STRENGTH_POINTS = {
    "very_strong": 8,
    "strong": 4,
    "moderate": 2,
    "supporting": 1,
    "stand_alone": 8,   # BA1 stand-alone benign behaves as -8 (auto-benign, handled below)
}

# Recognised codes and their natural direction (a sanity check on agent output).
PATHOGENIC_CODES = {"PVS1", "PS1", "PS2", "PS3", "PS4",
                    "PM1", "PM2", "PM3", "PM4", "PM5", "PM6",
                    "PP1", "PP2", "PP3", "PP4", "PP5"}
BENIGN_CODES = {"BA1", "BS1", "BS2", "BS3", "BS4", "BP1", "BP2", "BP3", "BP4", "BP5", "BP6", "BP7"}


@dataclass
class Criterion:
    code: str
    direction: str          # "pathogenic" | "benign"
    strength: str           # see _STRENGTH_POINTS
    applied: bool = True

    def points(self) -> int:
        if not self.applied:
            return 0
        mag = _STRENGTH_POINTS[self.strength]
        return mag if self.direction == "pathogenic" else -mag


def classify(criteria: Iterable[dict]) -> dict:
    """
    criteria: list of dicts with keys code, direction, strength, applied.
    Returns {classification, point_score, applied_codes, warnings}.

    Point thresholds (Tavtigian 2020):
        >= 10            Pathogenic
        6 .. 9           Likely pathogenic
        0 .. 5           Uncertain significance
        -1 .. -6         Likely benign
        <= -7            Benign
    BA1 (stand_alone benign) forces Benign regardless of other codes.
    """
    warnings: list[str] = []
    crits = [Criterion(c["code"], c["direction"], c["strength"], c.get("applied", True))
             for c in criteria]

    # Direction sanity check against the canonical sign of each code.
    for c in crits:
        if c.code in PATHOGENIC_CODES and c.direction != "pathogenic":
            warnings.append(f"{c.code} usually pathogenic but marked {c.direction}")
        if c.code in BENIGN_CODES and c.direction != "benign":
            warnings.append(f"{c.code} usually benign but marked {c.direction}")
        if c.code not in PATHOGENIC_CODES | BENIGN_CODES:
            warnings.append(f"unrecognised code {c.code}")

    # BA1 stand-alone override.
    if any(c.code == "BA1" and c.applied for c in crits):
        return {"classification": "Benign", "point_score": -8,
                "applied_codes": [c.code for c in crits if c.applied],
                "warnings": warnings, "override": "BA1 stand-alone"}

    score = sum(c.points() for c in crits)

    if score >= 10:
        cls = "Pathogenic"
    elif score >= 6:
        cls = "Likely pathogenic"
    elif score >= 0:
        cls = "Uncertain significance"
    elif score >= -6:
        cls = "Likely benign"
    else:
        cls = "Benign"

    return {"classification": cls, "point_score": score,
            "applied_codes": [c.code for c in crits if c.applied],
            "warnings": warnings, "override": None}


def max_credible_af(prevalence: float, max_allelic_contribution: float,
                    max_genotype_freq: float, penetrance: float = 1.0) -> float:
    """
    Whiffin et al. 2017 maximum credible population AF for a disease allele.

        max_AF = (prevalence * max_allelic_contribution * (1/penetrance)) / (2 * ...)

    Simplified, widely-used form for a dominant disorder:
        max_AF = (prevalence * max_allelic_contribution) / (2 * penetrance)
    For recessive, the genotype frequency (prevalence) maps to allele freq via sqrt.

    Args are kept explicit so the agent must supply justified inputs (else falls back to the
    coarse thresholds in config/thresholds.yaml).

    Returns an allele-frequency ceiling above which a variant is too common to be causal.
    """
    if prevalence <= 0 or penetrance <= 0:
        raise ValueError("prevalence and penetrance must be > 0")
    # Dominant approximation.
    af = (prevalence * max_allelic_contribution) / (2.0 * penetrance)
    # Never exceed the supplied per-disease max genotype frequency-derived ceiling.
    return min(af, max_genotype_freq)


def pp3_bp4_strength_from_score(tool: str, score: float) -> tuple[str, str] | None:
    """
    Map a continuous in-silico score to a calibrated PP3/BP4 strength
    (Pejaver et al. 2022 calibration, abbreviated thresholds).
    Returns (direction, strength) or None if indeterminate.
    """
    table = {
        # tool: [(threshold, direction, strength), ...] evaluated in order
        "revel": [(0.932, "pathogenic", "strong"), (0.773, "pathogenic", "moderate"),
                  (0.644, "pathogenic", "supporting"), (0.290, "benign", "supporting"),
                  (0.183, "benign", "moderate"), (0.016, "benign", "strong")],
        "alphamissense": [(0.99, "pathogenic", "strong"), (0.972, "pathogenic", "moderate"),
                          (0.564, "pathogenic", "supporting"), (0.34, "benign", "supporting")],
        "spliceai": [(0.5, "pathogenic", "moderate"), (0.2, "pathogenic", "supporting"),
                     (0.1, "benign", "supporting")],
    }
    rules = table.get(tool.lower())
    if rules is None or score is None:
        return None
    for thr, direction, strength in rules:
        if direction == "pathogenic" and score >= thr:
            return (direction, strength)
        if direction == "benign" and score <= thr:
            return (direction, strength)
    return None


if __name__ == "__main__":  # tiny self-test
    demo = [
        {"code": "PVS1", "direction": "pathogenic", "strength": "very_strong"},
        {"code": "PM2", "direction": "pathogenic", "strength": "supporting"},
        {"code": "PP3", "direction": "pathogenic", "strength": "moderate"},
    ]
    print(classify(demo))   # -> Likely pathogenic / Pathogenic depending on points (8+1+2=11 -> Pathogenic)
