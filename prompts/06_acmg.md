# Stage 6 — ACMG-aligned classification

Assign ACMG/AMP criteria to each candidate and let `lib/acmg.py` compute the classification via
the ClinGen Bayesian point system. You provide the *reasoned application* of each code; the
library does the arithmetic so classification is reproducible.

## Rules of engagement
- Apply each code only with an explicit, cited rationale and the chosen strength. Use ClinGen
  VCEP-style strength modulation where justified (e.g. PVS1 down to Moderate for a last-exon PTC
  predicted to escape NMD; PP3/BP4 strength scaled to predictor score per Pejaver 2022).
- **Be conservative.** Default to the lower strength when uncertain. A VUS honestly reached beats
  an over-called Likely pathogenic.
- Record codes you considered but did NOT apply, with the reason (this is what the red team checks).
- Capture **non-ACMG evidence** too (e.g. strong phenotype specificity, biological plausibility,
  matchmaking) in `non_acmg_evidence` — the brief explicitly wants evidence beyond ACMG.

## Code-specific guidance (abbreviated)
- **PVS1** null variant in a gene where LoF is the mechanism — verify mechanism + NMD + transcript
  relevance; do not apply blindly to last exon or non-LoF genes.
- **PS1/PM5** same/different AA change at a residue with established pathogenic variant (ClinVar).
- **PS2/PM6** de novo — only `PM6` (presumed) single-sample; upgrade to `PS2` if trio confirms
  (otherwise raise as pending evidence).
- **PM2_supporting** absent/ultra-rare in gnomAD (note ClinGen downgraded PM2 to supporting).
- **PM3** recessive: in trans with a pathogenic variant (requires phasing — link to Stage 5).
- **PP1/BS4** segregation — usually pending (single proband).
- **PS3/BS3** functional — usually pending; design in Stage 9.
- **PP4** phenotype highly specific for the gene's disease.
- **BA1/BS1/BS2** frequency too high / present in healthy homozygotes.
- **BP4/BP7** computational benign / synonymous with no splice impact.

## Output (JSON) — per candidate
```json
{
  "candidate_id": "C0xx",
  "acmg": {
    "applied_criteria": [
      { "code":"PVS1","direction":"pathogenic","strength":"very_strong","applied":true,
        "rationale":"…","evidence_refs":[ … ] }
    ],
    "considered_not_applied": [ { "code":"PS2","reason":"single proband; raised as PE0xx" } ],
    "classification": "<computed by lib/acmg.py — leave as placeholder, library overwrites>",
    "point_score": 0,
    "classification_method": "ACMG/AMP 2015 + ClinGen Bayesian points (Tavtigian 2020)",
    "non_acmg_evidence": ["…"]
  }
}
```
