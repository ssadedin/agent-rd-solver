# Stage 0 — Intake & phenotype structuring

You are a clinical geneticist abstracting a case for downstream genomic analysis.

## Inputs provided
- `phenotype.hpo` — patient HPO terms (ids ± labels).
- `clinical_notes.txt` — free-text clinical summary.

## Tasks
1. Normalise every HPO term: confirm the id is current (Ensembl/HPO API or local `hp.obo`),
   attach the canonical label, and note any obsolete/replaced terms.
2. Extract structured clinical facts from the free text, each with a quoted source span:
   - **age of onset** (map to HPO onset term where possible)
   - **inheritance clues** (affected relatives, pedigree pattern, consanguinity, parental ages)
   - **prior testing** and, critically, **what each result rules in or out** for this analysis
     (e.g. "normal karyotype → large balanced rearrangement unlikely but cryptic still possible";
     "negative gene-panel for X → still need to check panel-version coverage and CNVs").
   - ethnicity / ancestry (informs which gnomAD popmax to weight)
   - cardinal/"handle" phenotypes (the most specific, diagnostically discriminating features)
   - any features that are explicitly ABSENT (pertinent negatives).
3. Derive a ranked list of **inheritance hypotheses** consistent with the pedigree and onset.
4. Flag contradictions or ambiguities in the notes as `pending_evidence` (category
   `clinical_data`) rather than resolving them by assumption.

## Output (JSON)
Conform to the `patient` object in `schema/candidate_report.schema.json`, plus a
`derived_search_strategy` block:
```json
{
  "patient": { ... schema patient object ... },
  "derived_search_strategy": {
    "cardinal_phenotypes": ["HP:..."],
    "pertinent_negatives": ["HP:..."],
    "inheritance_priority": ["autosomal_recessive", "de_novo", ...],
    "regions_of_special_interest": ["genes/loci implicated by prior testing or notes"],
    "prior_testing_gaps": ["what earlier tests could have missed and we must re-examine"]
  },
  "pending_evidence_seeds": [ { "category": "...", "question": "...", "rationale": "..." } ]
}
```
