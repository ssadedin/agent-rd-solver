# Stage 8 — Synthesis & ranking

Assemble the final candidate set from everything that survived Stage 7 and order it.

## Tasks
1. Include only candidates with `survives_review: true`. Verify the excluded ones each have a
   recorded reason in `excluded_variants_summary`.
2. Rank by, in order: adversarial survival/confidence → ACMG classification → phenotype-match
   strength → ACMG point score (per `config/thresholds.yaml:reporting.rank_by`). Break ties with
   gene–disease validity and biological plausibility.
3. Mark at most one candidate `is_primary_hypothesis: true` per distinct phenotype explanation;
   note when two candidates jointly explain the phenotype (e.g. comp-het, or dual diagnosis).
4. Write each candidate's `overall_assessment.summary`: a clinician-readable paragraph stating
   what the variant is, why it is/ isn't likely causal, the single biggest remaining uncertainty,
   and what would resolve it.
5. Produce `run_summary` with counts and a `diagnostic_conclusion`:
   - `likely_solved` — ≥1 P/LP candidate, survives review, strong phenotype fit.
   - `promising_candidates` — strong VUS/LP needing one piece of evidence.
   - `candidates_pending_evidence` — best candidates blocked on pending items.
   - `uninformative` — nothing credible; recommend reanalysis path.

## Output (JSON)
A complete report object conforming to `schema/candidate_report.schema.json`
(`candidate_variants`, `excluded_variants_summary`, `run_summary`), merging all prior stage
outputs into each `candidate`. Do not introduce any new factual claim here — synthesis only.
