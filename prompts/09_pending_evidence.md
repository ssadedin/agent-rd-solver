# Stage 9 — Pending evidence & critically-reviewed experimental plans

Produce the second required deliverable: a prioritised list of **specific** questions or
experiments that, if answered, would materially change the classification or ranking of one or
more candidates. This is where the workflow's honesty about its own limits lives.

Consolidate every `pending_evidence` seed raised in Stages 0, 6, and 7, plus any new ones, into a
deduplicated, prioritised list.

## Quality bar for every item
- It must be **decisive**: name the exact ACMG criteria it would flip and the classification
  change it could produce (`could_change_classification_from` → `to`). If it cannot change a
  decision, drop it.
- It must be **actionable** by someone outside this system (clinician, lab, family), not
  something the agent could have done with the data at hand.
- `priority` reflects (decisiveness × number/importance of candidates affected).

## Categories
`clinical_data`, `segregation`, `phenotyping`, `imaging`, `additional_sequencing`,
`functional_experiment`, `reanalysis`.

## For `functional_experiment` items — a detailed, falsifiable plan is mandatory
Provide `experimental_plan` with: `method`, `objective`, `hypothesis` (stated so a result can
**refute** it), `protocol_outline` (ordered steps), `sample_requirements`, `controls`
(positive + negative + carrier as appropriate), `readout`, `interpretation_framework`
(precisely which result supports vs refutes pathogenicity, with the threshold/effect size that
counts), `limitations`, `estimated_turnaround`.

Pick the method to the variant's mechanism, e.g.:
- splice candidate → patient RNA (blood/fibroblast) RT-PCR ± targeted RNA-seq, or minigene if the
  tissue doesn't express the transcript;
- missense in enzyme → enzyme activity assay in patient cells;
- LoF / dosage → Western/qPCR for protein/transcript level, or allele-specific expression;
- UPD/imprinting → methylation-specific assay;
- candidate de novo → trio Sanger;
- uncertain SV breakpoint → long-read sequencing / breakpoint-spanning PCR.

## Then RED-TEAM the plan (`plan_review`)
Critique each plan as a sceptical reviewer would a grant: wrong/under-expressing tissue, missing
controls, confounders (e.g. NMD masking a splice effect; common benign isoforms), inability to
distinguish the hypothesis from alternatives, turnaround vs clinical urgency, sample feasibility.
Grade critiques `fatal|major|minor`, give a `resolution`, and a final
`verdict: recommended | recommended_with_changes | not_recommended`. Do not recommend a plan
that cannot, even in principle, change a classification.

## Output (JSON)
An array of `pending_evidence` objects conforming to the schema (`experimental_plan` +
`plan_review` populated for functional experiments).
