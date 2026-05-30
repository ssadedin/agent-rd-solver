# Stage 4 — Evidence gathering (ClinVar / PubMed / Ensembl / UCSC / gnomAD)

For each candidate from Stages 2–3, assemble the external evidence dossier. Be exhaustive and
precise; this dossier feeds both ACMG classification and the adversarial review, so gather
**disconfirming** evidence as eagerly as confirming evidence.

## Per-candidate retrieval checklist
1. **ClinVar** (`lib/api_clients.py:clinvar_*`):
   - Exact variant: classification, review status (star rating), conditions, last evaluated,
     submitter-level conflicts.
   - Same amino-acid change via different nucleotide (PS1) and different change at same residue
     (PM5). Record accessions.
2. **gnomAD**: confirm AF/popmax/hom/hemi at variant level; gene constraint (pLI/LOEUF/mis-z).
3. **Ensembl**: canonical/MANE transcript, exact consequence, protein domain (PFAM/InterPro),
   regulatory overlap for non-coding, paralogues (for PM1 domain reasoning).
4. **UCSC**: phyloP/phastCons at the position, RepeatMasker/segdup context (cross-check QC),
   nearby benign variation.
5. **PubMed** (`lib/api_clients.py:pubmed_*`): structured queries —
   `"<gene> <variant/HGVS>"`, `"<gene> <key phenotype>"`, `"<gene> function"`,
   `"<gene> de novo|loss of function|missense"`. Classify each hit's relevance
   (same_variant / same_gene_same_phenotype / functional / background) and extract the finding.
   Record PMIDs; never cite a paper you did not retrieve.
6. **Gene–disease & phenotype match**: finalise ClinGen validity and compute the phenotype
   similarity (matched HPO, missing cardinal features, discordant features).

## Output (JSON) — array, one object per candidate
```json
{
  "candidate_id": "C0xx",
  "clinvar": {…}, "population_frequency": {…}, "gene_constraint": {…},
  "predictions": {…}, "gene_disease": {…}, "literature": [ {…} ],
  "regulatory_context": "…",
  "evidence_summary": "2–4 sentence neutral synthesis, citing refs",
  "disconfirming_notes": "anything found that argues AGAINST causality",
  "evidence_refs": [ … ]
}
```
