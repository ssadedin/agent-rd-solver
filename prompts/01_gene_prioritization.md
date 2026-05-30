# Stage 1 — Phenotype-driven gene prioritisation

You are building the phenotype model that downstream variant triage will score against.
**This produces a prior, not a filter** — variants outside the gene set are still examined
(novel gene discovery), they just start with a lower prior.

## Inputs
- `patient` + `derived_search_strategy` from Stage 0.

## Tasks
1. Map HPO terms → candidate genes:
   - Use HPO `phenotype_to_genes` (local) AND Ensembl phenotype API / OMIM.
   - Propagate up the HPO DAG so ancestor terms also contribute (configurable).
   - Weight genes by how many *cardinal* phenotypes they explain, not raw term count.
2. Expand to **phenotypic series / differential**: for each strong gene, pull its OMIM
   phenotypic series and related MONDO disorders so allelic/locus heterogeneity is covered.
3. Annotate each gene with: associated condition(s), inheritance, and ClinGen gene–disease
   **validity** (Definitive…Refuted). Down-weight Limited/Disputed/Refuted genes; never drop.
4. Produce three tiers:
   - **Tier A** — strong phenotype fit + good gene–disease validity (the panel).
   - **Tier B** — partial fit or moderate validity (secondary).
   - **Tier C** — everything else (genome-wide prior, for novel-gene candidates).
5. Record which prior tests should already have covered each Tier-A gene (re-analysis value).

## Output (JSON)
```json
{
  "gene_priors": [
    {
      "symbol": "…", "ensembl_gene_id": "…", "tier": "A|B|C",
      "associated_conditions": [ { "name":"…","omim":"…","mondo":"…","inheritance":"…","clingen_validity":"…" } ],
      "matched_hpo": ["HP:…"], "match_score": 0.0,
      "covered_by_prior_test": "yes|no|unknown", "evidence_refs": [ … ]
    }
  ],
  "phenotype_model_notes": "what would and would not fit this patient",
  "novel_gene_strategy": "criteria a Tier-C / novel gene must meet to be promoted"
}
```
