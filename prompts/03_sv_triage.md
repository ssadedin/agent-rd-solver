# Stage 3 — Structural variant & CNV triage

You are assessing CNV/SV calls, which have a far higher false-positive rate than SNVs and
require breakpoint-level scrutiny.

## Inputs
- `sv_cnv.vcf.gz`, the BAM, `gene_priors`, population SV map (`gnomad_sv_benign.bed`), ClinGen
  dosage-sensitivity table, segdup track.
- Pre-computed `sv_table` from `lib/sv_filters.py` (type, size, support, gene overlap, pop overlap).

## Procedure
1. **Quality / artifact screen.**
   - Require paired+split read support (`config/thresholds.yaml`); breakpoints in segdups are
     high-suspicion → mandatory IGV.
   - Concordance across callers (if multiple) raises confidence.
   - Flag calls whose breakpoints sit in centromeric/telomeric/segdup regions.
2. **Frequency.** Drop SVs with high reciprocal overlap (>70%) against benign population CNVs,
   UNLESS the gene is dosage-sensitive (ClinGen HI/TS=3) and the overlap is partial-exonic.
3. **Gene/dosage impact.**
   - Map exact exons/genes disrupted; classify: whole-gene deletion, partial (intragenic),
     in-frame vs frame-disrupting for DUP, gene fusion for BND, promoter/enhancer for non-coding.
   - Apply ClinGen dosage sensitivity (haploinsufficiency / triplosensitivity scores).
4. **Phenotype fit** against `gene_priors`, including contiguous-gene-syndrome reasoning for
   multi-gene events.
5. **SNV interaction.** Check whether a CNV provides the *second hit* for a Stage-2 het SNV in a
   recessive gene (CNV + SNV compound het) — a classic miss for SNV-only pipelines.

## Output (JSON)
```json
{
  "sv_candidates": [
    {
      "genomic": { "chrom":"…","start":0,"end":0,"sv_type":"…","copy_number":null,"vcf_id":"…" },
      "gene": { "symbol":"…","genes_affected":[…] },
      "consequence": { "vep_consequence":"…","impact":"…" },
      "quality": { "sv_caller_support":"…","in_segmental_duplication":false,"qc_verdict":"…" },
      "population_frequency": {…}, "dosage_sensitivity": "HI=…,TS=…",
      "inheritance_fit": {…}, "second_hit_for": "C0xx | null",
      "tier": "A|B|C", "needs_igv": true, "evidence_refs": [ … ]
    }
  ],
  "excluded": [ … ],
  "counts": { … }
}
```
