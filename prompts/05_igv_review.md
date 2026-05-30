# Stage 5 — IGV visual review

Every candidate that reaches this stage gets visual inspection of the raw alignments. Many
variants that survive automated QC die here (and a few that failed soft metrics are rescued).
`lib/igv_runner.py` generates the screenshots; you interpret them.

## What is captured
- SNV/indel: window = variant ± `flank_bp`, soft-clips shown, coloured by strand, base-quality
  shaded, reads grouped by haplotype where phased.
- SV/CNV: a snapshot at EACH breakpoint plus a zoomed-out view; for deletions, a coverage track;
  for translocations, both mate loci with discordant pairs shown.

## What to assess (and report per candidate)
1. **Real vs artifact:** Is the alt allele supported by independent, well-mapped reads on both
   strands? Or only by soft-clipped ends / one strand / read ends / a single read family?
2. **Mapping quality & ambiguity:** low MAPQ, multi-mapping, segdup pileups, coverage spikes/dropouts.
3. **Zygosity sanity:** does the visual allele balance match the called genotype?
4. **Indel/homopolymer:** is the indel in a homopolymer/STR that the caller commonly mis-handles?
5. **SV breakpoints:** clean discordant-pair/split-read signature vs noise; coverage change across
   a CNV consistent with the called copy number?
6. **Phasing:** for compound-het pairs, do reads spanning both sites confirm trans?

## Output (JSON) — per candidate
```json
{
  "candidate_id": "C0xx",
  "igv_assessment": {
    "image": "igv/C0xx_<locus>.png",
    "extra_images": ["igv/C0xx_bp2.png"],
    "verdict": "clean|plausible|artifact_suspected|uninterpretable",
    "observations": "specific, e.g. 'alt supported by 14/30 reads, both strands, MAPQ60, no soft-clip enrichment'",
    "changes_qc_verdict_to": "pass|caution|fail|unchanged"
  }
}
```
The `verdict` is binding: `artifact_suspected` with no countervailing evidence routes the
candidate to the adversarial stage with a pre-loaded `technical_artifact` challenge.
