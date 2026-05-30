# Undiagnosed Rare Disease Diagnostic Workflow

An agentic, adversarially-reviewed workflow for identifying candidate causal variants
in an undiagnosed rare disease patient from whole-genome sequencing data.

It is designed to mirror, end-to-end, the reasoning a multidisciplinary research team
(clinical geneticist + bioinformatician + variant scientist) would perform when working
up a difficult case — but with an explicit *red-team* stage that actively tries to
**eliminate** every candidate before it is reported.

---

## 1. Inputs

| Input | Description |
|-------|-------------|
| `snv_indel.vcf.gz` | VEP-annotated SNV/indel calls from WGS (`CSQ` INFO field expected) |
| `sv_cnv.vcf.gz` | CNV + structural-variant calls (e.g. Manta/GRIDSS/CNVnator/cn.MOPS) |
| `proband.bam` | Aligned reads (+ `.bai`), used for visual QC of every candidate |
| `phenotype.hpo` | One HPO term per line (`HP:0001250`), or HPO id + label |
| `clinical_notes.txt` | Free text: age of onset, family history, prior testing, consanguinity, ethnicity, exam findings |

> The workflow is **single-proband** by default. If parental/sibling VCFs or BAMs are
> supplied (`--parents`, `--siblings`), inheritance segregation is upgraded from
> *inferred* (from clinical text) to *observed*.

## 2. Outputs

| Output | Description |
|--------|-------------|
| `report.json` | Structured report conforming to [`schema/candidate_report.schema.json`](schema/candidate_report.schema.json) |
| `report.html` | Self-contained human-friendly rendering of the JSON (no server needed) |
| `igv/*.png` | IGV screenshots referenced by candidates |
| `audit/` | Per-stage prompts, tool commands, raw API responses (full provenance) |

The report contains:
1. **Candidate variants** with structured, ACMG-aligned (but not ACMG-limited) evidence.
2. **Pending evidence** — specific questions/experiments that would materially change a
   classification, each with a justification and a *critically-reviewed* experimental plan.
3. **Excluded-variant summary** — what was filtered and *why* (auditability).

The workflow **never stops early to ask for information**. Missing information becomes a
`pending_evidence` item; the run always completes a full report on the data at hand.

---

## 3. Architecture

```
                    ┌─────────────────────────────────────────────────────────┐
                    │                     orchestrator.py                       │
                    │  (drives stages, owns state, enforces adversarial gate)   │
                    └─────────────────────────────────────────────────────────┘
        deterministic libs  ◄──────────────┴──────────────►  reasoning agents (LLM + tools)
        (lib/*.py)                                            (prompts/*.md)

  Stage 0  Intake & phenotype structuring ........... agent  → prompts/00_intake.md
  Stage 1  Phenotype-driven gene prioritisation ..... agent  → prompts/01_gene_prioritization.md
  Stage 2  SNV/indel triage & QC .................... lib+agent → lib/vcf_filters.py, prompts/02_snv_triage.md
  Stage 3  SV/CNV triage & QC ....................... lib+agent → lib/sv_filters.py, prompts/03_sv_triage.md
  Stage 4  Evidence gathering (ClinVar/PubMed/...) ... lib+agent → lib/api_clients.py, prompts/04_evidence.md
  Stage 5  IGV visual review ........................ tool+agent → lib/igv_runner.py, prompts/05_igv_review.md
  Stage 6  ACMG-aligned classification .............. lib+agent → lib/acmg.py, prompts/06_acmg.md
  Stage 7  ADVERSARIAL review (red team) ............ agent  → prompts/07_adversarial_review.md
  Stage 8  Synthesis & ranking ...................... agent  → prompts/08_synthesis.md
  Stage 9  Pending-evidence proposal + plan review .. agent  → prompts/09_pending_evidence.md
  Render   JSON → HTML .............................. lib    → render/render_html.py
```

### Division of labour: deterministic vs. agentic
- **Deterministic code** (`lib/`) does everything that must be reproducible and verifiable:
  VCF parsing, hard QC thresholds, frequency filtering, ACMG point arithmetic, API calls,
  IGV scripting. These produce *facts*.
- **Reasoning agents** (`prompts/`) do everything requiring judgement: phenotype matching,
  weighing conflicting evidence, the adversarial critique, and prose synthesis. Agents are
  **only allowed to assert facts that trace to a deterministic tool output or a cited API
  response** (see §5 anti-hallucination contract).

---

## 4. The adversarial design (why candidates must *survive*, not merely *score*)

A naive pipeline ranks variants by a score and reports the top N. That over-reports.
Real diagnostic rigour is **eliminative**: a variant is reported only if a determined
sceptic *fails* to disqualify it.

Two mechanisms enforce this:

1. **Per-candidate red team (Stage 7).** A separate agent persona, with no stake in the
   prior reasoning, receives each candidate and is instructed to *destroy* it. It must
   enumerate every disqualifying hypothesis across nine attack surfaces (see
   `prompts/07_adversarial_review.md`): technical artifact, mapping ambiguity, frequency,
   inheritance inconsistency, phenotype mismatch, gene–disease validity, ClinVar conflict,
   in-silico weakness, and "better explanation exists elsewhere". Each challenge is graded
   `fatal | major | minor` and must be explicitly `refuted | mitigated | unrefuted`.
2. **Survival gate.** The orchestrator demotes any candidate with an `unrefuted` `fatal`
   challenge to the excluded list; the proposer agent is given one rebuttal round per
   challenge (with tool access) before the gate closes. This is a structured debate, not a
   single pass.

The same critique discipline is applied to **proposed experiments** (Stage 9): every
suggested wet-lab method is itself red-teamed for feasibility, confounders, and whether it
can actually move the classification.

---

## 5. Anti-hallucination contract (binding on every agent)

Every agent prompt ends with this contract:
- Cite the source of every factual claim: a tool stdout line, a JSON field, or an API
  accession (`ClinVar VCV…`, `PMID:…`, `gnomAD v4 …`). Uncited claims are invalid.
- If a needed fact is unavailable, output `UNKNOWN` and (if it matters) raise a
  `pending_evidence` item — never guess.
- Never invent gene–disease relationships, frequencies, PMIDs, or ACMG codes.
- Distinguish *observed* (from data) from *inferred* (from clinical text) from *assumed*.

---

## 6. Assumed runtime environment

Provisioned and on `PATH` (invoked via Bash in the sandbox):
`bcftools`, `tabix/bgzip`, `samtools`, `bedtools`, `igv` (+ `xvfb-run`), `python3`.

Reachable APIs (keys read from env where applicable):
- **ClinVar / PubMed** via NCBI E-utilities (`NCBI_API_KEY` recommended for rate limits)
- **Ensembl REST** (`https://rest.ensembl.org`) — GRCh38 assumed; GRCh37 via `grch37.rest…`
- **UCSC** REST/DAS (`https://api.genome.ucsc.edu`) — repeats, segdups, conservation tracks
- **gnomAD** GraphQL (`https://gnomad.broadinstitute.org/api`) for population frequency

Reference data expected under `$WORKFLOW_DATA` (configurable in `config/tools.yaml`):
HPO `phenotype_to_genes.txt`, OMIM/MONDO maps, gnomAD constraint, MANE, ClinGen gene–disease
validity, RefSeq/Ensembl GTF, segdup/repeat BED tracks.

## 7. Running

```bash
python3 orchestrator.py \
  --case-id RD-2026-0142 \
  --snv-vcf  inputs/snv_indel.vcf.gz \
  --sv-vcf   inputs/sv_cnv.vcf.gz \
  --bam      inputs/proband.bam \
  --hpo      inputs/phenotype.hpo \
  --notes    inputs/clinical_notes.txt \
  --assembly GRCh38 \
  --out      results/

# results/report.json   results/report.html   results/igv/*.png   results/audit/
```

See [`config/thresholds.yaml`](config/thresholds.yaml) for every tunable QC / frequency cut.
