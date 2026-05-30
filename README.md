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

## 7. Setup

```bash
pip install -r requirements.txt          # requests, PyYAML, jsonschema, anthropic
# system tools (bcftools/samtools/bedtools/igv/bgzip/tabix) are provisioned per config/tools.yaml
```

The **gene BED** used by the SV stage is fetched automatically from an authoritative source —
no manual curation:

```bash
# GENCODE basic annotation (default, versioned, GRCh38/GRCh37), MANE-Select genes only:
python3 -m lib.fetch_gene_bed --assembly GRCh38 --mane-only --out $WORKFLOW_DATA/tracks/gene.bed.gz

# or UCSC MANE/ncbiRefSeqSelect track:
python3 -m lib.fetch_gene_bed --assembly GRCh38 --source ucsc --out tracks/gene.bed.gz

# optionally restrict to a PanelApp Australia panel (coordinates from GENCODE, genes from the panel):
python3 -m lib.fetch_gene_bed --assembly GRCh38 --panelapp 250 --out panel.bed.gz
```

Sources: **GENCODE** (EBI, authoritative reference gene model) and **UCSC** are reached directly;
**PanelApp Australia** (`panelapp-aus.org`) supplies the clinical gene whitelist + green/
amber/red confidence. If a source is unreachable the fetcher exits with explicit guidance and a
`--gene-list` fallback (supply a plain symbol list).

## 8. Validating inputs

Inputs are validated before any analysis runs; the workflow aborts with actionable, *sourcing-
aware* errors (how to annotate/compress/index/obtain each file) rather than failing deep in a
stage. Run the check standalone:

```bash
python3 -m lib.validate_inputs \
  --snv-vcf snv_indel.vcf.gz --sv-vcf sv_cnv.vcf.gz --bam proband.bam \
  --hpo phenotype.hpo --notes clinical_notes.txt --gene-bed tracks/gene.bed.gz
```

It checks: bgzip+tabix indexing, VEP `CSQ` presence, BAM sort/index, HPO id validity, non-empty
notes, and chr-naming concordance across inputs vs the declared assembly.

## 9. Running

```bash
python3 orchestrator.py \
  --case-id RD-2026-0142 \
  --snv-vcf  inputs/snv_indel.vcf.gz \
  --sv-vcf   inputs/sv_cnv.vcf.gz \
  --bam      inputs/proband.bam \
  --hpo      inputs/phenotype.hpo \
  --notes    inputs/clinical_notes.txt \
  --gene-bed tracks/gene.bed.gz \
  --assembly GRCh38 \
  --out      results/

# results/report.json   results/report.html   results/igv/*.png   results/audit/
```

The reasoning stages require an LLM backend wired into `AgentRunner.complete()` (see the
docstring in `orchestrator.py`); deterministic stages, validation, fixtures, gene-BED fetching
and rendering run without it.

## 10. Try it without real data (fixtures + smoke test)

```bash
python3 examples/make_fixtures.py     # toy VCF/BAM/HPO/notes in examples/fixtures/
python3 tests/smoke_test.py           # exercises validation, SV triage, ACMG, schema, render
```

The fixtures are valid-but-synthetic inputs to exercise the deterministic stages offline. A
worked end-to-end output (no LLM needed to view it) is in
[`examples/example_report.json`](examples/example_report.json) →
[`examples/example_report.html`](examples/example_report.html).

See [`config/thresholds.yaml`](config/thresholds.yaml) for every tunable QC / frequency cut.
