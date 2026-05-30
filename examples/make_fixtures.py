"""
Generate small, valid toy inputs so the deterministic stages can be exercised without real
patient data. These are NOT biologically meaningful genomes — they are the minimum valid
VCF/BAM/HPO/notes needed to run validation, SV triage, IGV scripting and rendering end-to-end.

Outputs (under examples/fixtures/):
  snv_indel.vcf.gz (+.tbi)   VEP-annotated SNV/indel
  sv_cnv.vcf.gz   (+.tbi)    SV/CNV calls
  proband.bam     (+.bai)    a handful of reads around each locus
  gene.bed.gz     (+.tbi)    tiny gene BED (so SV overlap runs offline)
  phenotype.hpo              HPO terms
  clinical_notes.txt         clinical free text

Requires bgzip, tabix, samtools on PATH (present in the workflow sandbox).
Run:  python3 examples/make_fixtures.py
"""
from __future__ import annotations
import shutil
import subprocess
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parent / "fixtures"

# GRCh38 contig lengths for the chromosomes we touch (so coordinates are in-range).
CONTIGS = {"chr1": 248956422, "chr2": 242193529, "chr7": 159345973,
           "chr15": 101991189, "chr20": 64444167}

CSQ_FORMAT = ("Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|"
              "MANE_SELECT|gnomADg_AF|REVEL|CADD_PHRED|SpliceAI_pred_DS_AG|LoF")

# (chrom, pos, ref, alt, GT, DP, GQ, AD, csq_fields...)
SNVS = [
    ("chr20", 63452690, "G", "A", "0/1", 42, 99, "22,20",
     "A|missense_variant|MODERATE|KCNQ2|ENSG00000075043|NM_172107.4|c.821C>T|p.Ala274Val|"
     "NM_172107.4|0|0.94|29.4|0.02|"),
    ("chr2", 165310406, "AG", "A", "0/1", 35, 90, "20,15",
     "-|frameshift_variant|HIGH|SCN2A|ENSG00000136531|NM_021007.3|c.4880del|p.Gly1627fs|"
     "NM_021007.3|0|||0|LC"),
    ("chr1", 11796321, "C", "T", "1/1", 55, 99, "0,55",
     "T|missense_variant|MODERATE|MTHFR|ENSG00000177000|NM_005957.5|c.665C>T|p.Ala222Val|"
     "NM_005957.5|0.31|0.08|12.1|0.01|"),
    ("chr7", 117559593, "G", "A", "0/1", 40, 99, "21,19",
     "A|stop_gained|HIGH|CFTR|ENSG00000001626|NM_000492.4|c.1657C>T|p.Arg553Ter|"
     "NM_000492.4|0.0001||35.0|0|HC"),
]

# (chrom, pos, end, svtype, svlen, GT, PR_ref,PR_alt, SR_ref,SR_alt)  – one good, one artifact-y
SVS = [
    ("chr15", 42680000, 42690000, "DEL", -10000, "0/1", "30,2", "40,1"),   # weak support -> caution
    ("chr2", 165200000, 165400000, "DUP", 200000, "0/1", "20,12", "25,9"),  # spans SCN2A region
]

GENE_BED = [  # tiny offline gene BED for SV overlap
    ("chr2", 165200000, 165500000, "SCN2A"),
    ("chr15", 42600000, 42750000, "GENEX"),
    ("chr20", 63400000, 63500000, "KCNQ2"),
]

HPO = """HP:0001250 Seizure
HP:0011097 Epileptic spasms
HP:0001518 Small for gestational age
HP:0001263 Global developmental delay
HP:0002187 Intellectual disability, profound
HP:0010851 EEG with burst suppression
"""

NOTES = """Male infant, neonatal-onset epileptic encephalopathy. Tonic seizures from day 3 of
life evolving to epileptic spasms; EEG burst-suppression; profound global developmental delay.
No dysmorphism. Brain MRI normal at 2 months.
Family history: non-consanguineous parents, one healthy older sister, no epilepsy/DD in family.
Ancestry: Northern European.
Prior testing: chromosomal microarray normal (no pathogenic CNV >50 kb); epilepsy gene panel
(v3, 110 genes) reported no pathogenic variants (SNV/indel reporting only).
"""


def _need(tool: str) -> None:
    if not shutil.which(tool):
        sys.exit(f"ERROR: '{tool}' not found on PATH. Install htslib/samtools (see config/tools.yaml).")


def write_snv_vcf(path: Path) -> None:
    lines = ["##fileformat=VCFv4.2"]
    for c, ln in CONTIGS.items():
        lines.append(f"##contig=<ID={c},length={ln}>")
    lines += [
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS mapping quality">',
        '##INFO=<ID=FS,Number=1,Type=Float,Description="Phred strand bias">',
        f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from '
        f'Ensembl VEP. Format: {CSQ_FORMAT}">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND",
    ]
    for c, pos, ref, alt, gt, dp, gq, ad, csq in SNVS:
        info = f"MQ=60.0;FS=1.2;CSQ={csq}"
        lines.append(f"{c}\t{pos}\t.\t{ref}\t{alt}\t500\tPASS\t{info}\tGT:DP:GQ:AD\t"
                     f"{gt}:{dp}:{gq}:{ad}")
    path.write_text("\n".join(lines) + "\n")


def write_sv_vcf(path: Path) -> None:
    lines = ["##fileformat=VCFv4.2"]
    for c, ln in CONTIGS.items():
        lines.append(f"##contig=<ID={c},length={ln}>")
    lines += [
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">',
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PR,Number=.,Type=Integer,Description="Paired-read support ref,alt">',
        '##FORMAT=<ID=SR,Number=.,Type=Integer,Description="Split-read support ref,alt">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPROBAND",
    ]
    for i, (c, pos, end, svt, svlen, gt, pr, sr) in enumerate(SVS, 1):
        info = f"SVTYPE={svt};END={end};SVLEN={svlen}"
        lines.append(f"{c}\t{pos}\tSV{i}\tN\t<{svt}>\t60\tPASS\t{info}\tGT:PR:SR\t{gt}:{pr}:{sr}")
    path.write_text("\n".join(lines) + "\n")


def write_gene_bed(path: Path) -> None:
    path.write_text("".join(f"{c}\t{s}\t{e}\t{n}\n" for c, s, e, n in GENE_BED))


def write_bam(path: Path) -> None:
    sam = ["@HD\tVN:1.6\tSO:coordinate"]
    for c, ln in CONTIGS.items():
        sam.append(f"@SQ\tSN:{c}\tLN:{ln}")
    rid = 0
    for c, pos, *_ in SNVS:
        for k in range(10):
            flag = 16 if k % 2 else 0           # mix strands
            start = pos - 50 + (k % 3)
            seq = "A" * 100
            qual = "I" * 100
            sam.append(f"r{rid}\t{flag}\t{c}\t{start}\t60\t100M\t*\t0\t0\t{seq}\t{qual}")
            rid += 1
    for c, pos, end, *_ in SVS:                  # a few reads at each breakpoint
        for p in (pos, end):
            for k in range(6):
                sam.append(f"r{rid}\t{0 if k%2 else 16}\t{c}\t{p-30}\t40\t100M\t*\t0\t0\t"
                           f"{'A'*100}\t{'I'*100}")
                rid += 1
    sam_text = "\n".join(sam) + "\n"
    tmp_sam = path.with_suffix(".sam")
    tmp_sam.write_text(sam_text)
    subprocess.run(f"samtools view -bS {tmp_sam} | samtools sort -o {path} -",
                   shell=True, check=True)
    subprocess.run(["samtools", "index", str(path)], check=True)
    tmp_sam.unlink()


def bgzip_tabix(path: Path, preset: str) -> None:
    subprocess.run(["bgzip", "-f", str(path)], check=True)
    subprocess.run(["tabix", "-f", "-p", preset, str(path) + ".gz"], check=True)


def main() -> None:
    for t in ("bgzip", "tabix", "samtools"):
        _need(t)
    OUT.mkdir(parents=True, exist_ok=True)
    write_snv_vcf(OUT / "snv_indel.vcf"); bgzip_tabix(OUT / "snv_indel.vcf", "vcf")
    write_sv_vcf(OUT / "sv_cnv.vcf"); bgzip_tabix(OUT / "sv_cnv.vcf", "vcf")
    write_gene_bed(OUT / "gene.bed"); bgzip_tabix(OUT / "gene.bed", "bed")
    write_bam(OUT / "proband.bam")
    (OUT / "phenotype.hpo").write_text(HPO)
    (OUT / "clinical_notes.txt").write_text(NOTES)
    print(f"Fixtures written to {OUT}/")
    for f in sorted(OUT.iterdir()):
        print("  ", f.name)


if __name__ == "__main__":
    main()
