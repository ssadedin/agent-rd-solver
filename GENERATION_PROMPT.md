# Prompt

Create a workflow designed to solve  an undiagnosed rare disease patient. The inputs will be:

* A VCF file containing VEP-annotated SNV/indel calls derived from whole genome sequencing and a separate VCF file containing CNV and SV calls, a BAM file containing aligned reads
* a list of HPO terms describing the patient phenotypes
* a text description describing clinical factors including age of onset, family history, results from prior testing, and other unstructred information about the patient

The workflow should take all reasonable steps that a set of human researchers would take, including identifying candidate variants, reviewing them for quality and other factors that would filter them out from consideration.  The code should assume a set of standard tools will be provisioned for the agent to execute using a Bash command within a sandboxed environment. It should also assume that  API access to ClinVar, PubMed, Ensembl and UCSC is available, and the IGV tool is available to run and generate screenshots of regions of interest from the alignments.

The workflow should ensure all results are rigorously reviewed and therefore employ adversarial review to ensure all possible flaws in hypothesised variants are closely examined for factors that would eliminate them  from consideration.

Outputs should be (a) proposed variants for consideration along with structured evidence supporting the reasoning. This evidence should align with, but not be restricted to that used in the ACMG guidelines for clinical variant classification (b) a list of pending evidence that could add significant clarity to one or more of the variants under consideration if a specific question was clarified (outside of the agentic's system abilitiy). This could be a request for clinical data, or an experimental method to be applied.  If an experimental method is suggested, it must be accompanied by a precise justification for exactly how it will improve the understanding and a detailed experimental plan, itself subject to critical review. Do NOT stop early to request information - fully work to completion with all the information you have to produce a comprehensive output report for all candidate variants.

The output should be written in (a) a JSON structured format that can be determined (develop a reasonable schema for it) accompanied by (b) an HTML file that displays the JSON in a human friendly manner.
