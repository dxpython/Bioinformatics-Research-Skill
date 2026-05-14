# Task router test cases

## 1. Purpose

Verify that routing in `SKILL.md` / `references/workflows/00_task_router.md` maps user requests to the correct workflows and supporting assets.

## 2. Test cases (20)

---

### Case 01: Bulk RNA-seq differential expression

**User input:**
> "I have a count matrix and group labels and want differentially expressed genes."

**Expected routing:**
- Primary: `references/workflows/bulk_rnaseq_workflow.md`
- Template: `references/templates/data_interpretation_template.md`

**Rationale:**
- Keywords: count matrix, groups, DEGs
- Data: bulk RNA-seq
- Goal: differential expression

---

### Case 02: scRNA-seq clustering

**User input:**
> "I have an h5ad file and want to explore cell clusters."

**Expected routing:**
- Primary: `references/workflows/single_cell_workflow.md`
- Template: `references/templates/figure_explanation_template.md` (UMAP, dot plots)
- Example: `references/examples/example_single_cell.md`

**Rationale:**
- Keywords: h5ad, clusters
- Data: scRNA-seq
- Goal: clustering / annotation

---

### Case 03: Literature-driven mechanism (multi-workflow)

**User input:**
> "I want to study the relationship between cuproptosis and autophagy in cervical cancer."

**Expected routing:**
- Primary: `references/workflows/literature_driven_mechanism_workflow.md`
- Secondary: `references/workflows/bulk_rnaseq_workflow.md` + `references/workflows/clinical_survival_workflow.md`
- Example: `references/examples/example_cancer_tcga.md`

**Rationale:**
- Keywords: cuproptosis, autophagy, cervical cancer
- Goal: curated gene sets → DE → survival → pathways; start with literature-derived gene lists

---

### Case 04: Microbiome

**User input:**
> "I have an OTU table and group labels."

**Expected routing:**
- Primary: `references/workflows/microbiome_workflow.md`
- Example: `references/examples/example_microbiome.md`

**Rationale:**
- Keywords: OTU table
- Data: 16S / microbiome
- Goal: community structure

---

### Case 05: Prognostic risk model

**User input:**
> "I want to build a prognostic risk model."

**Expected routing:**
- Primary: `references/workflows/clinical_survival_workflow.md` + `references/workflows/machine_learning_bioinformatics_workflow.md`
- Example: `references/examples/example_clinical_prediction.md`

**Rationale:**
- Keywords: prognosis, risk model
- Typical methods: LASSO-Cox, nomogram
- Ask for data source and outcome definition

---

### Case 06: Tumor immune microenvironment

**User input:**
> "I want to analyze the immune microenvironment in tumor samples."

**Expected routing:**
- Primary: `references/workflows/bulk_rnaseq_workflow.md` (immune deconvolution sections)
- Secondary: `references/workflows/clinical_survival_workflow.md`
- Example: `references/examples/example_cancer_tcga.md`

**Rationale:**
- Keywords: immune microenvironment, tumor
- Goal: CIBERSORT, ssGSEA, ESTIMATE, etc.

---

### Case 07: Pathway enrichment from a gene list

**User input:**
> "I already have a list of differential genes and want pathway enrichment."

**Expected routing:**
- Primary: `references/workflows/bulk_rnaseq_workflow.md` (enrichment sections)
- Template: `references/templates/data_interpretation_template.md`

**Rationale:**
- Keywords: gene list, pathway enrichment
- Input: DEG list
- Goal: GO/KEGG/GSEA

---

### Case 08: Spatial transcriptomics

**User input:**
> "I have 10x Visium data and want to analyze spatial expression patterns."

**Expected routing:**
- Primary: `references/workflows/spatial_transcriptomics_workflow.md`
- Example: `references/examples/example_spatial_transcriptomics.md`

**Rationale:**
- Keywords: Visium, spatial expression
- Data: spatial transcriptomics
- Goal: spatial domains, SVGs, deconvolution

---

### Case 09: Multi-omics integration

**User input:**
> "I have matched transcriptome and metabolome on the same samples and want joint analysis."

**Expected routing:**
- Primary: `references/workflows/multiomics_integration_workflow.md`
- Secondary: `references/workflows/bulk_rnaseq_workflow.md`
- Example: `references/examples/example_multiomics.md`

**Rationale:**
- Keywords: transcriptome, metabolome, joint analysis
- Data: paired multi-omics
- Goal: MOFA, DIABLO, pathway cross-talk

---

### Case 10: WGCNA

**User input:**
> "I want to use WGCNA to find gene modules related to a clinical trait."

**Expected routing:**
- Primary: `references/workflows/bulk_rnaseq_workflow.md` (WGCNA sections)
- Template: `references/templates/figure_explanation_template.md`

**Rationale:**
- Keywords: WGCNA, modules, clinical trait
- Goal: weighted gene co-expression networks

---

### Case 11: Pan-cancer gene profiling

**User input:**
> "I want expression and prognosis of one gene across 33 TCGA tumor types."

**Expected routing:**
- Primary: `references/workflows/bulk_rnaseq_workflow.md` + `references/workflows/clinical_survival_workflow.md`
- Example: `references/examples/example_cancer_tcga.md`

**Rationale:**
- Keywords: 33 cancers, pan-cancer
- Data: TCGA pan-cancer
- Goal: pan-tumor expression + survival

---

### Case 12: scRNA-seq + Visium integration

**User input:**
> "I have scRNA-seq and Visium from the same tumor sample and want integrated analysis."

**Expected routing:**
- Primary: `references/workflows/single_cell_workflow.md` + `references/workflows/spatial_transcriptomics_workflow.md`
- Strategy: cell2location / Tangram deconvolution
- Example: `references/examples/example_spatial_transcriptomics.md`

**Rationale:**
- Keywords: scRNA-seq, Visium, integration
- Data: single-cell + spatial
- Goal: deconvolution + cell mapping

---

### Case 13: Drug sensitivity association

**User input:**
> "I want to predict which genetic alterations are associated with chemotherapy sensitivity."

**Expected routing:**
- Primary: `references/workflows/machine_learning_bioinformatics_workflow.md`
- Resources: GDSC / CCLE when applicable

**Rationale:**
- Keywords: drug sensitivity, chemotherapy
- External pharmacogenomics data often required

---

### Case 14: Cell–cell communication

**User input:**
> "I want to analyze signaling between cell types in my scRNA-seq data."

**Expected routing:**
- Primary: `references/workflows/single_cell_workflow.md` (cell communication sections)
- Tools: CellChat / CellPhoneDB / NicheNet
- Example: `references/examples/example_single_cell.md`

**Rationale:**
- Keywords: single-cell, signaling, cell types
- Data: annotated scRNA-seq
- Goal: ligand–receptor inference

---

### Case 15: Microbiome + metabolome

**User input:**
> "I have 16S and metabolomics and want links between gut microbes and metabolites."

**Expected routing:**
- Primary: `references/workflows/microbiome_workflow.md` + `references/workflows/multiomics_integration_workflow.md`
- Example: `references/examples/example_microbiome.md`

**Rationale:**
- Keywords: 16S, metabolomics, microbiome, metabolites
- Goal: correlations, networks, mediation (optional)

---

### Case 16: Trajectory analysis

**User input:**
> "I want to trace T cell differentiation from naive to exhausted."

**Expected routing:**
- Primary: `references/workflows/single_cell_workflow.md` (pseudotime / trajectory)
- Tools: Monocle3 / scVelo / CytoTRACE

**Rationale:**
- Keywords: T cells, naive, exhausted, trajectory
- Data: scRNA-seq
- Goal: trajectory / RNA velocity

---

### Case 17: Nomogram construction

**User input:**
> "I want a nomogram to predict 1-, 3-, and 5-year survival."

**Expected routing:**
- Primary: `references/workflows/clinical_survival_workflow.md`
- Example: `references/examples/example_clinical_prediction.md`

**Rationale:**
- Keywords: nomogram, survival probability
- Goal: nomogram + calibration + DCA

---

### Case 18: Bulk RNA-seq from FASTQ

**User input:**
> "I have raw FASTQ and want a full RNA-seq pipeline from scratch."

**Expected routing:**
- Primary: `references/workflows/bulk_rnaseq_workflow.md` (upstream QC and quantification sections)
- Goal: QC → align → quantify → DE (link external best-practice pipelines as needed)

**Rationale:**
- Keywords: FASTQ, from scratch
- Data: raw sequencing

---

### Case 19: ML classification of subtypes

**User input:**
> "I want to use machine learning to classify tumor subtypes."

**Expected routing:**
- Primary: `references/workflows/machine_learning_bioinformatics_workflow.md`
- Methods: Random Forest / SVM / XGBoost + cross-validation

**Rationale:**
- Keywords: machine learning, classify, subtypes
- Goal: classification (not necessarily survival)
- Metrics: accuracy, AUC, F1

---

### Case 20: Reproducing paper figures

**User input:**
> "I want to reproduce figures from a paper: volcano plot, survival curve, and immune infiltration heatmap."

**Expected routing:**
- Template: `references/templates/figure_explanation_template.md`
- Workflows as needed: `references/workflows/bulk_rnaseq_workflow.md` + `references/workflows/clinical_survival_workflow.md`

**Rationale:**
- Keywords: reproduce, volcano, survival, immune heatmap
- Clarify raw data vs intermediate tables vs scores

---

## 3. Keyword → workflow summary

| Keywords | Workflow |
|----------|----------|
| count matrix / TPM / FPKM / DEG | `bulk_rnaseq_workflow.md` |
| h5ad / 10x / scRNA / clusters | `single_cell_workflow.md` |
| Visium / spatial / ST | `spatial_transcriptomics_workflow.md` |
| proteomics / protein abundance | `proteomics_workflow.md` |
| metabolomics / metabolite | `metabolomics_workflow.md` |
| OTU / ASV / 16S / microbiome | `microbiome_workflow.md` |
| methylation / ATAC / ChIP / epigenetic | `epigenomics_workflow.md` |
| mutation / MAF / VCF / variant | `genomics_variant_workflow.md` |
| survival / prognosis / KM / Cox / nomogram | `clinical_survival_workflow.md` |
| GSEA / GO / KEGG / enrichment | `bulk_rnaseq_workflow.md` (enrichment) |
| immune infiltration / CIBERSORT / TME | `bulk_rnaseq_workflow.md` (immune) |
| LASSO / Random Forest / prediction / ML | `machine_learning_bioinformatics_workflow.md` |
| multi-omics / integration / MOFA | `multiomics_integration_workflow.md` |
| mechanism / pathway hypothesis / literature | `literature_driven_mechanism_workflow.md` |

### Combined routes

| Scenario | Workflows |
|----------|-----------|
| Pan-cancer oncology | `bulk_rnaseq_workflow.md` + `clinical_survival_workflow.md` |
| Prognostic model | `clinical_survival_workflow.md` + `machine_learning_bioinformatics_workflow.md` |
| Literature-driven mechanism | `literature_driven_mechanism_workflow.md` + `bulk_rnaseq_workflow.md` + `clinical_survival_workflow.md` |
| Single-cell + spatial | `single_cell_workflow.md` + `spatial_transcriptomics_workflow.md` |
| Microbiome + metabolome | `microbiome_workflow.md` + `multiomics_integration_workflow.md` |

## 4. Pass / fail criteria

### Pass

1. Correct primary workflow(s).
2. Reasonable secondary workflows for multi-step questions.
3. Recommended templates match deliverables.
4. When inputs are ambiguous, the assistant asks for missing details.

### Fail

1. Wrong omics routing (e.g. scRNA routed only as bulk).
2. Missing essential modules (e.g. prognostic model without survival context).
3. Routing to files that do not exist in this repository.
4. No clarification on underspecified requests such as "analyze my data" with no format.
