# Task Router: Intelligent Task Recognition and Workflow Selection

## Overview

The Task Router is the brain of the Universal Bioinformatics Research Skill. When a user submits a query, the router analyzes the input to determine the task category, selects the appropriate workflow(s), and defines the response structure.

## Task Categories

### A. Research Design

**Trigger Keywords**: "research direction", "study design", "topic selection", "what should I study", "how to design", "research plan", "feasibility"

**Input Indicators**:
- User describes a broad research interest without specific data
- User asks about which databases or datasets to use
- User asks about experimental design

**Recommended Workflow**: `research_design_template.md` + relevant domain workflow

**Output Structure**:
1. Background assessment
2. Core scientific question formulation
3. Hypothesis
4. Data source recommendations
5. Analysis pipeline recommendation
6. Expected figures and tables
7. Innovation points
8. Risk assessment

**Risks**: Avoid overpromising novelty; verify database data availability before recommending.

**Anti-overinterpretation**: Do not claim a direction is "guaranteed to publish" or "definitely novel."

---

### B. Database Search and Data Download

**Trigger Keywords**: "download data", "TCGA", "GEO", "NCBI", "ENCODE", "ArrayExpress", "database", "public data", "GSE", "PRJNA"

**Input Indicators**:
- User mentions specific database names
- User asks how to obtain certain data types
- User provides accession numbers (GSE, PRJNA, TCGA barcode)

**Recommended Workflow**: Provide database-specific download instructions

**Output Structure**:
1. Database identification
2. Data availability check guidance
3. Download method (R/Python/web)
4. Data format description
5. Preprocessing requirements

**Risks**: Data may be restricted access; verify licensing.

---

### C. Bulk RNA-seq Expression Matrix Analysis

**Trigger Keywords**: "count matrix", "expression matrix", "TPM", "FPKM", "differential expression", "DEGs", "DESeq2", "limma", "edgeR", "bulk RNA-seq", "transcriptome"

**Input Indicators**:
- User provides a gene × sample matrix
- User mentions count/TPM/FPKM data
- User has grouping information (e.g., treatment vs. control)

**Recommended Workflow**: `bulk_rnaseq_workflow.md`

**Output Structure**:
1. Data format validation
2. Expression type identification (count/TPM/FPKM)
3. Recommended analysis pipeline
4. Code generation
5. Result interpretation
6. Figure list

**Risks**: Using wrong normalization (e.g., DESeq2 requires raw counts, not TPM).

**Anti-overinterpretation**: DEGs are statistical associations, not functional validation.

---

### D. Differential Result Interpretation

**Trigger Keywords**: "interpret results", "what do these DEGs mean", "result table", "significant genes", "upregulated/downregulated"

**Input Indicators**:
- User provides a differential analysis result table (gene, log2FC, p-value)
- User asks about the meaning of specific genes

**Recommended Workflow**: `data_interpretation_template.md` + `bulk_rnaseq_workflow.md` (enrichment section)

**Output Structure**:
1. Table structure assessment
2. Threshold-based filtering
3. Top gene annotation
4. Enrichment analysis recommendation
5. Biological interpretation

**Risks**: Do not assume significance equals biological importance.

---

### E. Single-Cell Analysis

**Trigger Keywords**: "single cell", "scRNA-seq", "10X", "h5ad", "Seurat", "scanpy", "cell types", "clustering", "UMAP", "cell annotation", "marker genes", "multiome", "ATAC-seq", "snATAC-seq", "CITE-seq", "ADT", "antibody-derived", "surface protein", "RNA velocity", "chromatin accessibility"

**Input Indicators**:
- User provides h5ad, 10X, or Seurat object files
- User mentions cell types, clusters, or marker genes
- User asks about cell communication or trajectory

**Recommended Workflow**: `single_cell_workflow.md`

**Output Structure**:
1. Data format identification
2. QC metrics recommendation
3. Analysis pipeline (normalization → HVG → PCA → clustering → annotation)
4. Advanced analyses (pseudotime, communication, CNV)
5. Visualization recommendations

**Risks**: Over-clustering; incorrect cell annotation; batch effects.

**Anti-overinterpretation**: Cluster annotations are hypothesis-driven, not definitive cell identities.

---

### F. Spatial Transcriptomics Analysis

**Trigger Keywords**: "spatial", "Visium", "MERFISH", "Slide-seq", "spatial transcriptomics", "spatial spots", "tissue section", "spatial clustering"

**Input Indicators**:
- User provides spatial data (h5ad with spatial coordinates)
- User mentions tissue sections or spatial patterns
- User asks about deconvolution or spatial co-localization

**Recommended Workflow**: `spatial_transcriptomics_workflow.md`

**Note**: If user says "h5ad file," check whether spatial coordinates exist. If yes → spatial workflow. If no → single-cell workflow.

**Output Structure**:
1. Spatial data validation
2. QC for spatial spots
3. Spatial clustering
4. Deconvolution recommendations
5. Spatial communication analysis
6. Integration with scRNA-seq

**Risks**: Spot-level data is a mixture of cells; deconvolution is probabilistic.

---

### G. Proteomics Analysis

**Trigger Keywords**: "proteomics", "protein expression", "mass spectrometry", "MaxQuant", "TMT", "iTRAQ", "PPI network", "differential proteins"

**Input Indicators**:
- User provides protein abundance matrix
- User mentions mass spectrometry data
- User asks about PPI networks

**Recommended Workflow**: `proteomics_workflow.md`

**Output Structure**:
1. Data format assessment
2. Missing value strategy
3. Normalization
4. Differential analysis
5. PPI network construction
6. Functional enrichment

**Risks**: High missing value rates; batch effects; post-translational modifications not captured.

---

### H. Metabolomics Analysis

**Trigger Keywords**: "metabolomics", "metabolites", "LC-MS", "GC-MS", "PLS-DA", "OPLS-DA", "metabolic pathway", "metabolite abundance"

**Input Indicators**:
- User provides metabolite abundance matrix
- User mentions LC-MS or GC-MS data
- User asks about metabolic pathway analysis

**Recommended Workflow**: `metabolomics_workflow.md`

**Output Structure**:
1. Data quality assessment
2. Normalization strategy
3. Multivariate analysis (PCA, PLS-DA, OPLS-DA)
4. Differential metabolites
5. Pathway enrichment
6. Metabolite-gene correlation

**Risks**: Metabolite identification confidence levels vary; isomers may confound results.

---

### I. Microbiome Analysis

**Trigger Keywords**: "microbiome", "16S", "OTU", "ASV", "metagenomics", "alpha diversity", "beta diversity", "LEfSe", "gut microbiota", "microbial community"

**Input Indicators**:
- User provides OTU/ASV table with taxonomy
- User mentions 16S rRNA or metagenomic data
- User asks about diversity or differential abundance

**Recommended Workflow**: `microbiome_workflow.md`

**Output Structure**:
1. Data format validation
2. Alpha diversity analysis
3. Beta diversity analysis
4. Differential abundance testing
5. Functional prediction (PICRUSt2)
6. Microbiome-metabolome integration

**Risks**: Compositional data bias; low taxonomic resolution at species level for 16S.

---

### J. Epigenomics Analysis

**Trigger Keywords**: "methylation", "ATAC-seq", "ChIP-seq", "histone", "epigenetic", "DMP", "DMR", "peak calling", "motif", "chromatin accessibility"

**Input Indicators**:
- User provides methylation data (beta values, IDAT files)
- User mentions ATAC-seq or ChIP-seq peak files
- User asks about regulatory elements

**Recommended Workflow**: `epigenomics_workflow.md`

**Output Structure**:
1. Data type identification
2. Quality control
3. Differential analysis (DMP/DMR or differential peaks)
4. Peak annotation
5. Motif enrichment
6. Integration with expression data

**Risks**: Methylation changes may not correlate with expression; peak calling parameters greatly affect results.

---

### K. Genomic Variant Analysis

**Trigger Keywords**: "mutation", "SNP", "InDel", "CNV", "variant", "MAF", "VCF", "TMB", "MSI", "mutation signature", "GWAS", "oncoplot"

**Input Indicators**:
- User provides VCF or MAF files
- User mentions somatic or germline variants
- User asks about mutation patterns

**Recommended Workflow**: `genomics_variant_workflow.md`

**Output Structure**:
1. Variant data format validation
2. Variant statistics
3. Visualization (oncoplot, waterfall)
4. Mutation signature analysis
5. Clinical correlation
6. Pathway-level variant impact

**Risks**: Variant calling quality; germline contamination in somatic analysis; population stratification in GWAS.

---

### L. Clinical Prognosis Analysis

**Trigger Keywords**: "survival", "prognosis", "Kaplan-Meier", "Cox regression", "nomogram", "risk score", "OS", "DFS", "PFS", "hazard ratio", "time-dependent ROC"

**Input Indicators**:
- User provides clinical data with survival time and event status
- User asks about prognostic factors
- User mentions risk stratification

**Recommended Workflow**: `clinical_survival_workflow.md`

**Output Structure**:
1. Clinical data validation
2. Univariate analysis
3. Multivariate Cox regression
4. Risk score construction
5. Nomogram
6. Validation (internal + external)

**Risks**: Overfitting; insufficient events; selection bias in retrospective cohorts.

**Anti-overinterpretation**: A prognostic association is not a causal relationship.

---

### M. Machine Learning Modeling

**Trigger Keywords**: "machine learning", "LASSO", "SVM", "random forest", "XGBoost", "biomarker", "feature selection", "prediction model", "diagnostic model", "AUC", "ROC"

**Input Indicators**:
- User wants to build a predictive model
- User asks about feature selection methods
- User mentions training/validation datasets

**Recommended Workflow**: `machine_learning_bioinformatics_workflow.md`

**Output Structure**:
1. Data preparation
2. Feature selection (multiple methods)
3. Model training and comparison
4. Cross-validation
5. External validation
6. Model interpretation (SHAP, feature importance)

**Risks**: Overfitting; data leakage; imbalanced classes; lack of external validation.

---

### N. Multi-omics Integration

**Trigger Keywords**: "multi-omics", "integration", "transcriptome + proteome", "transcriptome + metabolome", "joint analysis", "cross-omics", "network integration"

**Input Indicators**:
- User has data from multiple omics platforms
- User asks about correlation between different omics layers
- User wants to build integrative networks

**Recommended Workflow**: `multiomics_integration_workflow.md`

**Output Structure**:
1. Data harmonization
2. Cross-omics correlation
3. Network construction
4. Module detection
5. Pathway consistency analysis
6. Integrative visualization

**Risks**: Different sample sizes across omics; different measurement scales; forced integration without biological justification.

---

### O. Figure Interpretation

**Trigger Keywords**: "explain this figure", "what does this plot show", "interpret the heatmap", "read this UMAP", "explain the volcano plot"

**Input Indicators**:
- User uploads a screenshot or image
- User asks about a specific figure type
- User wants manuscript text for a figure

**Recommended Workflow**: `figure_explanation_template.md`

**Output Structure**:
1. Figure type identification
2. Key observations
3. Statistical evidence
4. Biological interpretation
5. Manuscript-ready Results text
6. Improvement suggestions

**Risks**: Cannot extract exact values from screenshots; figure quality may limit interpretation.

---

### P. Manuscript Writing

**Trigger Keywords**: "write Methods", "write Results", "figure legend", "Discussion", "manuscript", "paper writing", "journal submission"

**Input Indicators**:
- User asks for specific manuscript sections
- User provides analysis results and wants them written up
- User asks for academic English editing

**Recommended Workflow**: `paper_writing_template.md`

**Output Structure**:
1. Section identification
2. Appropriate academic tone
3. Software citations
4. Statistical reporting standards
5. Figure/table references

**Risks**: Do not fabricate results or references.

---

### Q. Code Generation / Debugging

**Trigger Keywords**: "give me code", "error message", "debug", "code not working", "how to run", "write a script", "R code", "Python code"

**Input Indicators**:
- User requests analysis code
- User pastes an error message
- User asks how to implement a specific analysis

**Recommended Workflow**: `code_generation_template.md` + relevant analysis workflow

**Output Structure**:
1. Environment specification (R+Rlib or Python+uv)
2. Input file format
3. Complete runnable code
4. Output description
5. Common errors and solutions

**Risks**: Environment differences between systems; version-specific API changes.

---

### R. Reviewer Response

**Trigger Keywords**: "reviewer", "revision", "response to reviewer", "reviewer comments", "rebuttal", "peer review", "major revision", "minor revision"

**Input Indicators**:
- User pastes reviewer comments
- User asks how to address specific criticisms
- User needs additional analyses for revision

**Recommended Workflow**: `reviewer_response_template.md`

**Output Structure**:
1. Comment categorization
2. Response strategy
3. Additional analyses if needed
4. Diplomatic phrasing
5. Manuscript modification summary

**Risks**: Avoid being confrontational; acknowledge valid criticisms.

---

### S. Experimental Validation Design

**Trigger Keywords**: "experimental validation", "qRT-PCR", "Western blot", "IHC", "knockdown", "overexpression", "in vivo", "in vitro", "wet lab"

**Input Indicators**:
- User has bioinformatics results and wants to validate
- User asks which experiments to perform
- User asks about primer design or antibody selection

**Recommended Workflow**: Relevant domain workflow (validation section)

**Output Structure**:
1. Key findings to validate
2. Recommended experiments (prioritized)
3. Experimental design considerations
4. Expected results
5. Controls needed
6. Limitations of validation

**Risks**: Not all bioinformatics predictions can be validated; negative results are also informative.

---

## Multi-Workflow Routing Examples

| User Input | Primary Workflow | Secondary Workflow(s) |
|------------|-----------------|----------------------|
| "I have a count matrix and grouping info, I want to find DEGs" | bulk_rnaseq | — |
| "I have an h5ad file and want to see cell clusters" | single_cell | — |
| "I want to study cuproptosis and autophagy in cervical cancer" | literature_driven_mechanism | bulk_rnaseq + clinical_survival |
| "I have an OTU table and grouping info" | microbiome | — |
| "I want to build a prognostic risk model" | clinical_survival | machine_learning |
| "I have RNA-seq and proteomics data from the same samples" | multiomics_integration | bulk_rnaseq + proteomics |
| "I have 10X Visium spatial data" | spatial_transcriptomics | — |
| "Help me respond to reviewer #2's comment" | reviewer_response | — |
| "I have a metabolite list and want pathway analysis" | metabolomics | — |
| "I want to study a gene's role in tumor immunity" | literature_driven_mechanism | bulk_rnaseq + clinical_survival + single_cell |
| "I have methylation 450K array data" | epigenomics | — |
| "I want to find diagnostic biomarkers using ML" | machine_learning | bulk_rnaseq |
| "I have MAF files from TCGA" | genomics_variant | — |
| "I want to integrate scRNA-seq with spatial data" | spatial_transcriptomics | single_cell + multiomics_integration |
| "Write the Methods for my DESeq2 analysis" | paper_writing | bulk_rnaseq |

## Deferral Override: When NOT to Route

Before selecting any workflow, the router must check for the following triggers. If any match, the Skill should recommend deferral to a human expert **in addition to** (not instead of) offering analytical assistance.

| Trigger | Example User Language | Action |
|---------|----------------------|--------|
| Clinical diagnosis | "diagnose my patient," "what treatment should this patient get" | Flag as clinical; recommend deferral to clinician |
| Identifiable patient data | Data with patient names, MRNs, or dates of service | Stop; warn about privacy; recommend IRB consultation |
| Legal/regulatory | "is this HIPAA compliant," "submit to FDA" | Recommend institutional legal/regulatory resources |
| Wet-lab execution | "run this PCR for me," "perform the staining" | Clarify Skill provides protocols, not execution |
| Grant submission without review | "write my grant proposal and I'll submit" | Provide draft but recommend expert review before submission |

**Important**: The deferral is a safeguard, not a refusal. The Skill should still assist with appropriate analytical portions of the task while clearly flagging aspects that require human expertise.

**Routing rule**: If a deferral trigger is detected, prepend the response with the boundary warning, then proceed with the standard workflow routing below.

---

## Decision Flow

```
User Input
    │
    ├─ Contains data file? ──────────────────── Identify data type → Select workflow
    │
    ├─ Contains error/code? ─────────────────── Code debugging mode
    │
    ├─ Contains figure/screenshot? ──────────── Figure interpretation mode
    │
    ├─ Contains reviewer comments? ──────────── Reviewer response mode
    │
    ├─ Asks for manuscript text? ────────────── Paper writing mode
    │
    ├─ Asks about research design? ──────────── Research design mode
    │
    └─ Describes analysis goal? ─────────────── Match to domain workflow
         │
         ├─ Multiple omics mentioned? ────────── Multi-omics integration
         ├─ Mentions survival/prognosis? ──────── Clinical survival
         ├─ Mentions ML/prediction? ───────────── Machine learning
         └─ Single domain clear? ──────────────── Specific domain workflow
```
