# Bioinformatics Research Skill

A comprehensive, modular bioinformatics research assistant supporting nearly all common bioinformatics scenarios — from bulk RNA-seq and single-cell analysis to spatial transcriptomics, proteomics, metabolomics, microbiome, epigenomics, genomics, survival analysis, machine learning, and multi-omics integration.

## What Is This?

This is a **Skill** built on the [Agent Skills open standard](https://agentskills.io/specification) — compatible with both **Claude Code** and **OpenAI Codex CLI**. It provides a structured knowledge base and prompt system for AI-assisted bioinformatics research:

- **Task recognition**: Automatically identifies what type of bioinformatics analysis the user needs
- **Workflow guidance**: Step-by-step analysis pipelines for 12+ omics domains
- **Code generation**: Runnable R and Python templates with proper environment management
- **Result interpretation**: Structured frameworks for explaining bioinformatics results
- **Manuscript support**: Templates for Methods, Results, Figure legends, Discussion, and reviewer responses

## Installation

### Requirements

- **Claude Code** (latest) or **OpenAI Codex CLI** (latest)
- R >= 4.2.0 (for R-based analyses)
- Python >= 3.9 + uv package manager (for Python-based analyses)

### Install for Claude Code

```bash
# Clone or symlink into the Claude Code skills directory
mkdir -p .claude/skills
ln -s /path/to/universal-bioinformatics-skill .claude/skills/bioinformatics

# Verify
claude --skills
```

The skill name in `SKILL.md` is `bioinformatics`, so the install directory must also be named `bioinformatics`.

### Install for Codex CLI

```bash
# Clone or symlink into the Codex agents directory
mkdir -p .agents/skills
ln -s /path/to/universal-bioinformatics-skill .agents/skills/bioinformatics

# Verify
codex skills list
```

### Install Globally (Available Across All Projects)

```bash
# Claude Code (global)
ln -s /path/to/universal-bioinformatics-skill ~/.claude/skills/bioinformatics

# Codex CLI (global)
ln -s /path/to/universal-bioinformatics-skill ~/.agents/skills/bioinformatics
```

### Using the Same Installation for Both

If you use both Claude Code and Codex, set up once and symlink:

```bash
# Store the skill in a shared location
git clone <this-repo> ~/shared-skills/bioinformatics

# Link into both platforms
mkdir -p ~/.claude/skills ~/.agents/skills
ln -s ~/shared-skills/bioinformatics ~/.claude/skills/bioinformatics
ln -s ~/shared-skills/bioinformatics ~/.agents/skills/bioinformatics
```

### How It Works (Progressive Disclosure)

Both platforms use a 3-level loading strategy to stay within context limits:

| Level | What | When |
|-------|------|------|
| 1 | Skill `name` + `description` (from frontmatter) | Always in context |
| 2 | Full `SKILL.md` body | Loaded when the skill is activated |
| 3 | Referenced files in `references/workflows/`, `references/templates/`, `references/environment/` | Loaded on demand as needed |

## Core Capabilities

| Domain | Key Analyses |
|--------|-------------|
| Bulk RNA-seq | DESeq2, limma, edgeR, GO/KEGG, GSEA, GSVA, WGCNA |
| Single-cell RNA-seq | Seurat, scanpy, clustering, annotation, CellChat, SCENIC, inferCNV |
| Spatial Transcriptomics | Spatial clustering, deconvolution, co-localization, spatial communication |
| Proteomics | Differential proteins, PPI networks, proteome-transcriptome integration |
| Metabolomics | PCA, PLS-DA, OPLS-DA, pathway enrichment |
| Microbiome | Alpha/beta diversity, LEfSe, PICRUSt2, microbiome-metabolome integration |
| Epigenomics | DNA methylation, ATAC-seq, ChIP-seq, DMP/DMR, motif enrichment |
| Genomics | Mutation analysis, TMB, MSI, mutation signatures, maftools, GWAS |
| Clinical Survival | KM, Cox regression, nomogram, timeROC, DCA, risk score |
| Machine Learning | LASSO, SVM-RFE, Random Forest, XGBoost, SHAP, external validation |
| Multi-omics | Cross-omics correlation, network integration, MOFA, pathway consistency |
| Literature-Driven | Hypothesis → public database mining → validation → publication design |

## Applicable Research Fields

- Oncology / Cancer biology
- Immunology
- Cardiovascular medicine
- Metabolic diseases
- Neuroscience
- Plant science / Agriculture
- Microbiology
- Pharmacology / Drug discovery
- Clinical cohort studies
- Animal model research

## Quick Start

### Environment Setup

See `references/environment/` for environment configuration guides:
- `R_environment.md` — R environment and package management
- `Python_uv_environment.md` — Python environment with uv package manager
- `uv_usage.md` — uv quick reference

## Workflow Selection Guide

| If the user says... | Use this workflow |
|---------------------|------------------|
| "I have a count matrix and grouping info" | `bulk_rnaseq_workflow.md` |
| "I have an h5ad file with cell data" | `single_cell_workflow.md` |
| "I have 10X Visium spatial data" | `spatial_transcriptomics_workflow.md` |
| "I have protein expression data" | `proteomics_workflow.md` |
| "I have metabolite abundance data" | `metabolomics_workflow.md` |
| "I have an OTU/ASV table" | `microbiome_workflow.md` |
| "I have methylation array data" | `epigenomics_workflow.md` |
| "I have mutation/MAF data" | `genomics_variant_workflow.md` |
| "I want to do survival analysis" | `clinical_survival_workflow.md` |
| "I want to build a prediction model" | `machine_learning_bioinformatics_workflow.md` |
| "I have multi-omics data" | `multiomics_integration_workflow.md` |
| "I want to study X pathway in Y disease" | `literature_driven_mechanism_workflow.md` |

See `references/workflows/00_task_router.md` for the complete routing logic.

## Example Usage Scenarios

### Scenario 1: TCGA Cancer Study
> "I want to study whether FOXM1 is associated with tumor immunity and prognosis in TCGA."

→ Uses: `literature_driven_mechanism_workflow.md` + `clinical_survival_workflow.md`

### Scenario 2: Single-Cell Analysis
> "I have scRNA-seq data from tumor and normal tissue. I want to identify immune cell subpopulations."

→ Uses: `single_cell_workflow.md`

### Scenario 3: Diagnostic Biomarker Model
> "I want to build a machine learning model to find diagnostic biomarkers."

→ Uses: `machine_learning_bioinformatics_workflow.md`

## Project Structure

```
universal-bioinformatics-skill/
│
├── SKILL.md                              # Core Skill specification
├── README.md                             # This file
├── CITATION.cff                          # Citation metadata
├── LICENSE                               # License file
│
├── references/                           # Reference materials (Agent Skills standard)
│   ├── environment/                      # Environment setup
│   │   ├── R_environment.md             # R + Rlib guide
│   │   ├── Python_uv_environment.md     # Python + uv guide
│   │   └── uv_usage.md                  # uv quick reference
│   │
│   ├── templates/                        # Reusable templates
│   │   ├── research_design_template.md  # Study design template
│   │   ├── data_interpretation_template.md
│   │   ├── figure_explanation_template.md
│   │   ├── paper_writing_template.md
│   │   ├── code_generation_template.md
│   │   ├── reviewer_response_template.md
│   │   ├── output_format_standards.md
│   │   └── upload_handling_guide.md
│   │
│   ├── workflows/                        # Analysis workflows
│   │   ├── 00_task_router.md
│   │   ├── bulk_rnaseq_workflow.md
│   │   ├── single_cell_workflow.md
│   │   ├── spatial_transcriptomics_workflow.md
│   │   ├── proteomics_workflow.md
│   │   ├── metabolomics_workflow.md
│   │   ├── microbiome_workflow.md
│   │   ├── epigenomics_workflow.md
│   │   ├── genomics_variant_workflow.md
│   │   ├── clinical_survival_workflow.md
│   │   ├── machine_learning_bioinformatics_workflow.md
│   │   ├── multiomics_integration_workflow.md
│   │   └── literature_driven_mechanism_workflow.md
│   │
│   └── examples/                         # Complete case studies (9)
│       ├── example_cancer_tcga.md
│       ├── example_clinical_prediction.md
│       ├── example_metabolic_disease.md
│       ├── example_metabolomics.md
│       ├── example_microbiome.md
│       ├── example_multiomics.md
│       ├── example_proteomics.md
│       ├── example_single_cell.md
│       └── example_spatial_transcriptomics.md
│
├── tests/                                # Validation
│   ├── test_skill_structure.md
│   └── test_task_router_cases.md
│
└── logs/                                 # Runtime logs
    └── .gitkeep
```

## FAQ

**Q: Can this Skill handle non-cancer research?**
A: Yes. It supports plant science, microbiology, metabolic diseases, animal models, and more. The task router automatically detects the research context.

**Q: Do I need to install all R packages?**
A: No. Install only what you need. The full package list is in `references/environment/R_environment.md`. Use `install_if_missing()` to install individual packages as needed.

**Q: Can I use pip instead of uv?**
A: Not recommended. uv provides faster installation, reproducible lock files, and better dependency resolution. See `references/environment/Python_uv_environment.md` for details.

**Q: What if a package installation fails?**
A: See `references/environment/R_environment.md` for troubleshooting common issues with system libraries (libcurl, libxml2, hdf5).

**Q: Can I combine multiple workflows?**
A: Yes. The task router supports multi-workflow routing. For example, a study might use `bulk_rnaseq_workflow.md` + `clinical_survival_workflow.md` + `machine_learning_bioinformatics_workflow.md`.

**Q: Can this write manuscripts?**
A: It provides structured templates for Methods, Results, Figure legends, and Discussion sections. Always review and adapt the generated text.

**Q: How do I add a new workflow?**
A: Create a new `.md` file in `references/workflows/` following the standard structure (10 sections). Add routing rules in `00_task_router.md`. Add an example case in `references/examples/`.

## Extending the Skill

### Adding a New Workflow

1. Create `references/workflows/new_workflow.md` with these sections:
   - Applicable Scenarios
   - Input Data Requirements
   - Data Validation Checklist
   - Recommended Analysis Steps
   - Recommended Statistical Methods
   - Recommended Figures
   - Result Interpretation Template
   - Manuscript Writing Template
   - Common Issues and Risks
   - Experimental Validation Suggestions

2. Add routing rules to `references/workflows/00_task_router.md`

3. Add an example case in `references/examples/`
---

*Universal Bioinformatics Research Skill v1.0*
