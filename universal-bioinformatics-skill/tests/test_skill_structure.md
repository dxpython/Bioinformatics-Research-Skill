# Project structure completeness checklist

## 1. Purpose

Verify that the `universal-bioinformatics-skill` directory layout is complete and expected files exist.

## 2. Directory overview

```
universal-bioinformatics-skill/
├── SKILL.md              # Core skill specification
├── README.md             # Project overview
├── references/           # Reference materials (Agent Skills standard)
│   ├── workflows/        # Workflow definitions (13)
│   ├── templates/        # Writing / interpretation templates (8)
│   ├── environment/      # Environment setup guides
│   └── examples/         # End-to-end case studies (9)
├── tests/                # Structure and router checks
└── logs/                 # Logs (e.g. installs)
```

## 3. Checklist

### 3.1 Root files

| # | Path | Status | Notes |
|---|------|--------|-------|
| 1 | `SKILL.md` | [ ] | Skill entry, core specification |
| 2 | `README.md` | [ ] | Project description |

### 3.2 references/workflows/

| # | Path | Status | Notes |
|---|------|--------|-------|
| 1 | `references/workflows/00_task_router.md` | [ ] | Task router |
| 2 | `references/workflows/bulk_rnaseq_workflow.md` | [ ] | Bulk RNA-seq |
| 3 | `references/workflows/single_cell_workflow.md` | [ ] | Single-cell RNA-seq |
| 4 | `references/workflows/spatial_transcriptomics_workflow.md` | [ ] | Spatial transcriptomics |
| 5 | `references/workflows/proteomics_workflow.md` | [ ] | Proteomics |
| 6 | `references/workflows/metabolomics_workflow.md` | [ ] | Metabolomics |
| 7 | `references/workflows/microbiome_workflow.md` | [ ] | Microbiome |
| 8 | `references/workflows/epigenomics_workflow.md` | [ ] | Epigenomics |
| 9 | `references/workflows/genomics_variant_workflow.md` | [ ] | Genomic variants |
| 10 | `references/workflows/clinical_survival_workflow.md` | [ ] | Clinical survival |
| 11 | `references/workflows/machine_learning_bioinformatics_workflow.md` | [ ] | ML biomarker modeling |
| 12 | `references/workflows/multiomics_integration_workflow.md` | [ ] | Multi-omics integration |
| 13 | `references/workflows/literature_driven_mechanism_workflow.md` | [ ] | Literature-driven mechanism |

### 3.3 references/templates/

| # | Path | Status | Notes |
|---|------|--------|-------|
| 1 | `references/templates/research_design_template.md` | [ ] | Study design |
| 2 | `references/templates/data_interpretation_template.md` | [ ] | Data interpretation |
| 3 | `references/templates/figure_explanation_template.md` | [ ] | Figure interpretation |
| 4 | `references/templates/paper_writing_template.md` | [ ] | Manuscript drafting |
| 5 | `references/templates/code_generation_template.md` | [ ] | Code generation |
| 6 | `references/templates/reviewer_response_template.md` | [ ] | Peer review responses |
| 7 | `references/templates/output_format_standards.md` | [ ] | Output format standards |
| 8 | `references/templates/upload_handling_guide.md` | [ ] | Upload handling guide |

### 3.4 references/environment/

| # | Path | Status | Notes |
|---|------|--------|-------|
| 1 | `references/environment/R_environment.md` | [ ] | R + Rlib |
| 2 | `references/environment/Python_uv_environment.md` | [ ] | Python + uv |
| 3 | `references/environment/uv_usage.md` | [ ] | uv quick reference |

### 3.5 references/examples/

| # | Path | Status | Notes |
|---|------|--------|-------|
| 1 | `references/examples/example_cancer_tcga.md` | [ ] | TCGA-style oncology |
| 2 | `references/examples/example_metabolic_disease.md` | [ ] | Metabolic disease |
| 3 | `references/examples/example_single_cell.md` | [ ] | scRNA-seq |
| 4 | `references/examples/example_spatial_transcriptomics.md` | [ ] | Spatial |
| 5 | `references/examples/example_microbiome.md` | [ ] | Microbiome |
| 6 | `references/examples/example_multiomics.md` | [ ] | Multi-omics |
| 7 | `references/examples/example_clinical_prediction.md` | [ ] | Clinical prediction |
| 8 | `references/examples/example_proteomics.md` | [ ] | Proteomics |
| 9 | `references/examples/example_metabolomics.md` | [ ] | Metabolomics |

### 3.6 tests/

| # | Path | Status | Notes |
|---|------|--------|-------|
| 1 | `tests/test_skill_structure.md` | [ ] | This checklist |
| 2 | `tests/test_task_router_cases.md` | [ ] | Router test cases |

### 3.7 logs/

| # | Path | Status | Notes |
|---|------|--------|-------|
| 1 | `logs/` | [ ] | Directory exists |

## 4. Example shell check script

```bash
#!/bin/bash
# Verify project structure

BASE_DIR="universal-bioinformatics-skill"
PASS=0
FAIL=0

check_file() {
  if [ -f "$BASE_DIR/$1" ]; then
    echo "[PASS] $1"
    ((PASS++))
  else
    echo "[FAIL] $1 — missing file"
    ((FAIL++))
  fi
}

check_dir() {
  if [ -d "$BASE_DIR/$1" ]; then
    echo "[PASS] $1/"
    ((PASS++))
  else
    echo "[FAIL] $1/ — missing directory"
    ((FAIL++))
  fi
}

echo "========================================="
echo "Structure check"
echo "========================================="

echo ""
echo "--- Directories ---"
check_dir "references"
check_dir "references/workflows"
check_dir "references/templates"
check_dir "references/environment"
check_dir "references/examples"
check_dir "tests"
check_dir "logs"

echo ""
echo "--- Root ---"
check_file "SKILL.md"
check_file "README.md"

echo ""
echo "--- Workflows ---"
for f in 00_task_router bulk_rnaseq_workflow single_cell_workflow spatial_transcriptomics_workflow proteomics_workflow metabolomics_workflow microbiome_workflow epigenomics_workflow genomics_variant_workflow clinical_survival_workflow machine_learning_bioinformatics_workflow multiomics_integration_workflow literature_driven_mechanism_workflow; do
  check_file "references/workflows/${f}.md"
done

echo ""
echo "--- Templates ---"
for f in research_design_template data_interpretation_template figure_explanation_template paper_writing_template code_generation_template reviewer_response_template output_format_standards upload_handling_guide; do
  check_file "references/templates/${f}.md"
done

echo ""
echo "--- Environment ---"
for f in R_environment Python_uv_environment uv_usage; do
  check_file "references/environment/${f}.md"
done

echo ""
echo "--- Examples ---"
for f in example_cancer_tcga example_metabolic_disease example_single_cell example_spatial_transcriptomics example_microbiome example_multiomics example_clinical_prediction example_proteomics example_metabolomics; do
  check_file "references/examples/${f}.md"
done

echo ""
echo "--- Tests ---"
check_file "tests/test_skill_structure.md"
check_file "tests/test_task_router_cases.md"

echo ""
echo "========================================="
echo "Result: PASS=$PASS, FAIL=$FAIL"
echo "Total checks: $((PASS + FAIL))"
echo "Pass rate: compute as PASS / (PASS+FAIL) if needed"
echo "========================================="
```

## 5. Acceptance criteria

### P0 (must)

- Required directories exist
- `SKILL.md` exists with full specification
- All workflow files exist and are non-empty
- All template files are complete

### P1 (should)

- All example files are complete
- Environment documentation is present

### P2 (optional)

- `README.md` present
- `logs/` directory exists
