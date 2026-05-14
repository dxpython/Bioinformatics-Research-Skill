# Literature-Driven Mechanism Analysis Workflow

## 1. Applicable Scenarios

- Investigating a specific pathway/mechanism in a disease using public databases
- "Disease + Pathway + Gene Set" hypothesis-driven mining
- Immune microenvironment characterization for a gene/pathway of interest
- Validating a molecular mechanism across multiple datasets
- Designing a complete publication-ready bioinformatics study from a hypothesis

## 2. Input Data Requirements

| Input | Description | Source |
|-------|-------------|--------|
| Research hypothesis | Disease + Mechanism/Pathway + Direction | User-defined |
| Gene set of interest | Pathway genes, signature genes | MSigDB, KEGG, literature |
| Public expression data | Bulk RNA-seq or microarray | TCGA, GEO, ICGC |
| Clinical data | Survival, stage, grade | TCGA, GEO |
| Single-cell data (optional) | For cell-type-level validation | GEO, HCA |

## 3. Data Validation Checklist

Before starting the mechanism-driven analysis, verify:
```
□ Public datasets are accessible and not withdrawn
□ Expression data are normalized (check distribution per sample)
□ Clinical annotation is complete (OS_time, OS_status, key covariates)
□ Gene symbols in the gene set are current (check against HGNC)
□ Sample size per group is adequate (typically > 10 per group)
□ Batch effects between datasets (if using multi-cohort) are assessed
□ Tissue type matches the research question (e.g., tumor vs. tumor, not tumor vs. normal for tumor-specific analysis)
□ Platform types are identified and cross-platform normalization planned
```

| Issue | Detection | Action |
|-------|-----------|--------|
| Gene set contains outdated symbols | Cross-reference HGNC | Update to current approved symbols |
| Small sample size | n < 20 per group | Acknowledge limitation; avoid overfitting |
| Clinical data incomplete | Missing survival or stage | Exclude samples; note as limitation |
| Cross-platform batch effects | PCA by dataset | Use ComBat; validate findings per dataset |
| Unbalanced groups | Unequal sample sizes per group | Use balanced resampling; weighted statistics |

## 4. Recommended Analysis Steps

### Step 1: Hypothesis Formulation and Gene Set Definition
1. Define the core hypothesis (e.g., "Ferroptosis-related genes are associated with prognosis and immune microenvironment in hepatocellular carcinoma")
2. Obtain gene sets: MSigDB hallmark/C2/C5, KEGG, Reactome, or curated from literature
3. Verify gene set quality (size, specificity, overlap with other sets)

### Step 2: Expression-Level Analysis
1. Download expression data (TCGA via TCGAbiolinks, GEO via GEOquery)
2. Calculate pathway activity score per sample:
   - **GSVA**: Gene Set Variation Analysis (per-sample enrichment score)
   - **ssGSEA**: Single-sample GSEA
   - **AUCell**: For single-cell data
   - **Average z-score**: Simple mean of standardized gene expression
3. Stratify patients by score (median or optimal cutpoint)
4. Compare expression of pathway genes between disease/normal groups

### Step 3: Clinical Correlation
1. Correlate pathway score with clinical variables (stage, grade, age, sex)
2. Survival analysis: KM curves + log-rank (high vs. low pathway score)
3. Cox regression: univariate and multivariate (adjusting for clinical covariates)
4. Subgroup analysis: stratify by stage, grade, or molecular subtype

### Step 4: GSEA Validation
1. Rank genes by pathway score correlation (Spearman r)
2. Run GSEA against hallmark, KEGG, Reactome gene sets
3. Identify co-activated and co-suppressed pathways
4. Validate pathway connections to the core mechanism

### Step 5: Immune Microenvironment Analysis
1. Estimate immune cell composition:
   - CIBERSORT (22 immune cell types)
   - MCPcounter
   - ESTIMATE (immune/stromal scores)
   - xCell, TIMER
2. Correlate pathway score with immune cell fractions
3. Compare immune checkpoint expression between high/low score groups
4. Assess immunotherapy response prediction (if applicable)

### Step 6: Key Gene Identification
1. Identify individual genes in the pathway with strongest clinical association
2. Univariate Cox regression for each pathway gene
3. LASSO or multivariate Cox for gene signature
4. Build risk score from selected genes

### Step 7: Multi-Dataset Validation
1. Validate findings in 2-3 independent GEO datasets
2. Check consistency of key gene expression patterns
3. Validate survival associations in external cohorts
4. Meta-analysis if sufficient datasets available

### Step 8: Single-Cell Level Validation (if available)
1. Locate pathway gene expression in specific cell types
2. Compare pathway activity between cell types using AUCell
3. Identify which cell population drives the bulk-level signal
4. Cell communication involving pathway genes

### Step 9: Experimental Validation Roadmap
1. Propose wet-lab experiments to validate key findings
2. Priority: qRT-PCR → WB → IHC → functional assays
3. Design knockdown/overexpression experiments for hub genes
4. Suggest animal model experiments if applicable

## 5. Recommended Statistical Methods

| Analysis | Method | Key Threshold |
|----------|--------|---------------|
| DEG | DESeq2 / limma | |log2FC| > 1, FDR < 0.05 |
| Pathway activity scoring | GSVA / ssGSEA | Normalized enrichment score; compare via t-test/Wilcoxon |
| Survival (univariate) | Kaplan-Meier + log-rank test | Cutoff: median, optimal, or tertile |
| Survival (multivariate) | Cox proportional hazards | Test PH assumption with Schoenfeld residuals |
| Risk score | LASSO-Cox (glmnet) | lambda.min or lambda.1se |
| Immune infiltration | CIBERSORT / ssGSEA / ESTIMATE | Spearman correlation with target |
| Gene set acquisition | Literature search + MSigDB + KEGG | Use curated, published gene sets only |
| Multi-dataset validation | Meta-analysis of HR (metafor) | Forest plot; test heterogeneity with I² |
| Correlation | Spearman (default), Pearson (if normal) | Report rho and FDR |

- GSVA/ssGSEA: validate pathway activity findings in at least 2 independent datasets
- Cox regression: always test proportional hazards assumption; report HR with 95% CI
- LASSO: stability of selected features should be assessed (e.g., bootstrap selection frequency)
- Multi-dataset: use meta-analysis rather than simple overlap for external validation evidence

## 6. Recommended Figures (Complete Figure Set for Publication)

| Figure | Content |
|--------|---------|
| Fig 1A | Study design flowchart |
| Fig 1B | Pathway gene expression heatmap (tumor vs. normal) |
| Fig 1C | GSVA score comparison between groups |
| Fig 2A | KM survival curve (high vs. low pathway score) |
| Fig 2B | Univariate Cox forest plot of pathway genes |
| Fig 2C | Multivariate Cox forest plot |
| Fig 3A | GSEA enrichment plots (top hallmark pathways) |
| Fig 3B | GSVA heatmap of co-activated pathways |
| Fig 4A | Immune cell correlation heatmap |
| Fig 4B | Immune checkpoint expression comparison |
| Fig 4C | ESTIMATE scores comparison |
| Fig 5A | LASSO coefficient path |
| Fig 5B | Risk score KM curve (training set) |
| Fig 5C | Risk score KM curve (validation set) |
| Fig 5D | Time-dependent ROC |
| Fig 6 | Nomogram + calibration curve |
| Fig 7 | Single-cell UMAP with pathway score |
| Supp | External validation, subgroup analyses |

## 7. Result Interpretation Template

> To investigate the role of [pathway] in [disease], we calculated [pathway] activity scores using GSVA across [N] [disease] samples from TCGA. The [pathway] score was significantly elevated in tumor compared to normal tissue (p = X.XX, Figure 1C). Patients with high [pathway] scores showed significantly worse overall survival (HR = X.XX, 95% CI: X.XX-X.XX, log-rank p = X.XX, Figure 2A). GSEA analysis revealed that high [pathway] activity was associated with [hallmark pathways], suggesting [mechanistic interpretation] (Figure 3A). Immune microenvironment analysis showed that high [pathway] scores correlated with increased infiltration of [immune cell types] (r = X.XX, p = X.XX) and elevated expression of immune checkpoints including PD-L1 and CTLA-4 (Figure 4). These findings were validated in [N] independent GEO cohorts (GSE XXXXX, GSE XXXXX).

## 8. Manuscript Writing Template

### Title Pattern
> [Pathway/Mechanism]-related [gene set/signature] predicts prognosis and [immune microenvironment/treatment response] in [disease]: a comprehensive bioinformatics analysis

### Methods
> [Pathway] gene set was obtained from [MSigDB/KEGG/literature]. Gene Set Variation Analysis (GSVA) was used to calculate per-sample pathway activity scores. Patients were stratified into high and low score groups based on the median value. Immune cell composition was estimated using CIBERSORT. Prognostic significance was assessed using Kaplan-Meier analysis and Cox proportional hazards regression. Findings were validated in independent cohorts from GEO.

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| Gene set too large/small | Optimal: 15-200 genes; too large reduces specificity |
| Multiple pathway overlaps | Check overlap; use unique genes if needed |
| TCGA bias | Validate in non-TCGA datasets |
| Immune estimation inconsistency | Use multiple methods; report concordant results |
| Over-interpreting GSVA scores | Scores are relative, not absolute pathway activity |
| No mechanism validation | State as hypothesis-generating; propose experiments |
| Circular reasoning | Don't use the same gene set for scoring AND enrichment |

## 10. Experimental Validation Suggestions

1. **qRT-PCR**: Validate key gene expression in clinical samples
2. **Western blot**: Protein-level validation
3. **IHC**: Spatial expression in tissue microarray
4. **Functional assays**: siRNA/shRNA knockdown of hub genes → proliferation, migration, apoptosis
5. **Flow cytometry**: Validate predicted immune cell composition changes
6. **Animal models**: Test gene/pathway manipulation in vivo
