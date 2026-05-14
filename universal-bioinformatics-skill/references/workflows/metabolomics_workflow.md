# Metabolomics Analysis Workflow

## 1. Applicable Scenarios

- Differential metabolite identification between groups
- Metabolic pathway activity analysis
- Biomarker discovery for disease diagnosis
- Metabolite-gene correlation analysis
- Integration with transcriptomics for metabolic mechanism studies
- Drug metabolism and pharmacometabolomics

## 2. Input Data Requirements

| Input | Description | Format |
|-------|-------------|--------|
| Metabolite abundance matrix | Metabolites (rows) × Samples (columns) | CSV/TSV/XLSX |
| Sample information | Sample IDs, group labels, covariates | CSV/TSV |
| Metabolite annotation | Metabolite IDs (HMDB, KEGG compound), names, m/z, RT | CSV/TSV |
| QC samples (optional) | Pooled QC samples for quality assessment | Included in matrix |

## 3. Data Validation Checklist

- [ ] Check for QC sample consistency (RSD of QC < 30%)
- [ ] Assess missing value pattern
- [ ] Identify and handle outlier samples
- [ ] Verify metabolite identification confidence levels

## 4. Recommended Analysis Steps

### Step 1: Quality Control
1. Evaluate QC sample clustering (should cluster tightly in PCA)
2. Filter metabolites with high RSD in QC samples (>30%)
3. Filter metabolites with excessive missing values (>50%)

### Step 2: Missing Value Imputation
- Replace with half of minimum value (common for metabolomics)
- KNN imputation
- Minimum value imputation per metabolite

### Step 3: Normalization
- Log2 transformation
- Pareto scaling (recommended for metabolomics)
- Auto scaling (mean-centered, unit variance)
- Total ion current (TIC) normalization

### Step 4: Multivariate Statistical Analysis

| Method | Purpose | Supervised? |
|--------|---------|-------------|
| PCA | Overview, outlier detection, QC assessment | No |
| PLS-DA | Discrimination between groups | Yes |
| OPLS-DA | Maximize group separation, identify key metabolites | Yes |

**PLS-DA/OPLS-DA validation:**
- Cross-validation (Q2 > 0.5 is good)
- Permutation test (200 permutations, p < 0.05)
- R2/Q2 values must be reported

### Step 5: Differential Metabolite Identification
- VIP score from (O)PLS-DA > 1
- Univariate test: t-test (normal) or Wilcoxon (non-normal)
- Fold change threshold (typically FC > 1.5 or FC > 2)
- Combined criteria: VIP > 1 AND p < 0.05 AND FC > 1.5

### Step 6: Metabolic Pathway Analysis
- KEGG metabolic pathway enrichment
- MetaboAnalyst pathway analysis (topology-based)
- SMPDB (Small Molecule Pathway Database)

### Step 7: Metabolite-Gene Correlation
- Correlate differential metabolites with gene expression
- Spearman correlation (non-parametric)
- Network visualization of metabolite-gene relationships
- Identify metabolic enzyme genes

## 5. Recommended Statistical Methods

| Analysis | Recommended Method | Alternative |
|----------|-------------------|-------------|
| QC assessment | RSD calculation across QC samples | PCA of QC samples |
| Normalization | Pareto scaling, auto scaling | log transformation, TIC |
| Multivariate analysis | PLS-DA / OPLS-DA (ropls, MetaboAnalystR) | PCA only for initial exploration |
| Model validation | Permutation test (n ≥ 1000), 7-fold CV | External test set validation |
| Differential metabolites | VIP > 1 AND p < 0.05 AND |FC| > 1.5 | Wilcoxon + VIP, t-test + VIP |
| Multiple testing | Benjamini-Hochberg FDR | Bonferroni (conservative) |
| Pathway enrichment | Hypergeometric test (MetaboAnalyst, KEGG) | Mummichog (for untargeted) |
| Metabolite-gene correlation | Spearman correlation | Kendall's tau (robust to outliers) |

- OPLS-DA requires validation: Q² > 0.5 and significant permutation test (p < 0.05)
- For untargeted metabolomics: report metabolite identification confidence level (Level 1-4 per MSI)
- Pathway analysis: account for the fact that not all metabolites in a pathway are detected
- Metabolite-gene integration: interpret correlation direction in the context of enzyme activity

## 6. Recommended Figures

| Figure | Description |
|--------|-------------|
| PCA score plot | Sample overview and QC check |
| PLS-DA score plot | Group discrimination |
| OPLS-DA score plot + S-plot | Key discriminating metabolites |
| Permutation test plot | Model validation |
| Volcano plot | Differential metabolites |
| Heatmap | Top differential metabolites |
| Pathway bubble plot | Enriched metabolic pathways |
| Correlation network | Metabolite-gene associations |
| Box plots | Individual metabolite abundances |

## 7. Result Interpretation Template

> Metabolomics profiling identified [N] metabolites across [M] samples. PCA showed clear separation between [Group1] and [Group2] (Figure XA). PLS-DA model achieved good discrimination (R2Y = X.XX, Q2 = X.XX), validated by permutation test (p < 0.05, Figure XB). A total of [N] differential metabolites were identified (VIP > 1, p < 0.05, FC > 1.5), including [X] increased and [Y] decreased metabolites (Figure XC). KEGG pathway analysis revealed significant enrichment in [pathway names] (Figure XD).

## 8. Manuscript Writing Template

### Methods

> **[Untargeted / Targeted] metabolomics analysis**
>
> Metabolic profiling was performed using [LC-MS / GC-MS] on a [instrument, e.g., Q-Exactive HF] coupled to a [column, e.g., ACQUITY UPLC BEH C18] column. Data preprocessing (peak detection, alignment, and annotation) was performed using [XCMS (vX.X) / Compound Discoverer (vX.X) / Progenesis QI (vX.X)]. Metabolites with relative standard deviation (RSD) > 30% in QC samples were excluded. Missing values were imputed using the [k-nearest neighbor / minimum value / half-minimum] method. Data were Pareto-scaled before multivariate analysis. Principal component analysis (PCA) was used to visualize overview trends. Partial least squares discriminant analysis (PLS-DA) and orthogonal PLS-DA (OPLS-DA) were performed using the ropls package (vX.X) in R or MetaboAnalystR (vX.X). Model validity was assessed by permutation testing (n = 1000) and 7-fold cross-validation. Differentially abundant metabolites were selected using VIP > 1, Wilcoxon rank-sum test p < 0.05, and |fold change| > 1.5. Pathway enrichment was performed against KEGG and SMPDB databases using MetaboAnalystR or the web-based MetaboAnalyst platform. Metabolite-gene correlations were assessed using Spearman's rank correlation with BH adjustment.

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| Poor QC clustering | Re-evaluate data quality; normalize by QC |
| Overfitting PLS-DA | Must validate with permutation test; use cross-validation |
| Metabolite misidentification | Report confidence levels (Level 1-4); validate with standards |
| Batch effects | Include QC-based correction (QC-RLSC) |
| Small sample size | Use non-parametric tests; avoid complex models |

## 10. Experimental Validation Suggestions

1. **Targeted metabolomics**: Validate differential metabolites with standards
2. **Enzymatic assays**: Measure key metabolic enzyme activities
3. **Isotope tracing**: Track metabolic flux in cell models
4. **Metabolic gene expression**: Correlate with qRT-PCR of metabolic enzymes
