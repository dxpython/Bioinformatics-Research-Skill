# Proteomics Analysis Workflow

## 1. Applicable Scenarios

- Differential protein expression between conditions
- Protein-protein interaction (PPI) network analysis
- Functional enrichment of differentially expressed proteins
- Integration of proteomics with transcriptomics
- Post-translational modification analysis
- Biomarker discovery at protein level

## 2. Input Data Requirements

| Input | Description | Format |
|-------|-------------|--------|
| Protein expression matrix | Proteins (rows) × Samples (columns) | CSV/TSV (from MaxQuant, Proteome Discoverer, DIA-NN) |
| Sample information | Sample IDs, group labels | CSV/TSV |
| Protein annotation (optional) | UniProt IDs, gene symbols | CSV/TSV |

## 3. Data Validation Checklist

- [ ] Check protein identification confidence (FDR < 1%)
- [ ] Assess missing value patterns (MCAR, MAR, MNAR)
- [ ] Check for contaminants and reverse sequences (remove)
- [ ] Verify sample-to-group mapping

## 4. Recommended Analysis Steps

### Step 1: Data Preprocessing
1. Remove contaminants, reverse hits, proteins identified by site only
2. Log2 transform intensity values
3. Assess missing value pattern

### Step 2: Missing Value Handling

| Strategy | When to Use |
|----------|-------------|
| Filter | Remove proteins with >50% missing across all samples |
| MinProb imputation | Missing Not At Random (below detection limit) |
| KNN imputation | Missing At Random |
| Mixed imputation | Combination (DEP package) |

### Step 3: Normalization
- Median centering
- Quantile normalization
- Variance stabilization (vsn)

### Step 4: Differential Protein Analysis
- limma (moderated t-test): recommended for most proteomics data
- Thresholds: |log2FC| > 1, adjusted p < 0.05
- Volcano plot and heatmap visualization

### Step 5: PPI Network Analysis
- Submit DEPs to STRING database
- Construct PPI network
- Identify hub proteins (degree, betweenness centrality)
- Module detection (MCODE, MCL)

### Step 6: Functional Enrichment
- GO (BP/MF/CC)
- KEGG pathways
- Reactome pathways

### Step 7: Proteomics-Transcriptomics Integration
- Correlate protein and mRNA levels for matched genes
- Four-quadrant plot (protein FC vs. mRNA FC)
- Identify post-transcriptionally regulated genes
- Joint pathway analysis

## 5. Recommended Statistical Methods

| Analysis | Recommended Method | Alternative |
|----------|-------------------|-------------|
| Missing value test | Little's MCAR test | Visual inspection of missingness patterns |
| Normalization | vsn (variance stabilization), quantile normalization | median polish |
| Differential proteins | limma with empirical Bayes | t-test, Wilcoxon (for small samples) |
| Multiple testing | Benjamini-Hochberg FDR | Storey's q-value |
| Enrichment | Hypergeometric test (clusterProfiler) | Fisher's exact test |
| Correlation (protein-mRNA) | Spearman correlation | Pearson (if normal) |
| PPI enrichment | STRINGdb enrichment p-value | Random network permutation |

- For differential analysis: limma is preferred due to robust variance estimation with small sample sizes common in proteomics
- Missing values: characterize the pattern (MCAR/MAR/MNAR) before imputation
- PPI network: always test for enrichment against a random network of same size
- For proteomics-transcriptomics integration: report both correlation coefficients and p-values with FDR correction

## 6. Recommended Figures

| Figure | Description |
|--------|-------------|
| Missing value heatmap | Pattern of missingness |
| PCA plot | Sample separation |
| Volcano plot | Differential proteins |
| Heatmap | Top differential proteins |
| PPI network | STRING interaction network |
| GO/KEGG plots | Enrichment results |
| Correlation scatter | Protein vs. mRNA expression |
| Four-quadrant plot | Protein FC vs. mRNA FC |

## 7. Result Interpretation Template

> Quantitative proteomics analysis identified [N] differentially expressed proteins (|log2FC| > 1, adjusted p < 0.05), including [X] upregulated and [Y] downregulated proteins in [Condition1] compared to [Condition2]. PPI network analysis revealed a core module containing [hub proteins], associated with [biological function]. GO analysis showed enrichment in [processes]. Integration with transcriptomic data showed [X]% concordance (same direction of change), with [Y] proteins showing post-transcriptional regulation.

## 8. Manuscript Writing Template

### Methods

> **Mass spectrometry-based proteomics analysis**
>
> Protein abundance data were obtained by [TMT / iTRAQ / label-free] mass spectrometry on a [instrument] platform. Raw spectra were searched against the [UniProt / SwissProt] database using [MaxQuant (vX.X) / Proteome Discoverer (vX.X)]. Protein quantification was based on [peptides with FDR < 1% / unique + razor peptides]. Proteins with > [X]% missing values across all samples were excluded. Remaining missing values were imputed using the [MinProb / KNN] method (impute R package, vX.X). Data were variance-stabilized and normalized using the vsn package (vX.X). Differential protein analysis was performed using limma (vX.X). Proteins with |log2FC| > [threshold] and BH-adjusted p-value < [threshold] were considered significantly changed. Protein-protein interaction networks were constructed using STRING (vX.X, confidence score > 0.4). GO and KEGG enrichment were performed using clusterProfiler (vX.X). For proteomics-transcriptomics integration, protein and mRNA fold changes for matched genes were correlated using Spearman's rank correlation. Post-transcriptionally regulated candidates were identified from the discordant quadrants of the four-quadrant plot.

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| High missing values (>50%) | Use stringent filtering + appropriate imputation |
| Batch effects | Include batch in limma design; use ComBat |
| Low overlap with transcriptomics | Expected; mRNA and protein often poorly correlated |
| Protein ID ambiguity | Use leading razor protein; map to gene symbols |

## 10. Experimental Validation Suggestions

1. **Western blot**: Validate key differential proteins
2. **ELISA**: Quantitative validation for secreted proteins
3. **IHC/IF**: Spatial protein expression
4. **Targeted proteomics (PRM/MRM)**: Precise quantification
