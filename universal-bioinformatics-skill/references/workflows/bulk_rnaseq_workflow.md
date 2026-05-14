# Bulk RNA-seq Analysis Workflow

## 1. Applicable Scenarios

- Differential gene expression analysis between groups (e.g., disease vs. control, treatment vs. vehicle)
- Gene expression profiling across conditions
- Identification of biomarkers or pathway-level changes
- Integration with clinical data for prognostic studies
- WGCNA for co-expression network analysis
- GSEA/GSVA for pathway-level inference

## 2. Input Data Requirements

| Input | Description | Format |
|-------|-------------|--------|
| Expression matrix | Genes (rows) × Samples (columns) | CSV/TSV/TXT |
| Sample information | Sample IDs, group labels, covariates | CSV/TSV |
| Gene annotation (optional) | Gene ID to symbol mapping | CSV/TSV |

### Expression Data Types

| Type | Description | Suitable For |
|------|-------------|-------------|
| **Raw counts** | Integer counts from featureCounts/HTSeq | DESeq2, edgeR (preferred) |
| **TPM** | Transcripts Per Million | Cross-sample comparison, GSVA, WGCNA |
| **FPKM/RPKM** | Fragments/Reads Per Kilobase per Million | Visualization only (not recommended for DE) |

**Critical**: DESeq2 and edgeR require **raw counts**. Do NOT use TPM/FPKM for these tools.

## 3. Data Validation Checklist

- [ ] Expression matrix has gene IDs as row names and sample IDs as column names
- [ ] Sample info sample IDs match expression matrix column names
- [ ] No duplicated gene IDs (if present, aggregate or filter)
- [ ] Count data contains non-negative integers (for DESeq2/edgeR)
- [ ] At least 2 replicates per group (3+ recommended)
- [ ] Check for outlier samples (PCA)
- [ ] Remove low-count genes (e.g., keep genes with >10 counts in ≥3 samples)

## 4. Recommended Analysis Steps

### Step 1: Data Preprocessing
1. Read expression matrix and sample information
2. Match sample IDs between expression and metadata
3. Filter low-expression genes
4. Assess data quality (library sizes, gene detection rates)

### Step 2: Exploratory Analysis
1. **PCA**: Visualize sample separation and detect batch effects
2. **Sample correlation heatmap**: Assess sample similarity
3. **Boxplot of library sizes**: Check for technical outliers

### Step 3: Differential Expression Analysis

#### Choosing the Right Tool

| Criterion | DESeq2 | edgeR | limma-voom |
|-----------|--------|-------|------------|
| Input | Raw counts | Raw counts | Raw counts (voom transforms) |
| Best for | Small-medium sample sizes | Small sample sizes | Large sample sizes (>20/group) |
| Normalization | Median-of-ratios | TMM | Quantile or TMM+voom |
| Dispersion | Shrinkage estimator | Tagwise/trended | Empirical Bayes |
| Speed | Moderate | Fast | Fast |

**Default recommendation**: DESeq2 for most studies. Use limma-voom for large sample sizes (>50 total).

### Step 4: Result Filtering

Standard thresholds (adjustable):
- |log2FoldChange| > 1 (2-fold change)
- Adjusted p-value < 0.05 (BH correction)

### Step 5: Visualization
1. **Volcano plot**: log2FC vs. -log10(p-value)
2. **MA plot**: log2FC vs. mean expression
3. **Heatmap**: Top DEGs across all samples (with clustering)
4. **Box/violin plots**: Individual gene expression by group

### Step 6: Functional Enrichment
1. **GO analysis**: Biological Process, Molecular Function, Cellular Component
2. **KEGG pathway analysis**: Metabolic and signaling pathways
3. **Reactome**: Reaction-level pathway analysis
4. **GSEA**: Ranked gene list enrichment (does not require DE cutoff)
5. **GSVA**: Per-sample pathway activity scores

### Step 7: Advanced Analyses (Optional)
1. **WGCNA**: Co-expression modules and hub genes
2. **PPI network**: Protein-protein interaction network of DEGs
3. **Transcription factor analysis**: Upstream regulators
4. **Immune cell estimation**: CIBERSORT, MCPcounter, ESTIMATE

### Step 8: External Validation
1. Validate DEGs in independent datasets (e.g., GEO)
2. Check consistency of enrichment results
3. Correlate key genes with clinical outcomes

## 5. Recommended Statistical Methods

| Analysis | Method | R Package |
|----------|--------|-----------|
| Differential expression | Negative binomial GLM | DESeq2 / edgeR |
| Multiple testing correction | Benjamini-Hochberg | Built into DESeq2/edgeR |
| Enrichment significance | Hypergeometric test | clusterProfiler |
| GSEA | Permutation-based | clusterProfiler |
| GSVA | Parametric kernel estimation | GSVA |
| Correlation | Pearson/Spearman | stats |

## 6. Recommended Figures

| Figure | Tool | Description |
|--------|------|-------------|
| PCA plot | ggplot2 | Sample-level dimensionality reduction |
| Volcano plot | ggplot2 + ggrepel | DEG overview with gene labels |
| Heatmap | ComplexHeatmap / pheatmap | Top DEG expression patterns |
| GO dotplot | enrichplot | Enrichment results |
| KEGG barplot | enrichplot | Pathway enrichment |
| GSEA enrichment plot | enrichplot | Running enrichment score |
| GSVA heatmap | ComplexHeatmap | Per-sample pathway scores |
| WGCNA module-trait | WGCNA | Module-phenotype correlations |
| Network plot | igraph / ggraph | PPI or co-expression network |


## 7. Result Interpretation Template

### DEG Results
> Differential expression analysis identified **X upregulated** and **Y downregulated** genes (|log2FC| > 1, adjusted p < 0.05) between [Group1] and [Group2]. The top upregulated genes included [Gene1] (log2FC = X.XX, padj = X.XXe-XX), [Gene2], and [Gene3]. The top downregulated genes included [Gene4], [Gene5], and [Gene6].

### Enrichment Results
> Gene Ontology analysis of upregulated genes revealed significant enrichment in biological processes related to [Process1] (padj = X.XX), [Process2], and [Process3]. KEGG pathway analysis further identified [Pathway1] and [Pathway2] as significantly enriched pathways.

### GSEA Results
> Gene Set Enrichment Analysis showed that [gene set name] was significantly enriched in [Group1] compared to [Group2] (NES = X.XX, padj = X.XX), suggesting [biological interpretation].

## 8. Manuscript Writing Template

### Methods
> Raw count data were obtained from [source]. Low-count genes (fewer than 10 counts across all samples) were removed. Differential expression analysis was performed using DESeq2 (v1.XX.X) with a design formula of ~[condition]. Genes with |log2FoldChange| > 1 and Benjamini-Hochberg adjusted p-value < 0.05 were considered significantly differentially expressed. Gene Ontology and KEGG pathway enrichment analyses were performed using clusterProfiler (v4.XX.X). Gene Set Enrichment Analysis (GSEA) was conducted using the fgsea algorithm implemented in clusterProfiler, with gene sets from MSigDB (v7.X).

### Results
> To identify transcriptomic changes associated with [condition], we performed differential expression analysis on [N] samples ([n1] [Group1] vs. [n2] [Group2]). Principal component analysis (PCA) revealed clear separation between groups along PC1 (Figure XA). A total of [N] differentially expressed genes (DEGs) were identified, including [X] upregulated and [Y] downregulated genes (Figure XB). Functional enrichment analysis of upregulated DEGs showed significant enrichment in [top pathways] (Figure XC).

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| Using TPM/FPKM with DESeq2 | Must use raw counts; re-quantify if needed |
| Too few DEGs | Relax thresholds or check for batch effects |
| Too many DEGs | Tighten thresholds or check data quality |
| Batch effects | Use ComBat or include batch in design formula |
| Unbalanced groups | Consider limma-voom; adjust for covariates |
| No biological replicates | Cannot perform reliable differential analysis |
| Gene ID conversion errors | Use AnnotationDbi or biomaRt for mapping |
| GSEA gives no results | Check gene ranking metric; try different gene sets |

## 10. Experimental Validation Suggestions

1. **qRT-PCR**: Validate top 5-10 DEGs in independent samples
2. **Western blot**: Confirm protein-level changes for key genes
3. **Immunohistochemistry (IHC)**: Validate spatial expression in tissue
4. **Functional assays**: siRNA/shRNA knockdown or overexpression of key genes
5. **Reporter assays**: For pathway activity validation
6. **External dataset validation**: Replicate findings in independent cohorts (GEO)
