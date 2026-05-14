# Multi-Omics Integration Workflow

## 1. Applicable Scenarios

- Transcriptomics + Proteomics joint analysis
- Transcriptomics + Metabolomics joint analysis
- Microbiome + Metabolomics + Host transcriptomics
- Bulk RNA-seq + Single-cell RNA-seq
- Single-cell + Spatial transcriptomics
- Clinical data + Multi-omics
- Cross-omics regulatory network construction

## 2. Input Data Requirements

| Input | Description | Format |
|-------|-------------|--------|
| Omics matrix 1 | Features × Samples (e.g., gene expression) | CSV/TSV |
| Omics matrix 2 | Features × Samples (e.g., protein/metabolite) | CSV/TSV |
| Sample metadata | Shared sample IDs across omics, group labels | CSV/TSV |
| Clinical data (optional) | Survival, clinical variables | CSV/TSV |

**Critical**: Sample IDs must be matchable across all omics layers.

## 3. Data Validation Checklist

Before integration, verify:
```
□ All omics datasets loaded and dimensions confirmed
□ Sample IDs are consistent across omics (same format, no typos)
□ Only shared samples are retained (document how many were dropped)
□ Each omics layer is independently normalized and QC-passed
□ Batch effects checked within each omics layer
□ Outlier samples identified in each layer (PCA-based)
□ Missing values characterized (MCAR/MAR/MNAR) per layer
□ Feature distributions are compatible for cross-omics correlation
□ Sample metadata is complete and aligned
□ Sample sizes per condition are adequate for multi-omics modeling (typically > 20)
```

| Issue | Detection | Action |
|-------|-----------|--------|
| Mismatched sample IDs | `intersect()` across sample lists | Standardize IDs; drop unmatched |
| Different scales across omics | PCA of concatenated data | Normalize per omics layer independently |
| Batch effects across layers | PCA grouped by omics type | ComBat within each layer before integration |
| Small sample overlap | < 15 shared samples | Reduce model complexity; use simpler methods |
| High missingness in one layer | > 30% missing per feature | Remove features; consider imputation |

## 4. Recommended Analysis Steps

### Step 1: Data Harmonization
1. Match samples across omics layers (keep only shared samples)
2. Normalize each omics layer independently
3. Handle missing values within each layer
4. Log-transform if needed
5. Standardize/scale features if cross-omics comparison is needed

### Step 2: Cross-Omics Correlation Analysis
- Pairwise Spearman/Pearson correlation between features across omics
- Filter by correlation threshold (|r| > 0.5, FDR < 0.05)
- Four-quadrant analysis (e.g., mRNA-up/protein-up, mRNA-up/protein-down)

### Step 3: Multi-Omics Factor Analysis

| Method | Package | Features |
|--------|---------|----------|
| MOFA/MOFA+ | MOFA2 (R/Python) | Unsupervised, finds shared variation |
| DIABLO (MINT) | mixOmics | Supervised, multi-block PLS |
| iCluster | iClusterPlus | Clustering across omics |
| SNF | SNFtool | Similarity network fusion |
| WGCNA (per-omics) | WGCNA | Co-expression/co-abundance modules |

### Step 4: Network Integration
1. Build correlation network across omics
2. Community/module detection (Louvain, MCL)
3. Identify hub nodes connecting omics layers
4. Pathway enrichment per network module
5. Visualize integrated network

### Step 5: Pathway Consistency Analysis
- Run enrichment analysis for each omics layer independently
- Compare enriched pathways across layers
- Identify consistently altered pathways (concordant) vs. discordant
- Upset plot or Venn diagram of shared pathways

### Step 6: Integration-Specific Analyses

#### Transcriptomics + Proteomics
- mRNA-protein correlation per gene
- Post-transcriptional regulation identification
- Concordance rate analysis
- Joint pathway analysis

#### Transcriptomics + Metabolomics
- Metabolite-gene correlation network
- Metabolic enzyme gene expression vs. metabolite levels
- Joint KEGG pathway mapping
- Mediation analysis (gene → enzyme → metabolite)

#### Microbiome + Metabolomics + Host
- Microbiome-metabolite Spearman correlation
- Metabolite-host gene correlation
- Tripartite network (microbe → metabolite → host gene)
- Procrustes analysis between ordinations
- Mantel test for distance matrix correlation

#### Bulk + Single-cell
- Deconvolution: estimate cell proportions in bulk using scRNA-seq reference
- CIBERSORTx: impute cell-type-specific gene expression from bulk
- Validate bulk DEGs in specific cell types
- Cell-type-specific differential expression

#### Single-cell + Spatial
- Map scRNA-seq cell types to spatial locations
- Validate cell communication predictions with spatial proximity
- Spatial enrichment of cell type signatures

## 5. Recommended Statistical Methods

| Analysis | Method | Notes |
|----------|--------|-------|
| Cross-omics correlation | Spearman + BH correction | Report |r| > 0.5 with p-adjusted |
| Shared variance | MOFA (MOFA2) | Report variance explained per factor and per omics |
| Supervised integration | DIABLO (mixOmics) | Report classification error rate and selected features per omics |
| Network significance | Permutation test of edge weights | Compare to random network of same size |
| Module-trait association | WGCNA module eigengene correlation | Spearman with BH correction |
| Pathway consistency | Hypergeometric overlap test | Test whether enriched terms from different omics overlap more than expected by chance |
| Integration concordance | Four-quadrant analysis + binomial test | Test whether same-direction changes exceed random expectation |

- MOFA: report the number of factors, variance explained per factor, and top-weighted features
- DIABLO: validate with cross-validation and permutation; report selected features per component
- Four-quadrant analysis: calculate concordance rate and test against null (50%)
- Always compare integrated results to single-omics baselines to demonstrate added value

## 6. Recommended Figures

| Figure | Description |
|--------|-------------|
| Multi-omics overview | Sankey diagram connecting omics layers |
| Correlation heatmap | Cross-omics feature correlations |
| Four-quadrant scatter | mRNA FC vs. Protein FC |
| Network graph | Integrated multi-omics network |
| MOFA factor plot | Shared variation across omics |
| Pathway comparison | Venn/UpSet of shared pathways |
| Circos plot | Cross-omics connections |
| Sankey diagram | Gene → Protein → Metabolite flow |
| Tripartite network | Microbe-metabolite-host |
| Procrustes plot | Ordination alignment between omics |

## 7. Result Interpretation Template

> To investigate the molecular mechanisms underlying [condition], we performed integrative analysis of transcriptomic and proteomic data from [N] matched samples. Cross-omics correlation analysis identified [N] gene-protein pairs with significant concordant changes (Spearman r > 0.5, FDR < 0.05). Overall mRNA-protein concordance rate was XX%. MOFA analysis identified [K] factors capturing shared variation, with Factor 1 (XX% variance) separating [Group1] from [Group2]. Joint pathway analysis revealed [Pathway1] and [Pathway2] as consistently altered at both transcriptomic and proteomic levels. Network integration identified [Gene/Protein] as a hub node connecting transcriptomic and proteomic changes (Figure X).

## 8. Manuscript Writing Template

### Methods

> **Multi-omics data integration**
>
> Matched [transcriptomic / proteomic / metabolomic] data from [N] samples were obtained from [source]. Each omics layer was preprocessed independently: RNA-seq data were normalized using [DESeq2 / TMM], proteomics data were variance-stabilized (vsn), and metabolomics data were Pareto-scaled. Only samples with data from all omics layers were included (n = [N]). Cross-omics pairwise Spearman correlation was performed between features from different omics layers, with BH adjustment (FDR < 0.05, |r| > 0.5 considered significant). Four-quadrant analysis was used to evaluate concordance between protein and mRNA fold changes. Multi-Omics Factor Analysis (MOFA2, vX.X) was applied to identify latent factors capturing shared variation across omics layers. Supervised integration was performed using DIABLO (mixOmics, vX.X) with [condition] as the outcome. Joint pathway enrichment was assessed by testing overlap of significantly enriched KEGG pathways from each omics layer using the hypergeometric test. Integrated multi-omics networks were constructed using correlation-based edge weights, and network modules were detected using the Louvain algorithm (igraph vX.X). Hub nodes were defined as features with degree > 95th percentile and high betweenness centrality. All analyses were performed in R (vX.X) and Python (vX.X).

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| Different sample sizes per omics | Keep only shared samples; imputation with caution |
| Different feature scales | Normalize each omics independently, then scale |
| Low mRNA-protein correlation | Expected (~0.4-0.6); focus on concordant changes |
| Spurious cross-omics correlations | Use FDR correction; validate with pathway consistency |
| Overly complex networks | Filter edges stringently; focus on modules |
| Forced integration | Only integrate if biological question warrants it |

## 10. Experimental Validation Suggestions

1. **Targeted assays**: Validate key cross-omics connections (qRT-PCR + WB + metabolite)
2. **Perturbation experiments**: Knock down hub gene, measure downstream omics changes
3. **Pathway validation**: Inhibitor/activator experiments for key pathways
4. **Independent cohort**: Validate integration findings in new samples
