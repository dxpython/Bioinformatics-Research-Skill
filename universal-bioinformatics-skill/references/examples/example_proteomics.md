# Case study: Quantitative proteomics of drug-resistant cancer cells

## 1. Research questions

**Main question:** What proteome-level changes underlie acquired drug resistance in cancer cells, and which protein networks drive the resistant phenotype?

**Sub-questions:**
1. Which proteins are differentially expressed between drug-sensitive and drug-resistant cells?
2. What protein-protein interaction networks are altered in resistant cells?
3. Which biological pathways are enriched among differentially expressed proteins?
4. Can proteomics identify actionable targets to overcome resistance?
5. How concordant are protein-level and mRNA-level changes (post-transcriptional regulation)?

## 2. Data and design

### Experimental design

| Group | Condition | n (biological replicates) | Samples |
|-------|-----------|--------------------------|---------|
| Sensitive | Parental cell line, vehicle-treated | 4 | Cell lysate |
| Resistant | Drug-selected resistant line | 4 | Cell lysate |
| Sensitive + drug | Parental treated with IC50 for 24h | 4 | Cell lysate |
| Resistant + drug | Resistant treated with IC50 for 24h | 4 | Cell lysate |

### Proteomics platform

| Parameter | Detail |
|-----------|--------|
| Platform | Q Exactive HF-X |
| Method | TMT 16-plex, LC-MS/MS |
| Digestion | Trypsin/LysC |
| Quantification | Reporter ion (MS2) |
| Database search | MaxQuant (v2.4) against UniProt human proteome |
| FDR threshold | 1% at PSM and protein level |

### Optional transcriptomics match

| Layer | Platform | Purpose |
|-------|----------|---------|
| RNA-seq | Illumina NovaSeq, polyA+ | Compare mRNA vs. protein changes |

## 3. Analysis workflow

### Step 1: Data preprocessing and QC

```r
# Load MaxQuant proteinGroups output
.libPaths(c("./Rlib", .libPaths()))
library(tidyverse)
library(vsn)
library(impute)

# Load protein data
protein_groups <- read.delim("proteinGroups.txt", stringsAsFactors = FALSE)
# Keep only identified by site, remove reverse and contaminant hits
proteins <- protein_groups %>%
  filter(Reverse != "+", Potential.contaminant != "+",
         Only.identified.by.site != "+")
# Extract reporter intensity columns
intensity_cols <- grep("Reporter.intensity.corrected", colnames(proteins), value = TRUE)
mat <- as.matrix(proteins[, intensity_cols])
colnames(mat) <- paste0("Sample_", 1:16)  # rename with actual sample labels

# QC: barplot of total intensities per sample
barplot(colSums(mat, na.rm = TRUE), las = 2, main = "Total reporter intensity per sample")
# QC: number of proteins identified per sample
barplot(colSums(!is.na(mat)), las = 2, main = "Proteins identified per sample")
```

### Step 2: Missing value handling

```r
# Assess missingness
missing_rate <- rowMeans(is.na(mat))
hist(missing_rate, main = "Missing value rate per protein")
# Filter proteins with > 50% missing
mat_filt <- mat[missing_rate <= 0.5, ]
# Impute remaining missing values (MinProb method — impute left-censored missing)
mat_imp <- impute.knn(mat_filt)$data  # or use DEP::impute() with MinProb

# If MNAR pattern suspected in certain groups:
# Use mixed imputation: left-censored for MNAR, knn for MAR
```

### Step 3: Normalization

```r
# Variance-stabilizing normalization
library(vsn)
mat_norm <- vsn::justvsn(mat_imp)
# Check normalization effect
meanSdPlot(mat_norm)
# PCA before and after normalization
pca_before <- prcomp(t(mat_imp), scale = TRUE)
pca_after <- prcomp(t(exprs(mat_norm)), scale = TRUE)
```

### Step 4: Differential protein analysis

```r
library(limma)

# Design matrix
group <- factor(c(rep("Sensitive", 4), rep("Resistant", 4),
                   rep("Sensitive_drug", 4), rep("Resistant_drug", 4)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
# Contrast of interest
contrast_matrix <- makeContrasts(
  Resistant_vs_Sensitive = Resistant - Sensitive,
  ResistantDrug_vs_SensitiveDrug = Resistant_drug - Sensitive_drug,
  levels = design
)
fit <- lmFit(exprs(mat_norm), design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2, trend = TRUE)
# Results
res <- topTable(fit2, coef = "Resistant_vs_Sensitive", number = Inf, adjust.method = "BH")
deg_proteins <- res %>%
  filter(adj.P.Val < 0.05, abs(logFC) > log2(1.5))
cat(sprintf("Differentially expressed proteins: %d (up: %d, down: %d)\n",
            nrow(deg_proteins), sum(deg_proteins$logFC > 0), sum(deg_proteins$logFC < 0)))
```

### Step 5: PPI network analysis

```r
library(STRINGdb)

string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
# Map proteins to STRING IDs
res$gene <- sapply(strsplit(res$Gene.names, ";"), `[`, 1)
mapped <- string_db$map(res, "gene", removeUnmappedRows = TRUE)
# Get PPI network
ppi <- string_db$get_interactions(mapped$STRING_id)
# Enrichment
enriched <- string_db$get_enrichment(mapped$STRING_id)
head(enriched[order(enriched$p_value), ], 20)
```

### Step 6: Functional enrichment

```r
library(clusterProfiler)
library(org.Hs.eg.db)

deg_genes <- bitr(deg_proteins$gene, fromType = "SYMBOL", toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
ego <- enrichGO(gene = deg_genes$ENTREZID, OrgDb = org.Hs.eg.db,
                ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05)
ekegg <- enrichKEGG(gene = deg_genes$ENTREZID, organism = "hsa",
                    pAdjustMethod = "BH", qvalueCutoff = 0.05)
dotplot(ego, showCategory = 15) + ggtitle("GO Biological Process")
```

### Step 7: Proteomics-transcriptomics integration (if RNA-seq available)

```r
# Load RNA-seq results, match by gene symbol
rna_res <- read.csv("rna_seq_degs.csv")  # columns: gene, log2FC, padj
merged <- merge(res, rna_res, by.x = "gene", by.y = "gene", suffixes = c("_protein", "_rna"))
# Four-quadrant analysis
ggplot(merged, aes(x = logFC_rna, y = logFC_protein)) +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  labs(x = "mRNA log2FC", y = "Protein log2FC",
       title = "Proteomics-Transcriptomics Concordance") +
  annotate("text", x = Inf, y = Inf, hjust = 1.2, vjust = 2,
           label = sprintf("Spearman rho = %.2f",
                          cor(merged$logFC_protein, merged$logFC_rna, method = "spearman",
                              use = "complete.obs")))
# Identify post-transcriptionally regulated candidates (discordant)
merged$discordant <- with(merged,
  (logFC_protein > 1 & logFC_rna < -1) | (logFC_protein < -1 & logFC_rna > 1))
```

## 4. Figure plan

| Figure | Content | Tool |
|--------|---------|------|
| Fig 1A | Experimental design schematic | BioRender / draw.io |
| Fig 1B | Protein identifications barplot per sample | base R / ggplot2 |
| Fig 1C | PCA score plot (all groups) | ggplot2 |
| Fig 1D | Missing value heatmap | ComplexHeatmap |
| Fig 2A | Volcano plot (Resistant vs. Sensitive) | EnhancedVolcano |
| Fig 2B | Heatmap of top 50 DEPs | ComplexHeatmap |
| Fig 2C | GO enrichment dotplot | clusterProfiler |
| Fig 2D | KEGG pathway enrichment barplot | clusterProfiler |
| Fig 3A | PPI network of DEPs (with log2FC color) | STRING / Cytoscape |
| Fig 3B | Network module enrichment | STRING |
| Fig 3C | Protein-mRNA correlation scatter | ggplot2 |
| Fig 3D | Four-quadrant plot | ggplot2 |
| Fig 4A | Boxplot of top 5 DEPs (validation targets) | ggplot2 |
| Fig 4B | Western blot validation (key targets) | N/A (wet-lab) |
| S1 | Drug treatment comparison volcano | EnhancedVolcano |
| S2 | Full enrichment results table | DT / kableExtra |

## 5. Results text templates

### DEG overview

> Quantitative TMT-based proteomics identified a total of [N] proteins across all samples (FDR < 1%). After filtering for proteins with > 50% valid values, [N_filt] proteins were retained for differential analysis. Comparison between resistant and sensitive cells revealed [N_dep] differentially expressed proteins (DEPs; |log2FC| > 0.585, adjusted p < 0.05), of which [X] were upregulated and [Y] downregulated in resistant cells (Figure 2A-B). The top upregulated proteins included [Protein1] (log2FC = X.XX, padj = X.XXe-XX) and [Protein2], while [Protein3] and [Protein4] were among the most downregulated.

### Pathway enrichment

> GO enrichment analysis of the DEPs revealed enrichment in biological processes related to [Process1] (padj = X.XX), [Process2], and [Process3] (Figure 2C). KEGG pathway analysis further identified [Pathway1] and [Pathway2] as significantly enriched (padj < 0.05), suggesting activation of [biological theme] in resistant cells (Figure 2D).

### PPI and integration

> STRING PPI network analysis revealed a densely connected module containing [N_module] proteins, with [Hub1] and [Hub2] identified as central hub nodes (degree > 95th percentile; Figure 3A). Integration with RNA-seq data showed moderate concordance between protein and mRNA fold changes (Spearman ρ = X.XX). A total of [N_discordant] proteins exhibited discordant regulation (protein change without corresponding mRNA change), indicating post-transcriptional regulation (Figure 3D). Notably, [Protein_X] was upregulated at the protein level despite unchanged mRNA, suggesting increased protein stability in resistant cells.

## 6. Practical notes

- **Missing values**: In TMT proteomics, MNAR (missing not at random) is common for low-abundance proteins. Use mixed imputation approaches.
- **Protein inference**: MaxQuant's proteinGroups may include protein groups. Report the leading razor protein or gene symbol for clarity.
- **FDR vs. p-value**: At the protein level, 1% FDR is standard. Use adjusted p-values throughout.
- **Batch effects**: With TMT plexes, each TMT channel is measured in the same MS run—batch effects are typically small, but still check PCA by TMT batch.
- **Concordance with RNA-seq**: ~40-60% of protein changes are concordant with mRNA (Spearman ρ ~0.4-0.6). Discordant genes are biologically interesting, not errors.
- **Validation**: Western blot of 3-5 key targets is the minimum validation expected. Include loading controls and quantification (densitometry).
- **Drug treatment contrast**: Comparing Resistant+drug vs. Sensitive+drug can reveal proteins specifically targeted or bypassed by the drug in resistant cells.
