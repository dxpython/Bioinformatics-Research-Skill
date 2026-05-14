# Case study: Transcriptome + metabolome integration in inflammatory metabolic reprogramming

## 1. Research questions

**Main question:** Under inflammatory states (e.g. sepsis, IBD, rheumatoid arthritis), how do transcriptomes and metabolomes co-reprogram to drive disease?

**Sub-questions:**
1. What are the major RNA and metabolite changes between disease and control?
2. Do DEGs and differential metabolites converge on shared pathways?
3. Which gene–metabolite pairs show strong covariation?
4. Can multi-omics networks nominate driver molecules?
5. What metabolic hallmarks emerge (Warburg shift, lipid rewiring, etc.)?

## 2. Data and design

### Cohort

| Group | Description | n | Data |
|-------|-------------|---|------|
| Control | Healthy | 20 | RNA-seq + untargeted metabolomics |
| Disease | Inflammatory phenotype (e.g. sepsis) | 20 | RNA-seq + untargeted metabolomics |

### Pairing

- **Same biospecimen batch:** RNA and metabolomics from the same patient/time point.
- **Unique IDs:** Perfect 1:1 sample mapping across omics.
- **Covariates:** Age, sex, BMI, medications, etc.

## 3. Analysis workflow

### Step 1: Transcriptome (standalone)

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(count_data, col_data, design = ~ condition)
dds <- DESeq(dds)
res_rna <- results(dds, contrast = c("condition", "Disease", "Control"))

deg <- subset(res_rna, abs(log2FoldChange) > 1 & padj < 0.05)
cat("DEG count:", nrow(deg), "\n")

library(clusterProfiler)
ego <- enrichGO(gene = rownames(deg), OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL", ont = "BP")
ekegg <- enrichKEGG(gene = entrez_ids, organism = "hsa")

gene_list <- sort(setNames(res_rna$log2FoldChange, rownames(res_rna)), decreasing = TRUE)
gsea_res <- gseKEGG(geneList = gene_list_entrez, organism = "hsa")
```

### Step 2: Metabolome (standalone)

```r
library(mixOmics)

metab_data[is.na(metab_data)] <- min(metab_data, na.rm = TRUE) / 5

metab_log <- log2(metab_data)
metab_scaled <- scale(metab_log, center = TRUE, scale = apply(metab_log, 2, sd)^0.5)

pca_metab <- pca(metab_scaled, ncomp = 5)
plotIndiv(pca_metab, group = metadata$condition, legend = TRUE)

plsda_metab <- plsda(metab_scaled, metadata$condition, ncomp = 3)
plotIndiv(plsda_metab, group = metadata$condition, ellipse = TRUE, legend = TRUE)

vip_scores <- vip(plsda_metab)
sig_metab <- names(which(vip_scores[, 1] > 1.0))

diff_pvals <- apply(metab_scaled, 2, function(x) {
  wilcox.test(x ~ metadata$condition)$p.value
})
diff_padj <- p.adjust(diff_pvals, method = "BH")
diff_fc <- apply(metab_data, 2, function(x) {
  log2(mean(x[metadata$condition == "Disease"]) /
       mean(x[metadata$condition == "Control"]))
})

sig_metabolites <- intersect(
  names(which(diff_padj < 0.05 & abs(diff_fc) > 0.5)),
  sig_metab
)
```

### Step 3: Cross-pathway enrichment

```r
gene_pathways <- enrichKEGG(gene = deg_entrez, organism = "hsa")

metab_pathways <- enrichKEGG(gene = kegg_compound_ids, organism = "hsa",
                              keyType = "kegg")

common_pathways <- intersect(
  gene_pathways@result$ID[gene_pathways@result$p.adjust < 0.05],
  metab_pathways@result$ID[metab_pathways@result$p.adjust < 0.05]
)

cat("Jointly enriched KEGG pathways:\n")
for (pw in common_pathways) {
  cat(" -", pw, "\n")
}

library(pathview)
for (pw_id in common_pathways) {
  pathview(
    gene.data = gene_fc_vector,
    cpd.data = metab_fc_vector,
    pathway.id = pw_id,
    species = "hsa",
    out.suffix = "multiomics"
  )
}
```

### Step 4: Multi-omics factor models

#### MOFA2

```r
library(MOFA2)

mofa_data <- list(
  RNA = as.matrix(vst_expr[deg_genes, ]),
  Metabolome = as.matrix(metab_scaled[, sig_metabolites])
)

mofa <- create_mofa(mofa_data)

model_opts <- get_default_model_options(mofa)
model_opts$num_factors <- 10

train_opts <- get_default_training_options(mofa)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42

mofa <- prepare_mofa(mofa, model_options = model_opts,
                     training_options = train_opts)

mofa <- run_mofa(mofa)

plot_variance_explained(mofa, plot_total = TRUE)
plot_weights(mofa, view = "RNA", factor = 1, nfeatures = 20)
plot_weights(mofa, view = "Metabolome", factor = 1, nfeatures = 20)

correlate_factors_with_covariates(mofa, covariates = clinical_data)
```

#### mixOmics DIABLO

```r
library(mixOmics)

data_list <- list(
  RNA = t(vst_expr[deg_genes, ]),
  Metabolome = metab_scaled[, sig_metabolites]
)

design <- matrix(0.1, nrow = 2, ncol = 2,
                 dimnames = list(names(data_list), names(data_list)))
diag(design) <- 0

sgccda_res <- block.splsda(
  X = data_list,
  Y = metadata$condition,
  ncomp = 3,
  design = design
)

plotIndiv(sgccda_res, ind.names = FALSE, legend = TRUE)
plotDiablo(sgccda_res, ncomp = 1)

circosPlot(sgccda_res, cutoff = 0.7, line = TRUE)
network(sgccda_res, cutoff = 0.6)
```

#### WGCNA + metabolites

```r
library(WGCNA)

datExpr <- t(vst_expr)
powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
softPower <- sft$powerEstimate

net <- blockwiseModules(datExpr, power = softPower,
                        TOMType = "unsigned",
                        minModuleSize = 30,
                        mergeCutHeight = 0.25)

module_eigengenes <- net$MEs
module_metab_cor <- cor(module_eigengenes, metab_scaled, use = "p")
module_metab_pval <- corPvalueStudent(module_metab_cor, nrow(datExpr))

lactate_cor <- module_metab_cor[, "Lactate"]
best_module <- names(which.max(abs(lactate_cor)))
```

### Step 5: Joint correlation network

```r
library(igraph)

cross_cor <- cor(t(vst_expr[deg_genes, ]), metab_scaled[, sig_metabolites],
                 method = "spearman")
cross_pval <- matrix(NA, nrow = nrow(cross_cor), ncol = ncol(cross_cor))

for (i in 1:nrow(cross_cor)) {
  for (j in 1:ncol(cross_cor)) {
    test <- cor.test(vst_expr[deg_genes[i], ],
                     metab_scaled[, sig_metabolites[j]],
                     method = "spearman")
    cross_pval[i, j] <- test$p.value
  }
}

cross_padj <- matrix(p.adjust(cross_pval, method = "BH"),
                     nrow = nrow(cross_pval))

sig_pairs <- which(abs(cross_cor) > 0.5 & cross_padj < 0.05, arr.ind = TRUE)

edges <- data.frame(
  gene = rownames(cross_cor)[sig_pairs[, 1]],
  metabolite = colnames(cross_cor)[sig_pairs[, 2]],
  correlation = cross_cor[sig_pairs],
  stringsAsFactors = FALSE
)

g <- graph_from_data_frame(edges, directed = FALSE)
V(g)$type <- ifelse(V(g)$name %in% deg_genes, "Gene", "Metabolite")
V(g)$color <- ifelse(V(g)$type == "Gene", "#4DBBD5", "#E64B35")

write.csv(edges, "network_edges.csv", row.names = FALSE)
```

## 4. Network layout (Cytoscape)

```
Core: highest-degree genes/metabolites
  ├── Gene nodes: circles, blue gradient, size ∝ |log2FC|
  ├── Metabolite nodes: diamonds, red gradient, size ∝ |log2FC|
  └── Edges:
       ├── Positive correlation: red solid
       └── Negative correlation: blue dashed
       └── Width ∝ |r|

Outer ring: pathway annotations / background shading per KEGG module

Legend: node types, direction of change, edge weight interpretation
```

### Cytoscape styling checklist

```
1. Import network_edges.csv + node_attributes.csv
2. Styles:
   - Node shape by type (gene = ellipse, metabolite = diamond)
   - Node color by signed log2FC (blue–white–red)
   - Node size by degree (20–80)
   - Edge color by sign of correlation
   - Edge width by |r|
3. Layout: Prefuse Force Directed
4. Annotate pathway clusters as network annotations
```

## 5. Sankey diagram (concept)

```r
library(networkD3)

nodes <- data.frame(
  name = c(
    "Module_turquoise (immune)", "Module_blue (metabolism)",
    "Module_brown (signaling)", "Module_green (oxidative stress)",
    "TCA cycle", "Glycolysis", "Fatty acid metabolism",
    "Amino acid metabolism", "Purine metabolism",
    "Inflammatory signaling",
    "Organic acids", "Amino acids", "Lipids", "Nucleotides", "Carbohydrates"
  )
)

links <- data.frame(
  source = c(0,0,1,1,1,2,3, 4,4,5,5,6,7,7,8,9),
  target = c(4,9,5,6,7,9,5, 10,12,10,14,12,11,13,13,10),
  value  = c(15,8,12,10,7,6,9, 10,5,8,6,12,9,4,7,11)
)

sankeyNetwork(
  Links = links, Nodes = nodes,
  Source = "source", Target = "target", Value = "value",
  NodeID = "name", fontSize = 14, nodeWidth = 30,
  colourScale = JS("d3.scaleOrdinal(d3.schemeCategory20)")
)
```

### How to read the Sankey

```
Left: gene co-expression modules (e.g. turquoise enriched for immune genes)
Middle: KEGG pathways bridging modules (e.g. glycolysis + TCA link immune and metabolic modules → immunometabolism)
Right: metabolite classes with largest flux changes (organic acids reflect bioenergetic reprogramming; amino acids reflect proteolysis)
```

## 6. Figure plan (14 main panels)

| ID | Figure | Type | Purpose |
|----|--------|------|---------|
| Fig 1A | RNA PCA | Scatter | QC |
| Fig 1B | Metabolome PCA | Scatter | QC |
| Fig 2A | RNA volcano | Volcano | DEGs |
| Fig 2B | Metabolite volcano | Volcano | DEMs |
| Fig 3A | GO bubble (genes) | Bubble | Gene pathways |
| Fig 3B | Pathway bubble (metabs) | Bubble | Metabolite pathways |
| Fig 3C | Pathway overlap Venn | Venn | Shared pathways |
| Fig 4 | KEGG pathview | Pathway map | Joint view |
| Fig 5A | MOFA variance | Bar | Factors |
| Fig 5B | MOFA loadings | Loading | Feature weights |
| Fig 6 | DIABLO circos | Circos | Cross-omics links |
| Fig 7 | Gene–metabolite network | Network | Interactions |
| Fig 8 | Sankey | Sankey | Multi-layer mapping |
| Fig 9 | Mechanism schematic | Illustration | Summary model |

## 7. Notes

### Integration pitfalls

1. **Pairing:** RNA and metabolomics must align sample-wise for causal claims.
2. **Scaling:** Normalize each omics layer before joint models.
3. **Feature imbalance:** Thousands of genes vs hundreds of metabolites—tune model weights.
4. **Causality:** Correlation ≠ mechanism; add time series or experiments when possible.
5. **Multiple testing:** Cross-omics correlation screens need strict FDR; expect power loss.

### Inflammatory metabolic signatures (examples)

| Change | Direction | Interpretation |
|--------|-----------|----------------|
| Glycolysis | ↑ | Warburg-like programs in activated leukocytes |
| OXPHOS | ↓ | Mitochondrial stress |
| Lactate | ↑ | Immunosuppressive niche metabolite |
| Tryptophan–kynurenine | ↑ | IDO1 axis, immune escape |
| Arginine depletion | ↓ | T-cell dysfunction |
| Fatty-acid oxidation | Mixed | Macrophage polarization states |
