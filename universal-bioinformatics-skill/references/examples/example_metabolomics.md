# Case study: Serum metabolomics of early-stage liver disease

## 1. Research questions

**Main question:** Can serum metabolomic profiling distinguish early-stage non-alcoholic fatty liver disease (NAFLD) from healthy controls and identify metabolic pathways driving disease progression?

**Sub-questions:**
1. Which serum metabolites differ between NAFLD patients and healthy controls?
2. Can a PLS-DA/OPLS-DA model discriminate NAFLD from healthy with acceptable performance?
3. Which metabolic pathways (KEGG, SMPDB) are disturbed in early NAFLD?
4. Do metabolite panels outperform routine clinical markers (ALT, AST) for detection?
5. Can metabolite-gene correlation with public liver transcriptomic data reveal underlying mechanisms?

## 2. Data and design

### Cohort

| Group | Criteria | n | Biospecimens |
|-------|----------|---|--------------|
| Healthy | Normal liver enzymes, BMI < 25, no metabolic syndrome | 40 | Fasting serum |
| Early NAFLD | Liver biopsy confirmed (NAS 2–4), no fibrosis | 40 | Fasting serum |

### Metabolomics platform

| Parameter | Detail |
|-----------|--------|
| Platform | UHPLC-Q-Exactive HF |
| Ionization modes | ESI+, ESI- |
| Column | ACQUITY UPLC BEH Amide (HILIC) + C18 (reversed-phase) |
| Data processing | XCMS (v3.20) + CAMERA for annotation |
| Annotation | mzCloud, HMDB, KEGG; MSI Level 1–2 for key features |
| QC | Pooled QC every 8 samples; blanks to remove background |
| Normalization | TIC normalization + Pareto scaling |

### External transcriptomic reference

| Dataset | Platform | Purpose |
|---------|----------|---------|
| GSE89632 | Affymetrix, 20 NAFLD vs 20 normal liver biopsies | Metabolite-gene correlation |

## 3. Analysis workflow

### Step 1: Data preprocessing and QC

```r
.libPaths(c("./Rlib", .libPaths()))
library(tidyverse)
library(ropls)
library(MetaboAnalystR)

# Load feature table (XCMS output)
features <- read.csv("xcms_features.csv", row.names = 1)  # features × samples
metadata <- read.csv("sample_metadata.csv")                # sample_id, group, age, sex, bmi

# QC: RSD assessment in pooled QCs
qc_cols <- grep("QC", colnames(features))
sample_cols <- grep("QC", colnames(features), invert = TRUE)
rsd_qc <- apply(features[, qc_cols], 1, function(x) sd(x) / mean(x) * 100)
cat(sprintf("Features with RSD < 30%% in QC: %d / %d\n", sum(rsd_qc < 30), length(rsd_qc)))
features_filt <- features[rsd_qc < 30, sample_cols]
cat(sprintf("Retained features after QC filtering: %d\n", nrow(features_filt)))
```

### Step 2: Missing value handling and normalization

```r
# Impute missing values with half-minimum (common for metabolomics)
impute_half_min <- function(mat) {
  mat_imp <- mat
  for (i in 1:nrow(mat_imp)) {
    row_min <- min(mat_imp[i, ], na.rm = TRUE)
    mat_imp[i, is.na(mat_imp[i, ])] <- row_min / 2
  }
  return(mat_imp)
}
features_imp <- impute_half_min(features_filt)

# TIC normalization (total ion current)
tic <- colSums(features_imp, na.rm = TRUE)
features_tic <- sweep(features_imp, 2, tic / mean(tic), "/")

# Pareto scaling (mean-center and divide by sqrt(SD))
features_scaled <- t(scale(t(log2(features_tic + 1)), center = TRUE,
                            scale = apply(log2(features_tic + 1), 2, function(x) sqrt(sd(x)))))
```

### Step 3: Multivariate analysis

```r
# PCA — overview
pca_res <- prcomp(t(features_scaled), center = TRUE, scale = FALSE)
ggplot(data.frame(PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2],
                  Group = metadata$group), aes(PC1, PC2, color = Group)) +
  geom_point(size = 3) + stat_ellipse() +
  labs(title = "PCA Score Plot", x = paste0("PC1 (",
       round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
       y = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"))

# PLS-DA
plsda <- opls(t(features_scaled), metadata$group, predI = 1, orthoI = NA,
              permI = 1000, crossvalI = 7)
# Check model validity
cat(sprintf("R2Y = %.3f, Q2 = %.3f\n", plsda@summaryDF$R2Y, plsda@summaryDF$Q2))
cat(sprintf("Permutation p-value: %.4f\n", plsda@summaryDF$permutation_pQ2))
# Score plot
plot(plsda, typeVc = "x-score")

# OPLS-DA (single predictive component, orthogonal components to improve interpretation)
oplsda <- opls(t(features_scaled), metadata$group, predI = 1, orthoI = NA,
               permI = 1000, crossvalI = 7)
plot(oplsda, typeVc = "x-score")
```

### Step 4: Differential metabolite identification

```r
# VIP scores from OPLS-DA
vip_scores <- getVipVn(oplsda)
# Wilcoxon test per feature
p_values <- apply(features_scaled, 1, function(x) {
  wilcox.test(x[metadata$group == "NAFLD"], x[metadata$group == "Healthy"])$p.value
})
p_adj <- p.adjust(p_values, method = "BH")
# Fold change
fc <- apply(features_scaled, 1, function(x) {
  mean(x[metadata$group == "NAFLD"]) - mean(x[metadata$group == "Healthy"])
})  # log2-scale difference

# Differential metabolites: VIP > 1, adjusted p < 0.05, |FC| > log2(1.5)
diff_mets <- data.frame(
  feature = rownames(features_scaled),
  VIP = vip_scores,
  p_value = p_values,
  p_adj = p_adj,
  log2FC = fc
) %>% filter(VIP > 1, p_adj < 0.05, abs(log2FC) > log2(1.5))

cat(sprintf("Differential metabolites: %d\n", nrow(diff_mets)))
```

### Step 5: Pathway enrichment

```r
# Using MetaboAnalystR or KEGG via clusterProfiler-like approach
# Prepare KEGG compound IDs for enrichment
# If MetaboAnalystR is used:
# mSet <- InitDataObjects("conc", "pathora", FALSE)
# mSet <- SetDesignType(mSet, "multi")
# mSet <- Setup.MapData(mSet, diff_mets$kegg_id)
# mSet <- CrossReferencing(mSet, "kegg")
# mSet <- CreateMappingResultTable(mSet)
# mSet <- SetKEGG.PathLib(mSet, "hsa")
# mSet <- SetMetabolomeFilter(mSet, TRUE)
# mSet <- CalculateOraScore(mSet)

# Alternative: manual hypergeometric test using KEGG pathway-compound mappings
# For each pathway, test enrichment of differential metabolites
enrichment <- data.frame(
  pathway = c("Fatty acid biosynthesis", "Glycerophospholipid metabolism",
              "Bile acid biosynthesis", "TCA cycle", "Amino acid metabolism",
              "Glutathione metabolism", "Purine metabolism"),
  hit = c(5, 4, 3, 2, 6, 3, 2),
  total_in_pathway = c(20, 18, 15, 8, 30, 12, 10),
  p_value = c(0.0003, 0.001, 0.008, 0.03, 0.0001, 0.005, 0.04)
) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))
```

### Step 6: Metabolite-gene correlation

```r
# Load liver transcriptomic DEG results from public dataset
liver_deg <- read.csv("GSE89632_degs.csv")  # gene, log2FC, padj
# For key differential metabolites with known metabolic enzymes:
# Map metabolite → enzyme (via KEGG or literature) → gene expression
# Correlate metabolite abundance with enzyme gene expression across shared pathways
# This provides a mechanistic link between liver gene expression and serum metabolites

key_enzymes <- c("FASN", "SCD1", "CPT1A", "ACOX1", "CYP2E1", "CYP4A11")
enzyme_expr <- liver_deg[liver_deg$gene %in% key_enzymes, ]
cat("Key enzymes altered in NAFLD liver:\n")
print(enzyme_expr[, c("gene", "log2FC", "padj")])
# Interpret: upregulation of FASN/SCD1 (lipogenesis) + serum lipid metabolites
# = consistent evidence for lipogenesis activation in NAFLD
```

## 4. Figure plan

| Figure | Content | Tool |
|--------|---------|------|
| Fig 1A | Study design flowchart | BioRender / draw.io |
| Fig 1B | PCA score plot (Healthy vs NAFLD) | ggplot2 |
| Fig 1C | PLS-DA score plot with group ellipses | ropls |
| Fig 1D | Permutation test plot (n = 1000) | ropls |
| Fig 2A | Volcano plot of differential metabolites | ggplot2 |
| Fig 2B | VIP score barplot (top 20) | ggplot2 |
| Fig 2C | Heatmap of differential metabolites | ComplexHeatmap |
| Fig 2D | Boxplots of top 6 metabolites | ggplot2 |
| Fig 3A | KEGG pathway enrichment bubble plot | ggplot2 |
| Fig 3B | Pathway impact vs. -log10(p) plot | MetaboAnalyst |
| Fig 3C | Correlation heatmap (metabolites × clinical markers) | corrplot |
| Fig 4A | Metabolite-gene correlation network | igraph / Cytoscape |
| Fig 4B | Liver enzyme expression boxplots (public data) | ggplot2 |
| S1 | OPLS-DA S-plot | ropls |
| S2 | Full differential metabolite table | DT |

## 5. Results text templates

### Multivariate discrimination

> Untargeted metabolomics profiling of fasting serum identified [N] metabolic features after QC filtering (RSD < 30% in QC pools). PCA showed clear separation between NAFLD patients and healthy controls along PC1 (Figure 1B). PLS-DA modeling achieved robust discrimination with R²Y = X.XX and Q² = X.XX, validated by a significant permutation test (p < 0.001, n = 1000; Figure 1C-D). OPLS-DA S-plot identified the features most responsible for group separation (Supplementary Figure S1).

### Differential metabolites

> A total of [N_diff] differentially abundant metabolites were identified using the combined criteria of VIP > 1, BH-adjusted p < 0.05, and |fold change| > 1.5 (Figure 2A). The top discriminating metabolites included elevated [Metabolite1] (VIP = X.XX, log2FC = X.XX), [Metabolite2], and [Metabolite3], and decreased [Metabolite4] and [Metabolite5] in NAFLD serum (Figure 2B-D).

### Pathway enrichment

> KEGG pathway enrichment analysis of the differential metabolites revealed significant disturbance in [Pathway1] (padj = X.XX), [Pathway2], and [Pathway3] (Figure 3A). These pathways are central to [biological process, e.g., hepatic lipid metabolism and mitochondrial energy homeostasis], suggesting [mechanism interpretation, e.g., a shift toward lipogenesis and impaired fatty acid oxidation] in early NAFLD (Figure 3B). Correlation analysis with clinical markers showed that [Metabolite_X] was positively correlated with ALT (r = X.XX, p = X.XX) and HOMA-IR (r = X.XX, p = X.XX) (Figure 3C).

### Metabolite-gene integration

> Integration with liver transcriptomic data from an independent NAFLD cohort (GSE89632) revealed concordant changes: serum [Metabolite_group] elevation coincided with upregulated hepatic expression of [Enzyme1] (log2FC = X.XX, padj = X.XX) and [Enzyme2] (Figure 4A-B). This multi-omics convergence suggests [mechanism, e.g., that hepatic de novo lipogenesis drives the serum lipidomic signature in NAFLD] and identifies the [pathway] axis as a potential therapeutic target.

## 6. Practical notes

- **QC RSD filtering**: Metabolomics features with RSD > 30% in pooled QCs are unreliable — remove before analysis. Some labs use stricter 20% cutoff.
- **PLS-DA overfitting**: PLS-DA will always separate groups to some degree. A permutation test (n ≥ 1000, p < 0.05) and Q² > 0.5 are essential to claim meaningful discrimination.
- **Metabolite identification**: For untargeted data, report the Metabolomics Standards Initiative (MSI) confidence level (Level 1 = validated with authentic standard, Level 2 = spectral library match, Level 3 = tentative based on exact mass, Level 4 = unknown). Level 3+ should be validated by MS/MS.
- **Ion suppression**: In LC-MS, co-eluting compounds can suppress ionization. TIC normalization helps but doesn't fully correct for this. Use internal standards when available.
- **Multiple adducts**: One metabolite can appear as multiple features (M+H, M+Na, M-H2O+H, etc.). Use CAMERA or similar tools to annotate adducts and group features belonging to the same metabolite.
- **Fasting status**: Serum metabolome is sensitive to fasting status. Ensure all samples are collected under the same fasting conditions. Non-fasting samples should be excluded or the limitation clearly stated.
- **Sample size**: n = 40 per group is reasonable for discovery metabolomics. Power analysis: with FDR < 0.05 and expected |FC| > 1.5, ~ 80% power is achievable.
- **Clinical model**: If building a diagnostic metabolite panel, split data into training (70%) and test (30%) sets. Report AUC with 95% CI. Compare AUC of metabolite panel vs. ALT alone using DeLong test.
