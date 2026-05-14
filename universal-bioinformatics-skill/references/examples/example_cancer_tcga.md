# Case study: FOXM1 in the tumor immune microenvironment and prognosis

## 1. Research questions

**Main question:** Is FOXM1 (Forkhead Box M1) associated with remodeling of the tumor immune microenvironment (TIME) and with prognosis?

**Sub-questions:**
1. Is FOXM1 differentially expressed between tumor and normal tissue?
2. Is high FOXM1 expression linked to worse prognosis?
3. Which immune infiltration patterns correlate with FOXM1?
4. Which pathways are enriched in the FOXM1-high group?
5. Can we build a multi-gene prognostic signature that includes FOXM1?

## 2. Data sources

| Data | Source | Notes |
|------|--------|-------|
| RNA-seq matrix | TCGA-LUAD | ~500+ LUAD samples |
| Clinical | TCGA Clinical | OS, TNM, age, sex, etc. |
| Normal reference | GTEx + TCGA normal | Paired/near-normal controls |
| External validation | GEO GSE72094 | Independent LUAD cohort |

### Download / prep (R sketch)

```r
# Option A: TCGAbiolinks
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
data <- GDCprepare(query)

# Option B: UCSC Xena preprocessed matrices
# https://xenabrowser.net/datapages/
# Download HTSeq counts + phenotype for TCGA-LUAD
```

## 3. Analysis workflow

### Step 1: Tumor vs normal differential expression

**Goal:** Confirm FOXM1 upregulation in tumor and obtain global DEGs.

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

keep <- rowSums(counts(dds) >= 10) >= 10
dds <- dds[keep, ]

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"))

foxm1_res <- res["FOXM1", ]
cat("FOXM1 log2FC:", foxm1_res$log2FoldChange,
    "padj:", foxm1_res$padj, "\n")
```

**Output:** DEG table (|log2FC| > 1, padj < 0.05).

### Step 2: FOXM1 across clinical subgroups

**Goal:** Describe FOXM1 by stage, TNM, etc.

```r
library(ggplot2)
library(ggpubr)

# Compare FOXM1: tumor vs normal; by stage I–IV; by T/N/M

median_exp <- median(foxm1_expression)
clinical$FOXM1_group <- ifelse(
  foxm1_expression > median_exp, "High", "Low"
)
```

### Step 3: Survival analysis

**Goal:** Relate FOXM1 to OS/PFS.

```r
library(survival)
library(survminer)

fit <- survfit(Surv(OS.time, OS.status) ~ FOXM1_group, data = clinical)

ggsurvplot(fit,
  data = clinical,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#E7B800", "#2E9FDF"),
  xlab = "Time (days)",
  ylab = "Overall Survival Probability"
)

cox_uni <- coxph(Surv(OS.time, OS.status) ~ FOXM1_exp, data = clinical)

cox_multi <- coxph(
  Surv(OS.time, OS.status) ~ FOXM1_exp + age + gender + stage,
  data = clinical
)
```

### Step 4: Immune infiltration

**Goal:** Link FOXM1 to immune context.

```r
# CIBERSORT: TPM matrix in → 22-cell proportions

library(GSVA)
immune_signatures <- list(
  CD8_T = c("CD8A", "CD8B", "GZMB", "PRF1"),
  Macrophage_M1 = c("NOS2", "IL1B", "TNF", "IL6"),
  Macrophage_M2 = c("CD163", "MRC1", "MSR1", "ARG1"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "IKZF2"),
  NK = c("NCR1", "KLRD1", "NKG7", "GNLY")
)

ssgsea_scores <- gsva(
  expr_matrix, immune_signatures,
  method = "ssgsea", kcdf = "Gaussian"
)

library(estimate)
# ImmuneScore, StromalScore, ESTIMATEScore

# Cross-check with TIMER2.0 web results when applicable
```

### Step 5: GSEA

**Goal:** Pathways enriched in FOXM1-high biology.

```r
library(clusterProfiler)
library(org.Hs.eg.db)

gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)

gsea_go <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  keyType = "SYMBOL",
  pvalueCutoff = 0.05
)

gsea_kegg <- gseKEGG(
  geneList = gene_list_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)

library(msigdbr)
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
gsea_hallmark <- GSEA(gene_list, TERM2GENE = hallmark[, c("gs_name", "gene_symbol")])
```

### Step 6: Co-expression and PPI

**Goal:** FOXM1-correlated genes and STRING/Cytoscape network.

```r
cor_genes <- apply(expr_matrix, 1, function(x) {
  cor.test(x, foxm1_expression, method = "spearman")
})

# Top 50 positive correlates → STRING → Cytoscape
```

### Step 7: ML prognostic signature

**Goal:** LASSO-Cox signature including FOXM1.

```r
library(glmnet)

x <- as.matrix(expr_selected_genes)
y <- Surv(clinical$OS.time, clinical$OS.status)

cv_fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
coef_lasso <- coef(cv_fit, s = "lambda.min")

selected_genes <- rownames(coef_lasso)[which(coef_lasso != 0)]
risk_score <- as.numeric(x[, selected_genes] %*% coef_lasso[selected_genes, ])

library(timeROC)
roc <- timeROC(
  T = clinical$OS.time,
  delta = clinical$OS.status,
  marker = risk_score,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE
)
```

## 4. Figure plan (12 main panels)

| ID | Figure | Type | Role | Tool |
|----|--------|------|------|------|
| Fig 1A | FOXM1 pan-cancer | Box | Pan-tumor view | ggplot2 |
| Fig 1B | Tumor vs normal FOXM1 | Paired box | Validation | ggpubr |
| Fig 1C | FOXM1 vs stage | Box | Clinical association | ggpubr |
| Fig 2A | Volcano | Volcano | Global DE | EnhancedVolcano |
| Fig 2B | Top 50 DEG heatmap | Heatmap | Patterns | ComplexHeatmap |
| Fig 3A | KM curves | Survival | Prognosis | survminer |
| Fig 3B | Uni/multivariable Cox forest | Forest | Independence | forestplot |
| Fig 4A | Immune scores vs FOXM1 | Lollipop / heatmap | TIME | ggplot2 |
| Fig 4B | Immune cells high vs low FOXM1 | Box | Contrast | ggpubr |
| Fig 5A | GSEA ridge | Ridge | Pathways | enrichplot |
| Fig 5B | Hallmark bubble | Bubble | Overview | ggplot2 |
| Fig 6 | LASSO path + ROC | Composite | Model | glmnet / timeROC |

### Supplementary

| ID | Content |
|----|---------|
| Fig S1 | PCA QC |
| Fig S2 | CIBERSORT stacked bar |
| Fig S3 | PPI network export |
| Fig S4 | GEO validation KM |
| Fig S5 | Nomogram |

## 5. Example Results paragraphs (English)

### 3.1 FOXM1 is elevated in LUAD

> FOXM1 expression was significantly upregulated in LUAD tumor tissues compared with adjacent normal tissues (log2FC = 2.34, adjusted P < 0.001, Fig. 1A-B). Pan-cancer analysis across 33 TCGA cancer types revealed that FOXM1 was overexpressed in 25 tumor types (Fig. 1A). Furthermore, FOXM1 expression was positively correlated with advanced clinical stage (P for trend < 0.01, Fig. 1C), suggesting its potential role in tumor progression.

### 3.2 High FOXM1 and poor prognosis

> Kaplan-Meier survival analysis demonstrated that patients in the FOXM1-high group had significantly shorter overall survival (OS) than those in the FOXM1-low group (median OS: 38.2 vs. 62.7 months, log-rank P = 0.003, Fig. 3A). Univariate Cox regression identified FOXM1 as a significant prognostic factor (HR = 1.85, 95% CI: 1.32-2.59, P < 0.001). After adjusting for age, gender, and TNM stage, multivariate Cox regression confirmed FOXM1 as an independent prognostic biomarker (HR = 1.62, 95% CI: 1.14-2.31, P = 0.008, Fig. 3B).

### 3.3 FOXM1 and immune remodeling

> To explore the relationship between FOXM1 and the tumor immune microenvironment (TIME), we performed CIBERSORT and ssGSEA analyses. FOXM1-high tumors exhibited significantly higher infiltration of M2 macrophages (P < 0.001) and regulatory T cells (Tregs, P = 0.005), while showing reduced CD8+ T cell infiltration (P = 0.012, Fig. 4A-B). These findings suggest that FOXM1 overexpression may contribute to an immunosuppressive microenvironment.

### 3.4 Pathway signatures in FOXM1-high tumors

> GSEA revealed that the FOXM1-high group was significantly enriched in cell cycle-related pathways (NES = 2.45, FDR < 0.001), DNA repair (NES = 2.12, FDR < 0.001), and MYC targets (NES = 1.98, FDR = 0.002, Fig. 5A-B). Conversely, immune-related pathways such as interferon-gamma response and allograft rejection were downregulated in the FOXM1-high group, consistent with the immunosuppressive phenotype observed in immune infiltration analysis.

### 3.5 Multi-gene prognostic model

> Using LASSO-Cox regression, we identified a 7-gene signature including FOXM1, AURKA, CCNB1, CDK1, TOP2A, BIRC5, and MKI67 (Fig. 6). The risk score effectively stratified patients into high-risk and low-risk groups with distinct survival outcomes (P < 0.001). Time-dependent ROC analysis showed AUC values of 0.78, 0.75, and 0.72 for 1-year, 3-year, and 5-year OS prediction, respectively. The model was further validated in an independent GEO cohort (GSE72094), yielding comparable performance (3-year AUC = 0.71, Fig. S4).

## 6. Pitfalls and reviewer-facing notes

1. **Batch effects:** When merging TCGA + GTEx, apply batch correction (e.g. ComBat-seq).
2. **Survival cutpoints:** Besides the median, consider `surv_cutpoint()` for optimal splits (validate externally).
3. **Immune deconvolution:** Use at least two methods (e.g. CIBERSORT + ssGSEA) for consistency.
4. **Multiple testing:** Apply FDR across all high-throughput tests.
5. **External validation:** Always test signatures in independent cohorts to limit overfitting.

### Common reviewer requests

| Request | Mitigation |
|---------|------------|
| External validation | Hold out GEO cohort(s) |
| Rationale for FOXM1 | Literature + unbiased screening |
| Immune validation | Multi-method deconvolution + IHC where possible |
| Overfitting concerns | 10-fold CV + external validation + calibration |
| Mechanism depth | Perturbation experiments or focused literature synthesis |
