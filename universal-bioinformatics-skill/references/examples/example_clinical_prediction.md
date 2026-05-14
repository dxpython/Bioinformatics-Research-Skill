# Case study: Multi-gene prognostic model for tumor risk stratification

## 1. Research questions

**Main question:** Can we build a robust multi-gene prognostic model from transcriptomics to stratify tumor patients?

**Sub-questions:**
1. Which genes associate independently with prognosis?
2. How do we select an optimal feature subset from candidates?
3. What is in-sample performance?
4. Does the model replicate in external cohorts?
5. How do we package results into a clinically interpretable tool (nomogram)?

## 2. Data plan

| Dataset | Source | Role | Approx. n | Notes |
|---------|--------|------|-----------|-------|
| Training | TCGA-LUAD | Model fitting | ~500 | Primary cohort |
| Internal test | Random 50% split of TCGA-LUAD | Internal hold-out | ~250 | Same platform |
| External 1 | GEO GSE72094 | Validation | ~400 | Independent |
| External 2 | GEO GSE31210 | Validation | ~200 | Independent |

### Inclusion / exclusion

**Include:**
- Available RNA-seq or microarray expression
- Complete survival (time + event)
- Follow-up > 30 days

**Exclude:**
- Missing key clinical variables (stage, age)
- Non-positive survival times
- Neoadjuvant therapy (if the analysis requires treatment-naïve tumors)

## 3. Analysis workflow

### Step 1: Prepare expression and survival tables

```r
library(TCGAbiolinks)
library(GEOquery)
library(survival)
library(survminer)

query <- GDCquery(
  project = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query)
tcga_data <- GDCprepare(query)

expr_tcga <- assay(tcga_data, "tpm_unstrand")
clinical_tcga <- as.data.frame(colData(tcga_data))

clinical_tcga <- clinical_tcga %>%
  filter(!is.na(days_to_death) | !is.na(days_to_last_follow_up)) %>%
  mutate(
    OS.time = ifelse(!is.na(days_to_death), days_to_death, days_to_last_follow_up),
    OS.status = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(OS.time > 30)

gse72094 <- getGEO("GSE72094", GSEMatrix = TRUE)
expr_geo1 <- exprs(gse72094[[1]])
clinical_geo1 <- pData(gse72094[[1]])

common_genes <- intersect(rownames(expr_tcga), rownames(expr_geo1))
```

### Step 2: Candidate gene screening

```r
library(survival)

uni_cox_results <- data.frame()
for (gene in rownames(expr_tcga)) {
  tryCatch({
    cox_fit <- coxph(Surv(OS.time, OS.status) ~ expr_tcga[gene, ],
                     data = clinical_tcga)
    s <- summary(cox_fit)
    uni_cox_results <- rbind(uni_cox_results, data.frame(
      gene = gene,
      HR = s$coefficients[, "exp(coef)"],
      pvalue = s$coefficients[, "Pr(>|z|)"],
      CI_lower = s$conf.int[, "lower .95"],
      CI_upper = s$conf.int[, "upper .95"]
    ))
  }, error = function(e) NULL)
}

candidate_genes <- uni_cox_results$gene[uni_cox_results$pvalue < 0.05]
cat("Univariate Cox candidates:", length(candidate_genes), "\n")

candidate_genes <- intersect(candidate_genes, deg_genes)
```

### Step 3: LASSO-Cox

```r
library(glmnet)

x_train <- t(expr_tcga[candidate_genes, ])
y_train <- Surv(clinical_tcga$OS.time, clinical_tcga$OS.status)

set.seed(42)
cv_fit <- cv.glmnet(
  x = x_train,
  y = y_train,
  family = "cox",
  alpha = 1,
  nfolds = 10,
  type.measure = "C"
)

cat("lambda.min:", cv_fit$lambda.min, "\n")
cat("lambda.1se:", cv_fit$lambda.1se, "\n")

coef_optimal <- coef(cv_fit, s = "lambda.min")
selected_genes <- rownames(coef_optimal)[which(coef_optimal != 0)]
selected_coefs <- coef_optimal[selected_genes, ]

cat("Genes after LASSO:", length(selected_genes), "\n")
cat(paste(selected_genes, collapse = ", "), "\n")

plot(cv_fit)
title("LASSO-Cox cross-validation")

plot(cv_fit$glmnet.fit, xvar = "lambda")
abline(v = log(cv_fit$lambda.min), lty = 2, col = "red")
abline(v = log(cv_fit$lambda.1se), lty = 2, col = "blue")
```

### Step 4: Risk score

```r
calculate_risk_score <- function(expr_matrix, genes, coefficients) {
  expr_subset <- expr_matrix[genes, ]
  risk_score <- as.numeric(t(coefficients) %*% expr_subset)
  return(risk_score)
}

risk_train <- calculate_risk_score(expr_tcga, selected_genes, selected_coefs)

median_risk <- median(risk_train)
clinical_tcga$risk_score <- risk_train
clinical_tcga$risk_group <- ifelse(risk_train > median_risk, "High", "Low")

risk_geo1 <- calculate_risk_score(expr_geo1, selected_genes, selected_coefs)
clinical_geo1$risk_score <- risk_geo1
clinical_geo1$risk_group <- ifelse(risk_geo1 > median_risk, "High", "Low")
```

### Step 5: Kaplan–Meier

```r
fit_train <- survfit(Surv(OS.time, OS.status) ~ risk_group, data = clinical_tcga)
p_train <- ggsurvplot(fit_train,
  data = clinical_tcga,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#E64B35", "#4DBBD5"),
  title = "Training (TCGA-LUAD)",
  xlab = "Time (days)",
  ylab = "Overall survival"
)

fit_geo1 <- survfit(Surv(OS.time, OS.status) ~ risk_group, data = clinical_geo1)
p_geo1 <- ggsurvplot(fit_geo1,
  data = clinical_geo1,
  pval = TRUE,
  risk.table = TRUE,
  palette = c("#E64B35", "#4DBBD5"),
  title = "Validation (GSE72094)",
  xlab = "Time (days)",
  ylab = "Overall survival"
)
```

### Step 6: Time-dependent ROC

```r
library(timeROC)

roc_train <- timeROC(
  T = clinical_tcga$OS.time,
  delta = clinical_tcga$OS.status,
  marker = clinical_tcga$risk_score,
  cause = 1,
  weighting = "marginal",
  times = c(365, 1095, 1825),
  iid = TRUE
)

cat("Training AUC (1/3/5y):\n")
cat("  1y:", round(roc_train$AUC[1], 3), "\n")
cat("  3y:", round(roc_train$AUC[2], 3), "\n")
cat("  5y:", round(roc_train$AUC[3], 3), "\n")

plot(roc_train, time = 365, col = "#E64B35", title = "Time-dependent ROC")
plot(roc_train, time = 1095, col = "#4DBBD5", add = TRUE)
plot(roc_train, time = 1825, col = "#00A087", add = TRUE)
legend("bottomright",
       legend = c(paste0("1y AUC = ", round(roc_train$AUC[1], 3)),
                  paste0("3y AUC = ", round(roc_train$AUC[2], 3)),
                  paste0("5y AUC = ", round(roc_train$AUC[3], 3))),
       col = c("#E64B35", "#4DBBD5", "#00A087"),
       lwd = 2)

roc_geo1 <- timeROC(
  T = clinical_geo1$OS.time,
  delta = clinical_geo1$OS.status,
  marker = clinical_geo1$risk_score,
  cause = 1,
  times = c(365, 1095, 1825),
  iid = TRUE
)
```

### Step 7: Multivariable Cox + forest plot

```r
multi_cox <- coxph(
  Surv(OS.time, OS.status) ~ risk_score + age + gender + stage,
  data = clinical_tcga
)

summary(multi_cox)

library(forestplot)
cox_summary <- summary(multi_cox)
forest_data <- data.frame(
  Variable = c("Risk score", "Age", "Gender (male)", "Stage II", "Stage III", "Stage IV"),
  HR = cox_summary$conf.int[, 1],
  Lower = cox_summary$conf.int[, 3],
  Upper = cox_summary$conf.int[, 4],
  P = cox_summary$coefficients[, 5]
)
```

### Step 8: Nomogram

```r
library(rms)

ddist <- datadist(clinical_tcga)
options(datadist = "ddist")

cox_nom <- cph(
  Surv(OS.time / 365, OS.status) ~ risk_score + age + gender + stage,
  data = clinical_tcga,
  x = TRUE, y = TRUE, surv = TRUE, time.inc = 3
)

nom <- nomogram(cox_nom,
  fun = list(
    function(x) surv(365/365, x),
    function(x) surv(1095/365, x),
    function(x) surv(1825/365, x)
  ),
  fun.at = c(0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
  funlabel = c("1-year survival", "3-year survival", "5-year survival"),
  lp = FALSE
)

plot(nom, xfrac = 0.3)
```

### Step 9: Calibration

```r
cal <- calibrate(cox_nom, cmethod = "KM", method = "boot",
                 u = 3, m = 50, B = 1000)

plot(cal,
  xlab = "Nomogram-predicted 3-year OS probability",
  ylab = "Observed 3-year OS (Kaplan-Meier)",
  main = "Calibration (3-year OS)"
)

for (time_point in c(1, 3, 5)) {
  cal_t <- calibrate(cox_nom, cmethod = "KM", method = "boot",
                     u = time_point, m = 50, B = 1000)
  plot(cal_t, main = paste0("Calibration (", time_point, "-year OS)"))
}
```

### Step 10: Decision curve analysis (DCA)

```r
library(ggDCA)

model_risk <- coxph(Surv(OS.time, OS.status) ~ risk_score, data = clinical_tcga)
model_stage <- coxph(Surv(OS.time, OS.status) ~ stage, data = clinical_tcga)
model_full <- coxph(Surv(OS.time, OS.status) ~ risk_score + age + gender + stage,
                    data = clinical_tcga)

dca_res <- dca(model_risk, model_stage, model_full,
               times = c(365, 1095, 1825))

ggplot(dca_res)
```

### Step 11: External validation metrics

```r
library(Hmisc)

c_train <- concordance.index(
  x = clinical_tcga$risk_score,
  surv.time = clinical_tcga$OS.time,
  surv.event = clinical_tcga$OS.status
)
cat("Training C-index:", round(c_train$c.index, 3),
    "(95% CI:", round(c_train$lower, 3), "-", round(c_train$upper, 3), ")\n")

c_geo1 <- concordance.index(
  x = clinical_geo1$risk_score,
  surv.time = clinical_geo1$OS.time,
  surv.event = clinical_geo1$OS.status
)
cat("Validation C-index:", round(c_geo1$c.index, 3),
    "(95% CI:", round(c_geo1$lower, 3), "-", round(c_geo1$upper, 3), ")\n")
```

## 4. Figure plan (14 main panels)

| ID | Figure | Type | Purpose |
|----|--------|------|---------|
| Fig 1A | Univariate Cox forest | Forest | Candidate genes |
| Fig 1B | LASSO coefficient path | Path | Shrinkage |
| Fig 1C | CV error vs lambda | CV curve | Lambda choice |
| Fig 2A | KM training | Survival | Stratification |
| Fig 2B | KM external 1 | Survival | Replication |
| Fig 2C | KM external 2 | Survival | Second cohort |
| Fig 3A | tdROC training | ROC | Discrimination |
| Fig 3B | tdROC external | ROC | External AUC |
| Fig 4 | Multivariable forest | Forest | Independence |
| Fig 5 | Nomogram | Nomogram | Clinical tool |
| Fig 6A | Calibration training | Calibration | Agreement |
| Fig 6B | Calibration external | Calibration | External agreement |
| Fig 7 | DCA | DCA | Net benefit |
| Fig 8 | Risk heatmap trio | Heatmap + tracks | Signature display |

### Risk “trio” plot layout

```
Top: risk score ranked samples (color = high/low)
Middle: survival time vs event/censor markers
Bottom: heatmap of signature genes (z-scored expression)
```

## 5. Validation hierarchy

```
Internal:
├── Random 70/30 split
├── 10-fold CV
└── Bootstrap optimism (B = 1000)

External (same assay family):
├── GEO same cancer type
├── Freeze training coefficients + cutoff
└── Recompute KM / ROC / calibration / DCA

External (cross-platform):
├── Microarray vs RNA-seq harmonization
├── Quantile or z-score alignment
└── Stress-test generalization

Prospective (gold standard):
├── Prospective enrollment + predefined lab panel
└── Compare predicted vs observed management/outcomes
```

### Cross-platform normalization sketch

```r
library(preprocessCore)

expr_geo_norm <- normalize.quantiles.use.target(
  x = as.matrix(expr_geo1),
  target = apply(expr_tcga, 1, mean)
)

zscore_normalize <- function(expr_matrix) {
  t(scale(t(expr_matrix)))
}

expr_tcga_z <- zscore_normalize(expr_tcga)
expr_geo1_z <- zscore_normalize(expr_geo1)
```

## 6. Modeling pitfalls and reviewer FAQ

### Pitfalls

1. **Leakage:** Feature selection must occur inside training folds only.
2. **Overfitting:** Keep signature size modest (often 5–20 genes) and validate externally.
3. **Cutpoints:** Median split is common; `surv_cutpoint` or X-tile require validation.
4. **Endpoints:** Match OS vs PFS vs DFS between training and validation.
5. **Platforms:** Harmonize RNA-seq vs microarray scales before scoring.

### Reviewer table

| Comment | Mitigation |
|---------|------------|
| Too many/few genes | Sensitivity analysis across lambda grid |
| Missing external data | At least one independent GEO cohort |
| Compare to existing models | C-index / AUC / NRI / IDI head-to-head |
| Clinical utility | Nomogram + DCA for net benefit |
| Interpretability | Biological annotation per gene |
| Subgroup performance | Stage / age / sex stratified validation |

### C-index benchmarks

| Range | Interpretation |
|-------|----------------|
| 0.5 | Chance |
| 0.5–0.6 | Poor |
| 0.6–0.7 | Fair |
| 0.7–0.8 | Good |
| 0.8–0.9 | Excellent |
| > 0.9 | Suspiciously high—check overfitting |
