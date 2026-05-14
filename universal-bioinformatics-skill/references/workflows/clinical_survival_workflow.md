# Clinical Survival Analysis Workflow

## 1. Applicable Scenarios

1. **Single-gene / multi-gene prognosis**: Relate expression or mutation status to OS/DFS/PFS.
2. **Risk score models**: Multi-gene prognostic signatures.
3. **Nomograms**: Integrate clinical and molecular predictors for individualized prognosis.
4. **Independent prognostic factors**: Univariate and multivariate Cox regression.
5. **Model evaluation**: timeROC, calibration curves, DCA (Decision Curve Analysis).
6. **Training–validation design**: Internal + external validation.
7. **Immunotherapy / chemotherapy response**: Link molecular markers to treatment outcome.

## 2. Input Data Requirements

| Data type | Format | Description |
|-----------|--------|-------------|
| Clinical table | CSV / TSV | sample_id, OS_time, OS_status, DFS_time, DFS_status, age, sex, stage, grade, etc. |
| Expression matrix | CSV / TSV | Rows = genes, columns = samples; TPM / FPKM / normalized counts |
| Molecular features | CSV / TSV | TMB, MSI, mutation status, CNV, methylation scores, etc. |
| Immune scores (optional) | CSV / TSV | CIBERSORT, ssGSEA, TIMER, etc. |

**Key requirements**:
- Consistent time units (days/months); months common for reporting.
- Event coding: 0 = censored (alive/lost to follow-up), 1 = event (death/recurrence).
- At least OS; add DFS/PFS/RFS per study aims.
- Training set: suggest ≥100 samples and ≥20 events.
- External validation should be independent (ideally different center/platform).
- For TCGA, prefer harmonized clinical data from UCSC Xena or GDC.

## 3. Data Validation Checklist

### 3.1 Clinical completeness
```
Checks:
1. Survival missingness:
   - Exclude samples missing OS_time or OS_status
   - Exclude or fix OS_time <= 0
   - Recode OS_status if not 0/1

2. Follow-up:
   - Median follow-up (reverse KM)
   - Distribution: sufficient follow-up to observe events
   - Event rate = events / N; suggest > 20%

3. Covariates:
   - Missingness per variable
   - Categorical level counts (avoid very small groups)
   - Continuous distributions (transform or bin if needed)
```

### 3.2 Align expression with clinical
```
1. Match sample IDs between expression and clinical tables
2. Duplicate biospecimens: keep most recent surgery or average per patient
3. Confirm numeric expression; avoid NA in key genes
```

### 3.3 Train/test split checks
```
1. Baseline Table 1: training vs validation
2. Key variables (age, stage, sex) should be balanced (or use external cohorts)
3. Prefer independent external cohorts when possible
```

## 4. Recommended Analysis Steps

### Step 1: Kaplan–Meier survival
```
1.1 Single-gene KM
    - Stratify by high/low expression (median / optimal cutoff / tertiles)
    - Plot KM curves
    - Log-rank test
    - Report median OS and HR (95% CI)
    - Optimal cutoff: surv_cutpoint() (survminer; maximally selected rank statistics)

1.2 Multi-gene / risk score KM
    - Split by median or optimal risk score
    - KM + log-rank p-value
    - Number at risk table

1.3 Subgroup analysis
    - Stratify by stage (e.g. I–II vs III–IV), age (≤60 vs >60), etc.
    - KM within each subgroup
    - Consistency of risk score across subgroups
```

### Step 2: Cox regression
```
2.1 Univariate Cox
    - Test each clinical/molecular feature vs outcome
    - Variables with p < 0.05 or p < 0.1 for multivariable model
    - Report HR, 95% CI, p-value
    - Forest plot for all univariate results

2.2 Multivariate Cox
    - Include univariate-significant variables
    - Independent prognostic value
    - Risk score significance after adjusting for clinical covariates
    - Check multicollinearity (VIF < 5)

2.3 Proportional hazards (PH) checks
    - cox.zph() and Schoenfeld residual plots
    - If PH violated: stratified Cox or time-dependent covariates
```

### Step 3: Risk score construction
```
3.1 Gene selection
    - Candidates: univariate Cox p < 0.05
    - Or combine DEG screening + Cox

3.2 Model fitting
    Method A: LASSO Cox (glmnet, family = "cox")
    - 10-fold CV for lambda.min or lambda.1se
    - Extract nonzero genes and coefficients

    Method B: Stepwise Cox (forward/backward/both)
    - AIC-based model selection

    Method C: Random survival forest
    - Variable importance; top N for Cox model

3.3 Risk score
    - Risk Score = Σ (coefficient_i × expression_i)
    - Median split high/low risk
    - Apply same formula/coefficients in training and validation
```

### Step 4: Time-dependent ROC (timeROC)
```
4.1 Time-dependent AUC
    - timeROC: 1-, 3-, 5-year AUC (survivalROC alternative)
    - Rule of thumb: AUC > 0.7 acceptable, > 0.8 good

4.2 ROC curves
    - Overlay multiple time points
    - Compare training vs validation

4.3 vs clinical-only models
    - AUC of risk score vs stage/age alone
    - Joint model (risk score + clinical) vs single predictors
```

### Step 5: Nomogram
```
5.1 Build nomogram
    - rms: cph() → nomogram()
    - Include independent factors from multivariable Cox
    - Predict 1-, 3-, 5-year survival probabilities

5.2 Calibration
    - Predicted vs observed survival
    - calibrate(), bootstrap = 1000
    - Ideal: curve along 45° line

5.3 C-index
    - Discrimination
    - C-index > 0.7 often acceptable
    - rcorr.cens() or concordance()
```

### Step 6: DCA (Decision Curve Analysis)
```
6.1 Concept
    - Net benefit across threshold probabilities
    - Compare to "treat all" and "treat none"

6.2 Plotting
    - stdca() / rmda / ggDCA
    - Compare risk score vs clinical vs combined models

6.3 Clinical impact curves (optional)
    - Number treated vs true high-risk at each threshold
```

### Step 7: Training / validation framework
```
7.1 Internal validation
    - Random split 7:3 or 6:4
    - Or K-fold CV
    - Bootstrap (rms::validate) recommended

7.2 External validation
    - Independent GEO / ICGC / in-house cohorts
    - Apply training coefficients to validation
    - Compare KM, Cox, timeROC, calibration

7.3 Meta-validation (optional)
    - Pool HR/AUC across multiple external datasets
```

## 5. Recommended Statistical Methods

| Goal | Method | R package |
|------|---------|-----------|
| Survival curves | Kaplan–Meier | survival, survminer |
| Group comparison | Log-rank | survival |
| Univariate prognosis | Univariate Cox | survival |
| Multivariable prognosis | Multivariate Cox | survival |
| Risk score | LASSO Cox | glmnet |
| Time-dependent ROC | Time-dependent AUC | timeROC, survivalROC |
| Nomogram | Cox → nomogram | rms |
| Calibration | Bootstrap calibration | rms |
| C-index | Harrell's C | survival, rms |
| DCA | Net benefit | ggDCA, rmda |
| Optimal cutoff | Maximally selected rank | survminer, maxstat |
| PH test | Schoenfeld residuals | survival (cox.zph) |
| Competing risks | Fine–Gray | cmprsk, tidycmprsk |

## 6. Recommended Figures

| Figure | Purpose | Tool |
|--------|---------|------|
| **KM curves** | Compare survival | survminer::ggsurvplot |
| **Forest plot** | Uni/multivariable Cox | forestplot, ggforest |
| **Risk score distribution** | Risk score + status + expression | ggplot2 |
| **timeROC** | 1/3/5-year discrimination | timeROC, ggplot2 |
| **Nomogram** | Individualized prediction | rms::nomogram |
| **Calibration** | Calibration performance | rms::calibrate |
| **DCA** | Clinical utility | ggDCA |
| **Heatmap** | Signature genes by risk group | ComplexHeatmap |
| **Subgroup forest** | Subgroup HR summary | forestplot |
| **C-index barplot** | Model comparison | ggplot2 |
| **Sankey** | Risk group → outcomes | ggalluvial |

---


## 7. Result Interpretation Template

### KM
```
Kaplan–Meier analysis showed significantly worse OS in the high [gene/risk score] group vs low
(median OS: [X] vs [Y] months, HR = [Z], 95% CI: [a]–[b], log-rank p = [P]).
5-year survival rates were [X%] and [Y%].
```

### Cox
```
Univariate Cox identified [variables] associated with OS (p < 0.05). Multivariable Cox confirmed
independent prognostic value for [Risk score] (HR = [X], 95% CI: [a]–[b], p = [P]), [Stage]
(HR = [X], p = [P]), and [Age] (HR = [X], p = [P]).
```

### Risk score
```
LASSO Cox selected [M] genes from [N] candidates. Risk Score = [formula]. In the training set,
high- vs low-risk HR = [X] (95% CI: [a]–[b], p < 0.001), replicated in external validation
(HR = [Y], 95% CI: [c]–[d], p = [P]).
```

### timeROC
```
Time-dependent ROC showed AUCs at 1, 3, and 5 years of [X], [Y], and [Z], outperforming
[clinical variables] alone (AUC = [a], [b], [c]).
```

### Nomogram
```
The nomogram integrating risk score and independent clinical factors predicted 1-, 3-, and
5-year survival well (C-index = [X], 95% CI: [a]–[b]). Calibration curves showed good agreement
between predicted and observed survival. DCA indicated higher net benefit vs treat-all/treat-none
for threshold probabilities between [X%] and [Y%].
```

## 8. Manuscript Writing Template

### Methods
```
Survival analysis and prognostic model construction:

Kaplan–Meier curves with log-rank tests compared OS between groups. Univariate and multivariate
Cox models identified independent prognostic factors. Proportional hazards were assessed with
Schoenfeld residuals.

A LASSO-penalized Cox model (glmnet, v[version]) selected prognostic genes and constructed a risk
score. The penalty (lambda) was chosen by 10-fold cross-validation. Patients were classified by
the median risk score into high- and low-risk groups.

Time-dependent ROC (timeROC) evaluated 1-, 3-, and 5-year discrimination. A nomogram was built
with rms. Performance used C-index, calibration curves (1000 bootstrap samples), and DCA.

The model was developed in [training cohort] (n=[N]) and validated in [validation cohort] (n=[M]).
Analyses used R (v[version]); two-sided p < 0.05 was significant.
```

## 9. Common Issues and Risks

| Issue | Risk | Mitigation |
|-------|------|------------|
| Small sample (events < 10 per variable) | High | Fewer covariates; LASSO shrinkage |
| Overfitting (train good, test poor) | High | Bootstrap/CV; fewer variables; external validation |
| PH assumption violated | High | Stratified/time-varying Cox; nonparametric alternatives |
| Cutpoint cherry-picking | Medium | surv_cutpoint; confirm in validation |
| Imbalanced train/test | Medium | Stratified split; baseline tables |
| Multiple testing (many KMs) | Medium | Adjust p or label exploratory |
| Immortal time bias | High | Landmark designs |
| Missing covariates | Medium | MICE; sensitivity analyses |
| Ignoring competing risks | Medium | Fine–Gray when non-cancer death is common |
| Optimistic C-index | Medium | Report bootstrap-corrected C-index |

## 10. Experimental Validation Suggestions

### Retrospective validation
```
1. Independent retrospective cohorts
   - In-house archives; public GEO/TCGA/ICGC
   - Validate risk stratification

2. Meta-analysis across datasets
   - Pool HRs; assess heterogeneity

3. Subgroup consistency
   - Stage, age, sex subgroups; define applicability
```

### Prospective validation
```
4. Prospective cohorts
   - Pre-specify risk score at enrollment; follow outcomes
   - Real-world predictive accuracy

5. Decision impact
   - Did risk score change management?
   - Outcomes with intensified treatment in high-risk patients
```

### Mechanistic validation
```
6. Functional validation of signature genes
   - qPCR, WB, IHC for expression
   - Perturbation experiments for biological role

7. Pathway linkage
   - Relate risk score to immunity, EMT, proliferation, etc.
   - Biological interpretation of the signature
```
