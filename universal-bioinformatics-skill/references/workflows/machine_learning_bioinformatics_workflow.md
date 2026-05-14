# Machine Learning Bioinformatics Workflow

## 1. Applicable Scenarios

- Diagnostic biomarker discovery
- Prognostic prediction model construction
- Treatment response prediction
- Multi-gene signature development
- Feature selection from high-dimensional omics data
- Model comparison and benchmarking

## 2. Input Data Requirements

| Input | Description | Format |
|-------|-------------|--------|
| Feature matrix | Genes/proteins/metabolites × Samples | CSV/TSV |
| Labels | Binary (disease/normal) or continuous (survival) | CSV/TSV |
| Validation set (strongly recommended) | Independent dataset for testing | Same format |

## 3. Data Validation Checklist

Before modeling, verify:
```
□ Feature matrix is complete (no NA, no Inf)
□ All features are numeric (encode categorical if needed)
□ Sample size is adequate: at least [10 × number of features] after feature selection
□ Class balance: minority class ≥ 10–15% of total for classification
□ Labels and features are aligned by sample ID
□ No constant or nearly-zero-variance features
□ Training/testing split is stratified by key variables (group, batch, outcome)
□ For survival: at least [number of features] events in the training set
□ External validation set is truly independent (different source/cohort)
```

| Issue | Detection | Action |
|-------|-----------|--------|
| Missing values | `is.na()` / `df.isnull().sum()` | Remove features with > 20% missing; impute rest |
| Class imbalance | `table(labels)` | SMOTE, undersampling, or class weights |
| Data leakage | Feature selection on full dataset | Re-do selection within CV only |
| High correlation | `cor()` > 0.9 between features | Remove one of each pair |
| Outliers | Boxplots, Z-score > 3 | Winsorize or remove |
| Wrong data type | `str()` / `df.dtypes` | Convert factors/strings to numeric |

## 4. Recommended Analysis Steps

### Step 1: Data Preparation
1. Handle missing values (impute or remove)
2. Remove zero-variance features
3. Log-transform if needed
4. Split data: training (70%) / testing (30%) — **before any feature selection**
5. For class imbalance: SMOTE, undersampling, or class weights

**Critical**: Feature selection MUST be done within training data only. Data leakage occurs if test data influences feature selection.

### Step 2: Feature Selection

| Method | Type | R Package | Python Package |
|--------|------|-----------|----------------|
| LASSO | Regularization | glmnet | sklearn |
| SVM-RFE | Wrapper | e1071/caret | sklearn |
| Random Forest Importance | Embedded | randomForest | sklearn |
| Boruta | Wrapper | Boruta | boruta_py |
| XGBoost Importance | Embedded | xgboost | xgboost |
| mRMR | Filter | mRMRe | mrmr_selection |

**Recommended strategy**: Use 2-3 methods and take intersection or union, then validate.

### Step 3: Model Training

| Model | Strengths | Best For |
|-------|-----------|----------|
| Logistic Regression | Interpretable, fast | Baseline model |
| LASSO Logistic | Built-in feature selection | High-dimensional data |
| SVM | Good with small samples | Small-medium datasets |
| Random Forest | Robust, handles nonlinearity | General purpose |
| XGBoost | High performance | Large datasets |
| LightGBM | Fast, memory-efficient | Very large datasets |
| CatBoost | Handles categorical features | Mixed feature types |

### Step 4: Cross-Validation
- 5-fold or 10-fold cross-validation (repeated 3-5 times)
- Stratified folds to maintain class balance
- Report mean ± SD of metrics across folds
- Nested cross-validation for unbiased estimates (outer: evaluation, inner: hyperparameter tuning)

### Step 5: Hyperparameter Tuning
- Grid search or random search
- For LASSO: lambda selection via cv.glmnet (lambda.min or lambda.1se)
- For Random Forest: ntree, mtry
- For XGBoost/LightGBM: learning_rate, max_depth, n_estimators, subsample

### Step 6: Model Evaluation

| Metric | When to Use |
|--------|-------------|
| AUC-ROC | Overall discrimination (primary metric) |
| Accuracy | Balanced datasets only |
| Sensitivity/Specificity | When false negatives/positives have different costs |
| Precision/Recall | Imbalanced datasets |
| F1 Score | Balance of precision and recall |
| PR-AUC | Imbalanced datasets (more informative than ROC-AUC) |
| Calibration | Probability estimates reliability |

### Step 7: External Validation
- Apply trained model to completely independent dataset
- Report same metrics as training
- Expected: some performance drop (AUC ~0.05-0.15 lower)
- If large drop: model may be overfitting

### Step 8: Model Interpretation

| Method | What It Shows |
|--------|---------------|
| Feature importance | Which features contribute most |
| SHAP values | Feature contribution direction and magnitude |
| Partial dependence plots | Feature effect on prediction |
| LIME | Local instance explanations |

### Step 9: Clinical Utility
- Decision curve analysis
- Net reclassification improvement (NRI)
- Integrated discrimination improvement (IDI)
- Compare with existing clinical models

## 5. Recommended Statistical Methods

| Analysis | Method | Key Metric / Threshold |
|----------|--------|------------------------|
| Feature selection stability | Bootstrap resampling of LASSO/RFE | Selection frequency > 50% |
| Model comparison | DeLong test or bootstrap | AUC comparison p-value |
| Classification | AUC-ROC (primary), sensitivity, specificity, PPV, NPV | AUC > 0.7 acceptable; > 0.8 good |
| Calibration | Hosmer-Lemeshow test, calibration curve | p > 0.05 (well-calibrated) |
| Survival models | C-index, time-dependent AUC | C-index > 0.65 acceptable |
| Overfitting check | Train vs. test AUC gap | Gap < 0.1 acceptable |
| Feature importance robustness | SHAP stability across CV folds | Consistent direction and ranking |
| Decision curve analysis | Net benefit across threshold probabilities | Net benefit > treat-all and treat-none |

- Nested cross-validation: outer loop for evaluation, inner loop for hyperparameter tuning — never tune on the outer test fold
- Report 95% confidence intervals for AUC via bootstrap (2000 iterations)
- Model interpretability is as important as performance: always provide SHAP or feature importance plots
- External validation in an independent dataset is mandatory before claiming clinical utility

## 6. Recommended Figures

| Figure | Description |
|--------|-------------|
| ROC curves | All models on same plot + AUC values |
| PR curves | For imbalanced data |
| Calibration plot | Predicted vs. observed probability |
| LASSO coefficient path | Lambda vs. coefficients |
| Feature importance barplot | Top features from RF/XGBoost |
| SHAP summary plot | Feature contributions |
| SHAP dependence plot | Individual feature effects |
| Model comparison table | Metrics across all models |
| Decision curve | Clinical net benefit |
| Confusion matrix heatmap | Classification performance |


## 7. Result Interpretation Template

> To identify diagnostic biomarkers for [disease], we applied [N] machine learning algorithms to the training cohort ([N1] samples). LASSO regression selected [K] features at the optimal lambda (lambda.1se). Random forest importance confirmed [X] of these features among its top 20. A [K]-feature [model] achieved an AUC of X.XX (95% CI: X.XX-X.XX) in the training set and X.XX (95% CI: X.XX-X.XX) in the independent validation set (Figure XA). SHAP analysis revealed [Feature1] and [Feature2] as the most influential predictors (Figure XB). Decision curve analysis demonstrated that the model provided greater net benefit than the treat-all or treat-none strategies across threshold probabilities of X% to Y% (Figure XC).

## 8. Manuscript Writing Template

### Methods

> **Machine learning-based [diagnostic / prognostic] model construction**
>
> A total of [N] samples with [N_features] features were randomly split into a training set ([70%]) and an internal test set ([30%]) using stratified sampling to preserve the outcome distribution. Feature selection was performed within 10-fold cross-validation on the training set only, using [LASSO regression (glmnet vX.X) / SVM-RFE (e1071 vX.X) / Boruta (Boruta vX.X)]. Features selected by at least [2/3] methods were retained. Five machine learning algorithms were trained and compared: [logistic regression / LASSO / random forest (randomForest vX.X) / SVM (e1071 vX.X) / XGBoost (xgboost vX.X)]. Hyperparameters were tuned using 5-fold cross-validation (inner loop of nested CV). Model performance was evaluated by area under the receiver operating characteristic curve (AUC), sensitivity, specificity, positive predictive value (PPV), and negative predictive value (NPV). The best-performing model was externally validated in [validation_cohort, n = N] from [GEO dataset GSEXXXXX / independent hospital cohort]. Model calibration was assessed using calibration curves and the Hosmer-Lemeshow test. SHAP values were computed to interpret feature contributions. Decision curve analysis was performed to assess clinical net benefit. All statistical analyses were two-sided with p < 0.05 considered significant, and 95% confidence intervals for AUC were estimated by bootstrap (2000 iterations).

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| Data leakage | Feature selection only within training set |
| Overfitting | Use cross-validation; external validation essential |
| Class imbalance | Use SMOTE, class weights, or PR-AUC instead of ROC-AUC |
| Too many features, too few samples | Aggressive feature selection; regularization |
| Non-reproducible results | Set random seeds; report full methodology |
| No external validation | Clearly state as major limitation |
| Black-box model | Use SHAP/LIME for interpretability |

## 10. Experimental Validation Suggestions

1. **Prospective cohort**: Validate model predictions prospectively
2. **Multi-center validation**: Test across different populations
3. **qRT-PCR panel**: Develop clinical-grade assay for selected biomarkers
4. **Functional validation**: Validate biological role of top features
