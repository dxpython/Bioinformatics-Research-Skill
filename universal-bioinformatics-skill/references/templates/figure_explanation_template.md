# Figure interpretation template

> Systematic guidance for interpreting common bioinformatics figures. Each figure type includes an interpretation checklist and Results-ready paragraph templates (English).

---

## General interpretation framework

Interpret any figure in seven layers:

### 1. What the figure shows (identify plot type)

```
Check:
  - Plot type (scatter / heatmap / bar / box / line / network, etc.)
  - X-axis meaning and units
  - Y-axis meaning and units
  - Color / shape / size encodings
  - Legend content
  - Panel labels (A, B, C, …)
```

### 2. Main pattern

```
Describe:
  - Overall trend (up/down, clusters, separation, flat)
  - Highlighted genes / samples / pathways
  - Direction and magnitude of group differences
  - Unusual or striking patterns
```

### 3. Statistical evidence

```
Report:
  - p-value / adjusted p / FDR
  - Effect size (fold change, Cohen's d, HR, AUC)
  - 95% CI
  - Sample size (n)
  - Statistical test used
```

### 4. Biological meaning

```
Consider:
  - What biological processes or pathways are implicated?
  - Consistency with prior literature?
  - Novel mechanistic hints?
```

### 5. Supported conclusions

```
State:
  - What the data directly support
  - What extra evidence would strengthen claims
  - How this panel fits the overall argument
```

### 6. Avoid over-interpretation

```
Common pitfalls:
  - Correlation ≠ causation
  - Statistical significance ≠ biological importance
  - Single-cohort results may not generalize
  - Bulk averages ≠ single-cell behavior
  - In silico predictions ≠ experimental proof
  - Small samples → cautious language
```

### 7. Results paragraph

> See type-specific templates below.

---

## Example 1: Volcano plot

### Figure features

- X: log2(fold change)
- Y: −log10(adjusted p-value)
- Colors: up (red), down (blue), not significant (gray)
- Thresholds: typically |log2FC| > 1 and padj < 0.05

### Interpretation checklist

| Layer | Content |
|-------|---------|
| Type | Volcano plot summarizing DE |
| Trend | [X] up, [Y] down; balance of directions |
| Labels | Genes at extremes (high |FC| and significance) |
| Stats | |log2FC| > [cutoff], adj.p < [cutoff], [BH/Bonferroni] |
| Biology | Link to enrichment results |
| Caveats | Large |FC| but nonsignificant genes; tiny |FC| but significant hits |

### Results paragraph template

> To identify differentially expressed genes (DEGs) between [group A] and [group B], we performed differential expression analysis using [DESeq2/limma-voom/edgeR]. With the criteria of |log2FC| > [threshold] and adjusted p-value < [threshold], a total of [N] DEGs were identified, including [X] upregulated and [Y] downregulated genes (Fig. [X]A). Among the most significantly upregulated genes were [Gene1], [Gene2], and [Gene3], while [Gene4] and [Gene5] were among the most significantly downregulated. Notably, [Gene1] showed the largest fold change (log2FC = [value], adjusted p = [value]), suggesting its potential role in [biological process].

---

## Example 2: Heatmap

### Figure features

- Rows: genes (often top DEGs or a gene set)
- Columns: samples (often grouped)
- Colors: expression (e.g. red high, blue low; often Z-scored)
- Clustering: row/column dendrograms
- Top annotation: group or clinical covariates

### Interpretation checklist

| Layer | Content |
|-------|---------|
| Type | Heatmap of [top N DEGs / gene set] |
| Trend | Separation by group; blocks of high/low expression |
| Clusters | Co-expression modules; misclassified outliers |
| Biology | Modules may share function |
| Caveats | Z-score colors are not absolute expression; clustering depends on distance/linkage |

### Results paragraph template

> Hierarchical clustering analysis of the top [N] DEGs revealed distinct expression patterns between [group A] and [group B] (Fig. [X]B). The heatmap demonstrated a clear separation of samples based on their group identity, with [group A] samples characterized by upregulation of [gene cluster description] and downregulation of [gene cluster description]. Two major gene clusters were identified: Cluster I (n = [number]) was enriched for genes involved in [biological process 1], while Cluster II (n = [number]) was predominantly associated with [biological process 2]. Notably, [number] out of [total] samples were correctly classified by unsupervised clustering, indicating a robust transcriptional distinction between the two groups.

---

## Example 3: UMAP / t-SNE

### Figure features

- X/Y: embedding coordinates (no physical units)
- Color: cell type, cluster, gene expression, sample
- Each point: one cell

### Interpretation checklist

| Layer | Content |
|-------|---------|
| Type | UMAP/t-SNE of single-cell heterogeneity |
| Trend | [K] clusters; composition shifts by condition |
| Annotation | Clusters mapped to [cell types] via markers |
| Between groups | Proportion changes of selected clusters |
| Biology | Subtype heterogeneity and states |
| Caveats | Distances are not strictly metric; resolution and seed affect layout |

### Results paragraph template

> To delineate the cellular heterogeneity within [tissue/condition], we performed single-cell RNA sequencing on [N] cells from [M] samples. After quality control and filtering, [N'] high-quality cells were retained for downstream analysis. UMAP visualization revealed [K] distinct cell clusters (Fig. [X]A), which were annotated as [cell type 1] (Cluster [N1], [percentage]%), [cell type 2] (Cluster [N2], [percentage]%), and [remaining cell types] based on canonical marker genes (Fig. [X]B). Compared to [control group], [disease group] showed a significant increase in [cell type] proportion ([X]% vs [Y]%, p = [value]) and a decrease in [cell type] proportion ([X]% vs [Y]%, p = [value]) (Fig. [X]C), suggesting a remodeling of the [tissue] cellular composition during [disease/condition].

---

## Example 4: Kaplan–Meier survival curve

### Figure features

- X: time (months/years/days)
- Y: survival probability (0–1 or %)
- Curves: typically high vs low risk
- Number at risk table
- Log-rank p, median survival

### Interpretation checklist

| Layer | Content |
|-------|---------|
| Type | KM curves comparing survival |
| Trend | Separation; worse survival in high-risk arm |
| Stats | Log-rank p; median survival; HR and 95% CI |
| Landmarks | 1/3/5-year survival rates |
| Biology | Link grouping variable to outcome |
| Caveats | Heavy censoring; unstable tails at low n; validate externally |

### Results paragraph template

> To evaluate the prognostic value of the [risk score / gene expression / molecular subtype], patients were stratified into high-risk (n = [N1]) and low-risk (n = [N2]) groups based on the [median / optimal cutoff] of [variable]. Kaplan-Meier analysis demonstrated that patients in the high-risk group had significantly worse overall survival (OS) compared to the low-risk group (median OS: [X] months vs [Y] months; HR = [value], 95% CI: [lower]-[upper]; log-rank p [= / <] [value]) (Fig. [X]A). The 1-year, 3-year, and 5-year survival rates were [X]%, [X]%, and [X]% in the high-risk group versus [Y]%, [Y]%, and [Y]% in the low-risk group, respectively. This survival difference was consistently observed in the validation cohort ([dataset name], log-rank p = [value]) (Fig. [X]B).

---

## Example 5: ROC curve

### Figure features

- X: 1 − specificity (FPR)
- Y: sensitivity (TPR)
- Diagonal: AUC = 0.5 reference
- AUC annotated; possibly multiple curves or time points

### Interpretation checklist

| Layer | Content |
|-------|---------|
| Type | ROC for [classifier / biomarker / model] |
| Metric | AUC: 0.5 random; 0.7–0.8 acceptable; 0.8–0.9 good; >0.9 excellent (check overfitting) |
| Compare | Models or time points; DeLong for AUC differences |
| Cutoff | If shown, sensitivity/specificity at best threshold |
| Caveats | Class imbalance; training-set AUC optimistic |

### Results paragraph template

> The discriminative ability of the [risk model / biomarker] was assessed using time-dependent ROC analysis. The AUC values for predicting [1-year / 3-year / 5-year] [overall survival / disease-free survival] were [value], [value], and [value], respectively (Fig. [X]A), indicating [acceptable / good / excellent] predictive performance. Compared with [clinical stage / other existing models], our [model name] demonstrated superior predictive accuracy (AUC: [value] vs [value], DeLong test p = [value]) (Fig. [X]B). In the external validation cohort ([dataset name]), the model maintained [good / stable] performance with AUC values of [value], [value], and [value] for [1-year / 3-year / 5-year] prediction, respectively (Fig. [X]C).

---

## Appendix: Quick notes for other common plots

### Bubble / dot plot (enrichment)

```
Check:
  - Bubble size = gene count
  - Color = p or adjusted p
  - X = gene ratio / rich factor / fold enrichment
  - Pathways sorted by significance

Key Results sentence:
  "GO enrichment analysis revealed that upregulated DEGs were significantly
   enriched in biological processes related to [process 1], [process 2],
   and [process 3] (Fig. [X]). KEGG pathway analysis identified [pathway 1]
   (adjusted p = [value]) and [pathway 2] (adjusted p = [value]) as the
   most significantly enriched pathways."
```

### Box plot (group comparison)

```
Check:
  - Medians and spread (IQR)
  - Whiskers and outliers
  - Significance stars (* / ** / *** / ns)
  - Test used (t-test / Wilcoxon / ANOVA / Kruskal–Wallis)

Key Results sentence:
  "The expression level of [Gene] was significantly higher in [group A]
   compared to [group B] (median: [value] vs [value], Wilcoxon test
   p = [value]) (Fig. [X])."
```

### Forest plot (Cox regression)

```
Check:
  - HR and 95% CI per variable
  - CI crossing 1 → not significant
  - p-values
  - Consistency univariate vs multivariable

Key Results sentence:
  "Univariate Cox regression analysis identified [variable 1] (HR = [value],
   95% CI: [range], p = [value]), [variable 2], and [variable 3] as
   significant prognostic factors. In multivariate Cox regression,
   [variable 1] (HR = [value], 95% CI: [range], p = [value]) and
   [variable 2] (HR = [value], 95% CI: [range], p = [value]) remained
   as independent prognostic predictors after adjusting for clinicopathological
   variables (Fig. [X])."
```

### Correlation scatter plot

```
Check:
  - Trend (positive / negative / none)
  - Regression line and band
  - Pearson r or Spearman rho and p
  - Outlier leverage

Key Results sentence:
  "[Variable A] showed a significant positive/negative correlation with
   [Variable B] (Spearman rho = [value], p = [value]) (Fig. [X]),
   suggesting that [biological interpretation]."
```

---

*Template v1.0 | Covers volcano, heatmap, UMAP/t-SNE, survival, ROC, and common companion plots*
