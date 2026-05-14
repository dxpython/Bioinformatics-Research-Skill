# Data interpretation template

> Use this checklist to interpret omics tables systematically. Complete a full data assessment before modeling.

---

## 1. Identify the data type

### Decision flow

```
Receive a table
  ├── What is the file format? (.csv / .tsv / .txt / .xlsx / .h5ad / .rds / .mtx)
  ├── What do rows represent? (genes / samples / cells / metabolites / proteins)
  ├── What do columns represent? (samples / genes / variables / time points)
  ├── What are value types? (integer counts / continuous / categorical / mixed)
  └── What is the size? (#rows × #columns)
```

### Quick reference

| Data type | Typical shape | Rows | Columns | Values | Common source |
|-----------|---------------|------|---------|--------|----------------|
| Bulk RNA-seq matrix | ~20k × tens–hundreds | Genes | Samples | Non-negative integers (counts) or non-negative reals (TPM/FPKM) | TCGA / GEO |
| Microarray | ~20k–54k × tens–hundreds | Probes/genes | Samples | Continuous (often log2) | GEO |
| scRNA-seq matrix | ~20k × thousands–10⁵ | Genes | Cells (barcodes) | Sparse integer counts | 10x Genomics |
| Clinical table | tens–hundreds × tens | Patients/samples | Variables | Mixed numeric + categorical | TCGA / hospital |
| Metabolomics | tens–thousands × tens–hundreds | Metabolites | Samples | Non-negative real (area/concentration) | MS platforms |
| Proteomics | thousands × tens–hundreds | Proteins/peptides | Samples | Non-negative intensity | MS platforms |
| Mutations (MAF) | thousands–10⁵ × fixed columns | Mutation events | Fields (gene, position, type, …) | Mixed | TCGA / VCF-derived |

---

## 2. Row and column semantics

### Checklist

```
1. Row IDs
   - Identifier style? (Ensembl / symbol / Entrez / metabolite name / patient ID)
   - Duplicate row names? (must resolve)
   - Row count plausible? (human genes ~20–25k; probes can exceed 50k)

2. Column IDs
   - Column naming? (sample ID / TCGA barcode / GSM ID / custom)
   - Do names encode groups? (e.g. Tumor_01, Normal_01)
   - Duplicate column names?

3. Numeric properties
   - Range: min, max, median
   - Negative values? (possible after log transforms)
   - Missing values? (NA / NaN / blanks)
   - Distribution: approximately normal? need transform?

4. Metadata
   - Sidecar annotation? (gene map, sample groups, platform)
   - Sequencing/microarray platform? (HiSeq, Affymetrix U133, 10x Chromium)
   - Prior processing? (normalization, filtering rules)
```

---

## 3. Fit for the planned analysis

| Check | Pass criteria | If it fails |
|-------|---------------|-------------|
| Sample size | ≥3 per group for DE; total ≥50 for stable survival modeling | Merge cohorts or change aims |
| Units vs method | counts for DESeq2/edgeR; ranked lists for GSEA | Reprocess or convert IDs |
| Gene ID type | Matches tool expectations (e.g. Entrez for some clusterProfiler flows) | Convert with `org.Hs.eg.db` / biomaRt; see **§7 Gene ID conversion pitfalls** before enrichment |
| Group labels | Clear contrasts (tumor vs normal, drug vs control) | Obtain or engineer labels |
| Paired design | Pairing key present when needed | Confirm pairing |
| Batch labels | Required when multi-batch | PCA/harmony/sva as appropriate |
| Survival fields | Time + status for survival endpoints | Cannot run KM/Cox without them |

---

## 4. Information you may still need

```
□ Group membership (case/control or treatment arms)
□ Clinical covariates (age, sex, stage, survival, therapy)
□ Platform / kit / sequencing version
□ Normalization already applied (which algorithm)
□ Gene/probe mapping tables
□ Batch or processing date per sample
□ Provenance of upstream pipeline (tool versions, parameters)
□ QC artifacts (FastQC / MultiQC)
```

---

## 5. Preprocessing guidance

### 5.1 Missing values

| Missing % | Strategy | Notes |
|-----------|----------|-------|
| < 5% | Drop affected rows/columns | Usually minor impact |
| 5–20% | MICE, KNN imputation | Validate imputation assumptions |
| > 20% | Drop feature or sample | Heavy imputation can bias |
| Metabolomics special | Half-minimum or probabilistic minimum | Common for low-abundance features |

### 5.2 Outliers

```
Detection:
  1. Boxplot rule: outside Q3 + 1.5×IQR or Q1 − 1.5×IQR
  2. |Z| > 3
  3. PCA distance from main cloud
  4. Sample–sample correlation far below cohort median

Actions:
  - Rule out technical failure (degraded library, mislabeling)
  - Document exclusions in Methods
  - Sensitivity analysis with/without outliers
```

### 5.3 Normalization choices

| Data | Suggested normalization | R hints |
|------|-------------------------|---------|
| RNA-seq counts (DE) | DESeq2 median ratio or TMM | DESeq2::DESeq(); edgeR::calcNormFactors() |
| RNA-seq (viz/compare) | TPM or vst/rlog | DESeq2::vst() |
| Affymetrix microarray | RMA | affy::rma(); oligo::rma() |
| Two-color microarray | Quantile between arrays | limma::normalizeBetweenArrays() |
| Metabolomics | Pareto / log + auto scaling | MetaboAnalystR workflows |
| Proteomics | Median or quantile normalization | NormalyzerDE, etc. |
| scRNA-seq | SCTransform or LogNormalize | Seurat::SCTransform(); NormalizeData() |

### 5.4 Batch correction

```
Detection:
  1. Color PCA by batch — cluster by batch?
  2. PVCA for variance components

Correction:
  1. ComBat (sva) for known batch with continuous data
  2. ComBat-seq for RNA-seq counts
  3. limma::removeBatchEffect() for visualization only
  4. Harmony / MNN for single-cell integration
```

---

## 6. Analyses enabled by common inputs

| You have | Typical analyses | May also need |
|----------|------------------|---------------|
| Expression + groups | DE, enrichment, GSEA, WGCNA | Nothing else for baseline |
| Expression + survival | Cox, risk scores, KM | Nothing else if fields complete |
| Expression, no groups | Clustering, WGCNA, deconvolution | Phenotype labels improve interpretation |
| scRNA matrix | Clustering, annotation, trajectory, CCC | Optional reference atlases |
| Expression + mutations | Mutation–expression associations | Paired mutation calls |
| Metabolomics + groups | DE metabolites, pathways, PLS-DA | Chemical ID dictionaries |
| Multi-omics matched | MOFA, SNF, iCluster, DIABLO | Strict sample alignment |

---

## 7. Gene ID conversion pitfalls (Symbol / Ensembl / Entrez)

Many enrichment tools (`clusterProfiler`, `enrichplot`, some MSigDB workflows) expect **Entrez Gene IDs** (or a consistent keytype). Converting from symbols or Ensembl is routine but error-prone. Plan for the following **before** you merge mapping results back into expression or DE tables.

### 7.1 One-to-many and duplicate mappings

| Situation | What goes wrong | Practical mitigation |
|-----------|-----------------|----------------------|
| One **gene symbol** → several **Entrez** IDs | `merge()` or `left_join()` duplicates rows or inflates gene counts | Define a rule: e.g. keep **one Entrez per symbol** using NCBI “preferred” / lowest Entrez / **first** after sorting; **document** the rule in Methods. |
| One **Ensembl** → multiple Entrez (history / splits) | Same as above; inconsistent rows across releases | Prefer **Ensembl → Entrez** from the **same GENCODE / Ensembl release** as your quantification; cross-check with current `org.Hs.eg.db` (or species OrgDb). |
| **Duplicate symbols** in the matrix | Ambiguous rownames after conversion | Deduplicate input features first, or map at **Ensembl-ID level** if your quantification is Ensembl-based. |

**Check:** After mapping, report **how many input IDs** became **how many output rows**; flag symbols that expanded to >1 Entrez.

### 7.2 Deprecated, withdrawn, and release mismatch

| Situation | What goes wrong | Practical mitigation |
|-----------|-----------------|----------------------|
| **Withdrawn / merged** Entrez records | Old lists or supplementary tables use IDs no longer in OrgDb | Refresh annotations; use `mapIds(..., multiVals)` awareness; drop or remap via NCBI “replaced by” where possible. |
| **GTF / quant** from assembly A, **OrgDb** from another build | Silent NA or wrong mapping for edge loci | Align **genome build + GENCODE version + Bioconductor OrgDb** (same release cycle when feasible). |
| **Alias symbols** (e.g. historical names) | Maps to NA or wrong gene | Prefer **Ensembl IDs** from the aligner/quantifier as the stable key, then map once to Entrez. |

### 7.3 `bitr()` / `mapIds()` returning `NA`

`clusterProfiler::bitr()` and `AnnotationDbi::mapIds()` return **`NA`** when no row exists in the OrgDb for that identifier.

| Common cause | What to do |
|--------------|------------|
| **Wrong species** OrgDb (e.g. human symbols with `org.Mm.eg.db`) | Confirm organism; use correct `org.*.eg.db` or `EnsDb.*`. |
| **Features not in OrgDb** (many lncRNAs, some pseudogenes, vendor-specific probe names) | Report **mapping rate** (% mapped); exclude or use **ENSEMBL** / custom GMT with matching ID space. |
| **Casing / suffixes** (`CASC11` vs `casc11`, **PAR** genes, readthrough names) | Normalize symbols; strip version from Ensembl (`ENSG00000141510.17` → `ENSG00000141510`) if the keytype is without version. |
| **Mixed ID types** in one column | Split columns; map each type with the correct `fromType` / `keytype`. |

**Minimum QA after conversion:**

```
□ Report: input N, successfully mapped N, NA count, duplicate Entrez count
□ For enrichment: filter to mapped IDs OR impute policy (usually drop NA for ORA)
□ Reconcile row order: DE results joined to Entrez must not duplicate gene rows unless intentional
```

### 7.4 Quick reference (R)

- **`AnnotationDbi::mapIds()`**: control duplicates with `multiVals` (`"first"`, `"filter"`, `"asNA"`); inspect unmapped with `is.na(...)`.
- **`clusterProfiler::bitr()`**: same underlying maps; always check **`sum(is.na(result$ENTREZID))`** (or target column) and **`any(duplicated(result$ENTREZID))`** after resolving symbols.
- Prefer **one stable ID** in your primary assay object (e.g. Ensembl from GTF), then a **single documented** mapping step to Entrez for GO/KEGG.

---

## 8. Interpretation pitfalls

| Pitfall | Issue | Mitigation |
|---------|-------|------------|
| Correlation ≠ causation | Co-expression ≠ regulation | Say “associated”; validate mechanistically |
| p-value fixation | Tiny effects can be “significant” at large N | Report effect sizes + FDR |
| Multiple testing | Testing ~20k genes inflates false positives | FDR (Benjamini–Hochberg) or Bayesian shrinkage |
| Overfitting | Training-set metrics exaggerate generalization | Held-out data, nested CV, bootstrap optimism |
| Data leakage | Feature selection on full data before splitting | Fit selectors inside training folds only |
| Survivorship bias | Analyzing only “complete” cases skews results | Describe missingness; pre-specify exclusions |
| Ecological fallacy | Bulk averages mask single-cell states | Confirm critical claims with scRNA or IHC |

### Short interpretation memo template

```
Data interpretation memo
========================

1. Overview
   - Type: [expression / clinical / metabolomics …]
   - Dimensions: [rows] × [cols]
   - Rows: [genes / patients / metabolites]
   - Columns: [samples / variables]
   - Value range: [min]–[max], median = [median]

2. Quality
   - Missingness: [%], pattern [random / structured]
   - Outliers: [y/n], count = [n]
   - Distribution: [normal / skewed / multimodal]

3. Suitability
   - Appropriate analyses: [list]
   - Not appropriate (and why): [list]
   - Missing metadata: [list]

4. Preprocessing plan
   - Normalization: [choice]
   - Missing data: [strategy]
   - Batch correction: [needed? method?]
   - Gene ID mapping (if enrichment / ORA): input ID type; mapped vs NA count; duplicate Entrez after join; rule for one-symbol→many-Entrez (see §7)

5. Recommended next steps
   - Workflow: [ordered steps]
   - Caveats: [list]
```

---

## Example A: Bulk RNA-seq matrix

### Toy table

```
           Sample_01  Sample_02  Sample_03  Sample_04  Sample_05
TP53           1523       2341       1890       4521       3201
EGFR           8921       7234       9012       2103        891
BRCA1           432        567        489       1023        678
...
```

### Interpretation

| Item | Conclusion |
|------|------------|
| Type | Bulk RNA-seq abundance table |
| Rows | Genes (symbols) |
| Columns | Samples |
| Values | Non-negative integers → likely raw counts |
| Suitable for | DESeq2 DE, enrichment, WGCNA (with design info) |
| Caveats | DESeq2 requires raw counts, not TPM; need group column |
| Preprocessing | Filter low genes (e.g. expressed in ≥50% samples with count >10); DESeq2 size factors |

---

## Example B: Clinical table

### Toy table

```
Patient_ID  Age  Gender  Stage  OS_time  OS_status  Treatment
TCGA-AB-01   56  Male    III    24.5     1          Surgery+Chemo
TCGA-AB-02   43  Female  II     48.2     0          Surgery
TCGA-AB-03   67  Male    IV     12.1     1          Chemo+Immuno
...
```

### Interpretation

| Item | Conclusion |
|------|------------|
| Type | Clinical / outcome table |
| Rows | Patients |
| Columns | Covariates + outcomes |
| Key fields | OS_time (months), OS_status (1 = event per this toy; verify coding) |
| Suitable for | KM, Cox, subgroup analyses |
| Caveats | Confirm censoring code (0/1 meaning); harmonize stage labels |
| Preprocessing | Consistent factor coding; check non-negative times; therapy completeness |

---

## Example C: Metabolomics

### Toy table

```
Metabolite      Control_1  Control_2  Control_3  Disease_1  Disease_2  Disease_3
Glucose          1.23e6     1.45e6     1.12e6     2.34e6     2.56e6     2.11e6
Lactate          5.67e5     4.89e5     5.12e5     8.90e5     9.23e5     8.45e5
Citrate          2.34e5     2.56e5     NA         1.89e5     1.67e5     1.78e5
...
```

### Interpretation

| Item | Conclusion |
|------|------------|
| Type | Untargeted/targeted metabolomics |
| Rows | Metabolites |
| Columns | Samples (names encode group) |
| Values | Non-negative intensities / concentrations |
| Suitable for | DE tests (t/Wilcoxon), PLS-DA, pathway databases via MetaboAnalyst |
| Caveats | Wide dynamic range → log transform; missing value in Control_3 for Citrate |
| Preprocessing | Impute (half-min, KNN), log2 transform, Pareto scaling |

---

*Template v1.0 | For omics table QA and interpretation*
