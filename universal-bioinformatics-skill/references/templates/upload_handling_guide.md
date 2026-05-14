# User File Upload Handling Guide

> How to handle tables, screenshots, code, and manuscript text that users upload during a session.

---

## 1. Table Files (CSV, TSV, XLSX, TXT)

### Step-by-step handling

1. **Determine file type**: Expression matrix, clinical data, differential results, enrichment results, or other
2. **Check data completeness**: Missing values, outliers, column names, row identifiers
3. **Assess suitability**: For the user's intended analysis pipeline
4. **Flag issues**: If data has problems (missing columns, wrong format, not enough samples), point them out and suggest corrections
5. **Recommend methods**: Select appropriate analysis given the data type and quality
6. **Provide code**: Generate runnable code that loads this specific file format

### Common issues to detect

| Issue | Detection | Suggested Fix |
|-------|-----------|---------------|
| Missing column names | First row is numeric | Ask user to add headers |
| Gene IDs not recognized | Cross-ref against Ensembl/HGNC | Suggest ID conversion |
| Mixed data types in columns | `str()` shows factors for numeric cols | Clean and re-import |
| Wrong delimiter | Inconsistent column widths | Try tab/comma auto-detection |
| Encoding issues | Garbled text | Try UTF-8 or locale-specific encoding |

---

## 2. Screenshots / Figure Images

### Step-by-step handling

1. **Identify figure type**: Volcano plot, heatmap, UMAP/t-SNE, KM survival curve, ROC curve, boxplot, barplot, etc.
2. **Describe displayed content**: What each axis represents, which groups are shown, color coding
3. **Extract key information**: Number of significant features, clustering patterns, p-values visible in annotations, survival HR visible in the plot
4. **Provide biological interpretation**: What these results mean in context (tumor biology, pathway activation, immune response, etc.)
5. **Provide manuscript text**: A paragraph ready for the Results section describing the figure
6. **Point out potential improvements**: Color scheme (colorblind-friendly), font size, missing statistical annotations, axis label clarity

### Important caveat

- Cannot extract exact numeric values from raster images — only approximate
- If the figure is low-resolution, some details may be missed — ask user for clarification if needed

---

## 3. Code (R or Python)

### Step-by-step handling

1. **Read code logic**: Understand what the user's code is trying to do
2. **Assess method appropriateness**: Is the chosen method suitable for the data type and research question?
3. **Debug errors**: If the code has errors, analyze root causes and provide fixes
4. **Recommend alternatives**: If the method is inappropriate for the data, point out the issue and suggest better approaches
5. **Optimize quality**: Improve code readability, add comments, fix hardcoded paths, parameterize

### Common code issues to flag

- Hardcoded file paths instead of variables/parameters
- No data validation before analysis
- Missing multiple testing correction
- Using DESeq2 on TPM/FPKM instead of raw counts
- Data leakage in ML (feature selection before train/test split)
- No session info capturing for reproducibility

---

## 4. Manuscript Paragraphs

### Step-by-step handling

1. **Determine section**: Methods, Results, Discussion, Figure Legend, or Reviewer Response
2. **Check description accuracy**: Are the stated results consistent with what was analyzed?
3. **Check statistical completeness**: Are p-values, effect sizes, confidence intervals, sample sizes, and software versions reported?
4. **Provide revision suggestions**: Academic tone improvements, tense corrections, figure citation consistency
5. **Reviewer response mode**: If reviewer comments are included, enter reviewer response mode (see `reviewer_response_template.md`)

### Section-specific checks

| Section | Key Checks |
|---------|-----------|
| Methods | Software + version + parameters + statistical tests + thresholds |
| Results | Figure cross-references, specific values (not just "significant"), logical flow |
| Discussion | Literature comparison, cautious mechanistic interpretation, limitations |
| Figure Legend | Self-contained, axes explained, colors defined, stats annotated |
| Abstract | Word count, key numbers included, matches manuscript content |

---

*Template v1.0*
