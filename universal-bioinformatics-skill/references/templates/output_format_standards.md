# Output Format Standards

> Standard output formats for code, result interpretation, and manuscript writing responses.

---

## 1. Code Output

- R code in ```r blocks
- Python code in ```python blocks
- Shell commands in ```bash blocks
- Comments in English

## 2. Result Interpretation Output

Structured output using this format:

```
## Analysis Result Interpretation

### 1. Analysis Overview
[Brief description of the analysis performed]

### 2. Key Findings
[List key findings]

### 3. Statistical Results
[Specific values and statistics]

### 4. Biological Significance
[Biological interpretation]

### 5. Limitations
[Limitations of this analysis]

### 6. Recommendations
[Follow-up analysis or experimental validation suggestions]

### 7. Manuscript Text Reference
[Paragraphs ready for manuscript use]
```

## 3. Manuscript Writing Output

- Use academic English
- Label software versions and parameters
- Correct citation formatting
- Consistent figure/table referencing

---

## 4. Figure Interpretation Output

### Standard flow

1. **Figure type identification** (volcano, heatmap, UMAP, KM curve, ROC, etc.)
2. **Key observations** (describe what the figure shows)
3. **Statistical evidence** (extract p-values, HR, AUC, fold changes from annotations)
4. **Biological interpretation** (what these results mean biologically)
5. **Manuscript-ready Results text** (a paragraph for the Results section)
6. **Improvement suggestions** (color scheme, axis labels, statistical annotations)

---

## 5. Code Generation Output

### Environment specification

- R code must include Rlib preamble (see `references/environment/R_environment.md`)
- Python code must use argparse and pathlib; designed for `uv run python`
- All parameters at script top with comments explaining valid ranges

### Input/output clarity

- Specify input file format (delimiter, column names, encoding)
- Specify all output files and their formats
- Include session info (R: `sessionInfo()`, Python: `sys.version` + package versions)

### Error handling

- Include common error messages and their fixes in code comments
- For R: use `tryCatch` around error-prone steps
- For Python: use try/except with informative messages

---

*Template v1.0*
