# Code generation template

> Use this template to standardize bioinformatics analysis code in R or Python: inputs/outputs, environment, parameters, and error handling.

---

## 1. Input file specifications

### Common input types

Before generating code, define input formats clearly.

| File type | Format | Example |
|-----------|--------|---------|
| Expression matrix | CSV/TSV; rows = genes, columns = samples; first column or row names = gene IDs | `gene,sample1,sample2,...` |
| Clinical table | CSV/TSV; rows = patients; columns include ID, OS_time, OS_status, etc. | `patient_id,age,gender,stage,os_time,os_status` |
| Gene list | TXT; one gene per line; with or without header | `TP53\nEGFR\nBRCA1` |
| Sample groups | CSV/TSV; two columns: sample_id, group | `sample_id,group` |
| Single-cell | 10x Genomics (`matrix.mtx` + `barcodes.tsv` + `features.tsv`) or `.h5` | Standard 10x output folder |
| GMT | Gene sets; tab-separated | `pathway_name\tdescription\tgene1\tgene2\t...` |

### Input checklist

```
Before running or shipping code, confirm:
  □ File paths are correct
  □ Delimiter (comma / tab / space)
  □ Header present or not
  □ Row names present or not (R)
  □ Encoding (UTF-8 preferred; avoid GBK unless required)
  □ Missing values: NA / NaN / blank / literal "NA"
  □ Gene ID type: symbol / Ensembl / Entrez
  □ Whether values are already normalized
```

---

## 2. Runtime environment

### R

```
R environment:
  - R version: >= 4.2.0 (4.3.x or 4.4.x recommended)
  - Packages: project-local library ./Rlib/ via .libPaths()
  - Install:
      install.packages("pkg", lib = "./Rlib")
      BiocManager::install("pkg", lib = "./Rlib")

  At the top of scripts:
      .libPaths(c("./Rlib", .libPaths()))
```

### Python

```
Python environment:
  - Python: >= 3.9 (3.10 or 3.11 recommended)
  - Env: uv (preferred) or pip + venv
  - Create env:
      uv venv .venv
      source .venv/bin/activate   # macOS/Linux
      .venv\Scripts\activate      # Windows

  - Install:
      uv pip install package_name
      # or
      uv add package_name         # with pyproject.toml

  Optional version check:
      import sys
      print(f"Python version: {sys.version}")
```

---

## 3. Full script skeletons

### R script structure

```r
# ============================================================
# Script: [analysis_name].R
# Description: [one-line purpose]
# Author: [name]
# Date: [date]
# R version: [R version]
# ============================================================

# ---- 0. Environment ----
.libPaths(c("./Rlib", .libPaths()))

required_packages <- c("pkg1", "pkg2", "pkg3")
bioc_packages <- c("bioc_pkg1", "bioc_pkg2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, lib = "./Rlib", repos = "https://cloud.r-project.org")
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = "./Rlib")
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, lib = "./Rlib", update = FALSE)
  }
}

library(pkg1)
library(pkg2)
library(pkg3)

# ---- 1. Parameters (user-editable) ----
input_expression <- "data/expression_matrix.csv"   # expression matrix path
input_clinical   <- "data/clinical_data.csv"        # clinical table path
output_dir       <- "results/"                      # output directory
log2fc_threshold <- 1                               # log2FC cutoff
padj_threshold   <- 0.05                            # adjusted p-value cutoff
fig_width        <- 8                               # figure width (inches)
fig_height       <- 6                               # figure height (inches)
fig_dpi          <- 300                             # figure resolution

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---- 2. Load and check data ----
cat(">> Loading data...\n")

expr <- read.csv(input_expression, row.names = 1, check.names = FALSE)
clinical <- read.csv(input_clinical)

cat(sprintf("   Expression matrix: %d genes x %d samples\n", nrow(expr), ncol(expr)))
cat(sprintf("   Clinical table: %d patients x %d variables\n", nrow(clinical), ncol(clinical)))

stopifnot("NA in expression matrix" = !any(is.na(expr)))
stopifnot("clinical missing OS_time" = "OS_time" %in% colnames(clinical))
stopifnot("clinical missing OS_status" = "OS_status" %in% colnames(clinical))

# ---- 3. Preprocessing ----
cat(">> Preprocessing...\n")

# [preprocessing code]

# ---- 4. Core analysis ----
cat(">> Running analysis...\n")

# [analysis code]

# ---- 5. Plots ----
cat(">> Plotting...\n")

# [plotting code]

ggsave(
  filename = file.path(output_dir, "figure_name.pdf"),
  plot = p,
  width = fig_width,
  height = fig_height,
  dpi = fig_dpi
)

# ---- 6. Save outputs ----
cat(">> Saving results...\n")

write.csv(result_table, file.path(output_dir, "result_table.csv"), row.names = FALSE)

# ---- 7. Session info ----
cat(">> Done.\n")
cat(sprintf("   Outputs: %s\n", output_dir))
sink(file.path(output_dir, "session_info.txt"))
sessionInfo()
sink()
```

### Python script structure

```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script: [analysis_name].py
Description: [one-line purpose]
Author: [name]
Date: [date]
Python: [version]
"""

# ---- 0. Imports ----
import os
import sys
import warnings
import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)

plt.rcParams["figure.dpi"] = 300
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["font.size"] = 12
sns.set_theme(style="whitegrid")
warnings.filterwarnings("ignore")


# ---- 1. CLI parameters ----
def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="[analysis description]")

    parser.add_argument("--input_expression", type=str, required=True,
                        help="Path to expression matrix (CSV/TSV)")
    parser.add_argument("--input_clinical", type=str, required=True,
                        help="Path to clinical table (CSV/TSV)")
    parser.add_argument("--output_dir", type=str, default="results/",
                        help="Output directory (default: results/)")
    parser.add_argument("--log2fc_threshold", type=float, default=1.0,
                        help="log2FC threshold (default: 1.0)")
    parser.add_argument("--padj_threshold", type=float, default=0.05,
                        help="Adjusted p-value threshold (default: 0.05)")
    parser.add_argument("--fig_width", type=float, default=8.0,
                        help="Figure width in inches (default: 8.0)")
    parser.add_argument("--fig_height", type=float, default=6.0,
                        help="Figure height in inches (default: 6.0)")

    return parser.parse_args()


# ---- 2. Load and check ----
def load_data(input_expression: str, input_clinical: str) -> tuple:
    """Load expression matrix and clinical table."""
    logger.info("Loading data...")

    expr = pd.read_csv(input_expression, index_col=0, sep=None, engine="python")
    clinical = pd.read_csv(input_clinical, sep=None, engine="python")

    logger.info(f"  Expression: {expr.shape[0]} genes x {expr.shape[1]} samples")
    logger.info(f"  Clinical: {clinical.shape[0]} patients x {clinical.shape[1]} columns")

    assert not expr.isnull().any().any(), "NA in expression matrix"
    assert "OS_time" in clinical.columns, "clinical missing OS_time"
    assert "OS_status" in clinical.columns, "clinical missing OS_status"

    return expr, clinical


# ---- 3. Preprocess ----
def preprocess(expr: pd.DataFrame, clinical: pd.DataFrame) -> tuple:
    """Preprocess data."""
    logger.info("Preprocessing...")

    # [preprocessing code]

    return expr_processed, clinical_processed


# ---- 4. Analysis ----
def run_analysis(expr: pd.DataFrame, clinical: pd.DataFrame, **params) -> pd.DataFrame:
    """Core analysis."""
    logger.info("Running analysis...")

    # [analysis code]

    return results


# ---- 5. Plots ----
def plot_results(results: pd.DataFrame, output_dir: str, **params):
    """Generate figures."""
    logger.info("Plotting...")

    fig, ax = plt.subplots(figsize=(params.get("fig_width", 8),
                                     params.get("fig_height", 6)))

    # [plotting code]

    fig.savefig(
        os.path.join(output_dir, "figure_name.pdf"),
        bbox_inches="tight"
    )
    plt.close(fig)
    logger.info(f"  Saved: {output_dir}/figure_name.pdf")


# ---- 6. Main ----
def main():
    args = parse_args()

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    expr, clinical = load_data(args.input_expression, args.input_clinical)
    expr, clinical = preprocess(expr, clinical)
    results = run_analysis(
        expr, clinical,
        log2fc_threshold=args.log2fc_threshold,
        padj_threshold=args.padj_threshold,
    )
    plot_results(
        results, args.output_dir,
        fig_width=args.fig_width,
        fig_height=args.fig_height,
    )
    results.to_csv(
        os.path.join(args.output_dir, "result_table.csv"),
        index=False
    )

    logger.info("Done.")
    logger.info(f"  Outputs: {args.output_dir}")


if __name__ == "__main__":
    main()
```

---

## 4. Parameter documentation

Document every user-tunable parameter at the top of the script or in a separate file.

| Parameter | Type | Default | Description | Notes |
|-----------|------|---------|-------------|-------|
| `input_expression` | str | required | Path to expression matrix | CSV/TSV; rows = genes, columns = samples |
| `input_clinical` | str | required | Path to clinical table | Must include OS_time and OS_status |
| `output_dir` | str | `"results/"` | Output folder | Created if missing |
| `log2fc_threshold` | float | `1` | log2FC cutoff for DE | Common: 0.585 (~1.5×), 1 (2×), 2 (4×) |
| `padj_threshold` | float | `0.05` | Adjusted p-value cutoff | Stricter analyses: 0.01 |
| `fig_width` | float | `8` | Figure width (in) | Single column ~3.5 in; double ~7 in |
| `fig_height` | float | `6` | Figure height (in) | Adjust per plot type |
| `fig_dpi` | int | `300` | Resolution | Journals often require ≥ 300 DPI |

---

## 5. Output files

### Example output tree

```
results/
  ├── result_table.csv          # Main table (e.g., DE genes)
  ├── figure_volcano.pdf
  ├── figure_heatmap.pdf
  ├── figure_survival_KM.pdf
  ├── figure_ROC.pdf
  ├── enrichment_GO.csv
  ├── enrichment_KEGG.csv
  ├── risk_score.csv
  ├── model_coefficients.csv
  ├── session_info.txt          # R sessionInfo()
  └── log.txt                   # optional run log
```

### Format guidelines

| Type | Format | Requirements |
|------|--------|----------------|
| Tables | CSV, UTF-8 | Headers; sensible numeric precision; avoid messy scientific notation in Excel |
| Publication figures | PDF or TIFF | ≥ 300 DPI; embedded fonts; CMYK if journal requires |
| Preview figures | PNG | 150–300 DPI; RGB |
| Logs | TXT | Timestamp, parameters, key statistics |

---

## 6. Figure captions from code

For each figure, document in comments or a sidecar:

```
Figure note template:
  - Title: [short title]
  - Type: [volcano / heatmap / KM / ...]
  - X-axis: [meaning and units]
  - Y-axis: [meaning and units]
  - Colors: [what they map to]
  - Main finding: [one sentence]
  - Caveats: [what to watch when interpreting]
```

---

## 7. Common errors and fixes

### R

| Error | Likely cause | Fix |
|-------|--------------|-----|
| `there is no package called 'xxx'` | Package missing | `install.packages("xxx", lib = "./Rlib")` or `BiocManager::install` |
| `cannot open the connection` | Bad path / missing file | Check path; `file.exists()` |
| `differing number of rows` | Merge/join mismatch | Check dimensions; `merge()` on ID |
| `cannot allocate vector of size X Gb` | Memory | Subset data; `gc()`; more RAM |
| `argument is of length zero` | Empty vector in `if` | Check `length(x) > 0` |
| `contrasts ... only to factors with 2 or more levels` | One level in factor | Fix grouping column; drop single-level groups |
| `glm.fit: algorithm did not converge` | Separation / sparse data | Check for quasi-separation; adjust model / iterations |

### Python

| Error | Likely cause | Fix |
|-------|--------------|-----|
| `ModuleNotFoundError` | Missing package | `uv pip install xxx` |
| `FileNotFoundError` | Bad path | `os.path.exists()` |
| `could not convert string to float` | Non-numeric cells | `pd.to_numeric(errors="coerce")` |
| `MemoryError` | Large data | dtypes; chunked `read_csv` |
| `KeyError` | Wrong column name | `df.columns.tolist()`; watch spaces/case |
| `SettingWithCopyWarning` | Chained assignment | `.loc[...] =` or `.copy()` |
| `CUDA out of memory` | GPU OOM | Smaller batch; `torch.cuda.empty_cache()` |

---

## 8. Quick reference: tunable analysis parameters

### Differential expression

| Parameter | Meaning | Typical range | Default |
|-----------|---------|---------------|---------|
| `log2fc_threshold` | Fold-change cutoff | 0.585–2 | 1 |
| `padj_threshold` | Adjusted p | 0.01–0.1 | 0.05 |
| `min_count` | Low-count filter | 5–20 | 10 |
| `min_samples` | Min samples expressing | 10%–50% of total | context-dependent |

### Enrichment / GSEA

| Parameter | Meaning | Typical range | Default |
|-----------|---------|---------------|---------|
| `pvalueCutoff` | Enrichment p | 0.01–0.1 | 0.05 |
| `qvalueCutoff` | q-value | 0.05–0.25 | 0.2 |
| `minGSSize` | Min genes in set | 5–15 | 10 |
| `maxGSSize` | Max genes in set | 200–1000 | 500 |
| `nPerm` | GSEA permutations | 1000–10000 | 1000 |

### Survival

| Parameter | Meaning | Typical range | Default |
|-----------|---------|---------------|---------|
| `cutoff_method` | Risk stratification | median / optimal / tertile | median |
| `min_follow_up` | Min follow-up (days) | 0–90 | 30 |
| `time_points` | ROC time points (years) | 1, 3, 5 | `c(1, 3, 5)` |

### LASSO / elastic net

| Parameter | Meaning | Typical range | Default |
|-----------|---------|---------------|---------|
| `alpha` | EN mix (0=Ridge, 1=LASSO) | 0–1 | 1 |
| `nfolds` | CV folds | 5–10 | 10 |
| `lambda` | Strength | `lambda.min` / `lambda.1se` | `lambda.min` |

### Single-cell (Seurat-like)

| Parameter | Meaning | Typical range | Default |
|-----------|---------|---------------|---------|
| `min_genes` | Min genes per cell | 200–500 | 200 |
| `max_genes` | Max genes per cell | 2500–8000 | data-dependent |
| `max_mt_percent` | Max MT% | 5%–20% | 10% |
| `n_variable_features` | HVG count | 2000–5000 | 3000 |
| `n_pcs` | PCs used | 15–50 | 30 |
| `resolution` | Clustering resolution | 0.1–2.0 | 0.5 |

---

*Template v1.0 | Standardized R/Python bioinformatics scripts*
