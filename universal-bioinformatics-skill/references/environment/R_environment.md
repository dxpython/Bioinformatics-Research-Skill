# R Environment Guide: Rlib Local Package Cache

## Why Use Rlib?

Traditional R package management installs packages to a system-wide library (e.g., `/usr/local/lib/R/site-library`). This causes several problems:

1. **Version conflicts**: Different projects may need different package versions
2. **Permission issues**: System library may require admin/root privileges
3. **Reproducibility**: Hard to track exactly which packages a project uses
4. **Portability**: Cannot easily move the project to another machine

The **Rlib local package cache** solves all of these by installing packages into a project-level `Rlib/` directory, similar to Python's virtual environments or Node.js's `node_modules/`.

## Quick Start

### Step 1: Initialize Rlib

Run the following in R:

```r
# Create Rlib directory and configure .libPaths()
rlib_path <- file.path(getwd(), "Rlib")
if (!dir.exists(rlib_path)) dir.create(rlib_path, recursive = TRUE)
.libPaths(c(rlib_path, .libPaths()))

# Optionally, write to .Rprofile for automatic loading
rprofile_path <- ".Rprofile"
if (!file.exists(rprofile_path)) {
  cat('local({\n',
      '  rlib_path <- file.path(getwd(), "Rlib")\n',
      '  if (dir.exists(rlib_path)) {\n',
      '    .libPaths(c(rlib_path, .libPaths()))\n',
      '  }\n',
      '})\n',
      file = rprofile_path, sep = "")
  cat("Created .Rprofile for automatic Rlib loading.\n")
}

cat("Rlib initialized at:", rlib_path, "\n")
cat(".libPaths()[1]:", .libPaths()[1], "\n")
```

This will:
- Create the `Rlib/` directory in your project root
- Add `Rlib/` to `.libPaths()` (highest priority)
- Create/update `.Rprofile` so future R sessions automatically use Rlib

### Step 2: Install Required Packages

Run the following in R (make sure Step 1 was run first to set up Rlib):

```r
# Helper: install package if not already available
install_if_missing <- function(pkg, bioc = FALSE, lib = NULL) {
  if (is.null(lib)) lib <- .libPaths()[1]
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
      BiocManager::install(pkg, lib = lib, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org")
    }
  }
  cat(sprintf("  %-35s: %s\n", pkg,
    if (requireNamespace(pkg, quietly = TRUE)) "OK" else "FAILED"))
}

# List of CRAN packages to install
cran_packages <- c(
  "tidyverse", "data.table", "ggplot2", "pheatmap", "RColorBrewer",
  "ggrepel", "patchwork", "survival", "survminer", "glmnet",
  "randomForest", "e1071", "caret", "pROC", "timeROC", "rms",
  "ggpubr", "corrplot", "igraph", "ggraph", "FactoMineR", "factoextra",
  "VennDiagram", "UpSetR"
)

# List of Bioconductor packages to install
bioc_packages <- c(
  "DESeq2", "edgeR", "limma", "clusterProfiler",
  "org.Hs.eg.db", "org.Mm.eg.db", "enrichplot", "ReactomePA",
  "GSVA", "GSEABase", "msigdbr", "AnnotationDbi", "STRINGdb",
  "maftools", "TCGAbiolinks", "SummarizedExperiment",
  "SingleCellExperiment", "scater", "scran", "DOSE", "ComplexHeatmap",
  "WGCNA"
)

cat("Installing CRAN packages...\n")
for (pkg in cran_packages) install_if_missing(pkg, bioc = FALSE)

cat("\nInstalling Bioconductor packages...\n")
for (pkg in bioc_packages) install_if_missing(pkg, bioc = TRUE)

cat("\nDone.\n")
```

This will:
- Install all common CRAN and Bioconductor packages to `Rlib/`
- Skip packages that are already installed
- Continue even if individual packages fail

### Step 3: Verify Installation

```r
# List all installed packages in Rlib
rlib_pkgs <- installed.packages(lib.loc = file.path(getwd(), "Rlib"))
cat(sprintf("Packages in Rlib: %d\n", nrow(rlib_pkgs)))
print(rownames(rlib_pkgs))
```

## Loading Rlib in Your Scripts

Every R analysis script should begin with:

```r
# === Load Rlib local package cache ===
project_root <- getwd()
rlib_path <- file.path(project_root, "Rlib")
if (dir.exists(rlib_path)) {
  .libPaths(c(rlib_path, .libPaths()))
}
```

If you have `.Rprofile` configured (from Step 1), this happens automatically.

## Installing Additional Packages

First, ensure Rlib is loaded:
```r
.libPaths(c(file.path(getwd(), "Rlib"), .libPaths()))
```

### Single CRAN package
```r
install.packages("package_name", lib = file.path(getwd(), "Rlib"),
                 repos = "https://cloud.r-project.org")
```

### Single Bioconductor package
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = file.path(getwd(), "Rlib"))
BiocManager::install("package_name", lib = file.path(getwd(), "Rlib"),
                     update = FALSE, ask = FALSE)
```

### Batch install (CRAN)
```r
for (pkg in c("pkg1", "pkg2", "pkg3")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, lib = file.path(getwd(), "Rlib"),
                     repos = "https://cloud.r-project.org")
}
```

### GitHub package
```r
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools", lib = file.path(getwd(), "Rlib"))
devtools::install_github("author/package", lib = file.path(getwd(), "Rlib"))
```

## Troubleshooting

### Problem: Package installation fails with compilation error

**Cause**: Missing system-level dependencies.

**Solution**: Install the required system libraries:

```bash
# macOS (Homebrew)
brew install hdf5 libxml2 openssl curl

# Ubuntu/Debian
sudo apt-get install libcurl4-openssl-dev libxml2-dev libssl-dev libhdf5-dev

# CentOS/RHEL
sudo yum install libcurl-devel libxml2-devel openssl-devel hdf5-devel
```

### Problem: BiocManager::install() fails

**Cause**: Bioconductor version incompatible with R version.

**Solution**:
```r
BiocManager::install(version = "3.18")  # Match your R version
```

### Problem: Package loads from system library instead of Rlib

**Cause**: `.libPaths()` order issue.

**Solution**: Ensure Rlib is first:
```r
.libPaths(c(rlib_path, .libPaths()))
print(.libPaths())  # Rlib should be [1]
```

### Problem: Disk space

**Cause**: Bioconductor packages can be large (especially annotation databases).

**Solution**: Check sizes with `du -sh Rlib/` and remove unused packages:
```r
remove.packages("unused_package", lib = file.path(getwd(), "Rlib"))
```

## Migrating Rlib

### To another machine (same OS and R version)

```bash
# On source machine
tar -czf Rlib_backup.tar.gz Rlib/

# On target machine
tar -xzf Rlib_backup.tar.gz
```

### To a different OS or R version

Reinstall from scratch: run the package installation code from Step 2 above.

## Server Usage

### HPC / Cluster Environment

1. Make sure R is loaded (e.g., `module load R/4.3.0`)
2. Set the project directory as working directory
3. Run setup:

```bash
cd /path/to/project
Rscript -e 'dir.create("Rlib"); .libPaths(c("Rlib", .libPaths())); cat("Rlib ready\n")'
# Then run the package installation from Step 2 in R
```

### Docker

Add to your Dockerfile:

```dockerfile
COPY . /app/
WORKDIR /app
RUN Rscript -e 'dir.create("Rlib"); .libPaths(c("Rlib", .libPaths())); \
  install.packages(c("tidyverse", "survival", "glmnet"), lib = "Rlib"); \
  BiocManager::install(c("DESeq2", "clusterProfiler", "limma"), lib = "Rlib")'
```

### Shared Server

Each user should have their own project copy with their own Rlib:

```bash
cp -r /shared/universal-bioinformatics-skill ~/my-analysis
cd ~/my-analysis
# Run the Rlib setup and package installation from Steps 1-2 in R
```

## Package List

### CRAN Packages

| Package | Purpose |
|---------|---------|
| tidyverse | Data wrangling and visualization |
| data.table | Fast data manipulation |
| ggplot2 | Grammar of graphics plotting |
| pheatmap | Heatmap visualization |
| ComplexHeatmap | Advanced heatmaps (Bioconductor) |
| RColorBrewer | Color palettes |
| ggrepel | Non-overlapping text labels |
| patchwork | Combine ggplots |
| survival | Survival analysis |
| survminer | Survival visualization |
| glmnet | LASSO/elastic net regression |
| randomForest | Random forest models |
| e1071 | SVM and other ML algorithms |
| caret | ML model training framework |
| pROC | ROC curve analysis |
| timeROC | Time-dependent ROC |
| rms | Regression modeling strategies |
| ggpubr | Publication-ready plots |
| corrplot | Correlation visualization |
| igraph | Network analysis |
| ggraph | Network visualization |
| WGCNA | Weighted gene co-expression |
| FactoMineR | Multivariate analysis |
| factoextra | Multivariate visualization |
| VennDiagram | Venn diagrams |
| UpSetR | UpSet plots |

### Bioconductor Packages

| Package | Purpose |
|---------|---------|
| DESeq2 | Differential expression (count data) |
| edgeR | Differential expression (count data) |
| limma | Differential expression (microarray/RNA-seq) |
| clusterProfiler | GO/KEGG enrichment |
| org.Hs.eg.db | Human gene annotation |
| org.Mm.eg.db | Mouse gene annotation |
| enrichplot | Enrichment visualization |
| ReactomePA | Reactome pathway analysis |
| GSVA | Gene set variation analysis |
| GSEABase | GSEA infrastructure |
| msigdbr | MSigDB gene sets |
| AnnotationDbi | Annotation database interface |
| STRINGdb | STRING PPI network |
| maftools | Mutation analysis |
| TCGAbiolinks | TCGA data access |
| SummarizedExperiment | Data container |
| SingleCellExperiment | Single-cell data container |
| scater | Single-cell QC and visualization |
| scran | Single-cell normalization |
| DOSE | Disease ontology enrichment |
