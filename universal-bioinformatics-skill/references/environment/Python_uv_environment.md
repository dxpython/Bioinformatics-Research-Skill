# Python Environment Guide: uv Package Manager

## What is uv and Why Use It?

**uv** is an extremely fast Python package and project manager written in Rust. It replaces pip, pip-tools, virtualenv, and pyenv in a single tool with 10-100x faster performance.

### Why uv instead of pip + requirements.txt?

| Feature | pip + requirements.txt | uv + pyproject.toml |
|---------|----------------------|---------------------|
| Speed | Slow (resolves dependencies sequentially) | Extremely fast (parallel Rust resolver) |
| Lock file | No native lock file | `uv.lock` for reproducible builds |
| Virtual env | Separate tool (venv/virtualenv) | Built-in (`uv venv`) |
| Dependency resolution | Fragile, often conflicts | Robust resolver, clear error messages |
| Project standard | Ad-hoc | PEP 621 compliant (pyproject.toml) |
| Reproducibility | Poor (no lock, no hash verification) | Excellent (lock file with hashes) |

## Quick Start

### Step 1: Install uv

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# macOS (Homebrew)
brew install uv

# Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# Verify installation
uv --version
```

### Step 2: Initialize the Project

```bash
cd universal-bioinformatics-skill

# Create virtual environment (uses Python specified in pyproject.toml)
uv venv

# Sync all dependencies from pyproject.toml
uv sync
```

### Step 3: Run Scripts

```bash
# Run a Python analysis script
uv run python your_analysis.py

# Start Jupyter Lab
uv run jupyter lab

# Run any Python command
uv run python -c "import scanpy; print(scanpy.__version__)"
```

## Common Operations

### Add a New Package

```bash
# Add a package (updates pyproject.toml and uv.lock)
uv add scanpy

# Add a development dependency
uv add --dev pytest

# Add with version constraint
uv add "pandas>=2.0,<3.0"

# Add an optional dependency group
uv add --optional gpu rapids-singlecell
```

### Remove a Package

```bash
uv remove package-name
```

### Update Packages

```bash
# Update all packages
uv sync --upgrade

# Update a specific package
uv add package-name --upgrade
```

### Export Dependencies

```bash
# Export to requirements.txt (for compatibility)
uv export --format requirements-txt > requirements.txt

# Export with hashes for security
uv export --format requirements-txt --hash > requirements.txt
```

### Check Environment

```bash
# List installed packages
uv pip list

# Show package info
uv pip show scanpy

# Check for dependency conflicts
uv sync --dry-run
```

## Reproducing the Environment on Another Machine

### Method 1: Using uv (Recommended)

```bash
# On the new machine
git clone <project-repo>
cd universal-bioinformatics-skill
uv sync  # Installs exact versions from uv.lock
```

### Method 2: Using requirements.txt (Fallback)

```bash
# Export on source machine
uv export --format requirements-txt > requirements.txt

# On target machine
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Working with Jupyter

```bash
# Start Jupyter Lab
uv run jupyter lab

# Register kernel for use in external Jupyter
uv run python -m ipykernel install --user --name bioinformatics --display-name "Bioinformatics (uv)"
```

## Working Alongside R

This project uses both R (with Rlib) and Python (with uv). They are independent:

- R scripts: Run with `Rscript` or `source()` in R, using the Rlib cache
- Python scripts: Run with `uv run python`, using the uv virtual environment
- Data exchange: Use CSV/TSV files or h5ad format for interoperability

### Typical Mixed Workflow

```bash
# Step 1: Run R differential analysis
Rscript analysis/deg_analysis.R

# Step 2: Use Python for ML on the R results
uv run python analysis/ml_biomarker.py

# Step 3: Visualize with R
Rscript analysis/enrichment_plots.R
```

## Server / HPC Usage

### With module system

```bash
module load python/3.11
curl -LsSf https://astral.sh/uv/install.sh | sh
cd /path/to/project
uv sync
```

### In Docker

```dockerfile
FROM python:3.11-slim
COPY --from=ghcr.io/astral-sh/uv:latest /uv /usr/local/bin/uv
WORKDIR /app
COPY pyproject.toml uv.lock ./
RUN uv sync --frozen
COPY . .
CMD ["uv", "run", "python", "your_script.py"]
```

### In CI/CD

```yaml
# GitHub Actions example
- uses: astral-sh/setup-uv@v4
- run: uv sync
- run: uv run pytest
```

## Troubleshooting

### Problem: Python version not found

```bash
# Install the required Python version
uv python install 3.11
```

### Problem: Package build fails (needs system deps)

```bash
# macOS
brew install hdf5 libomp

# Ubuntu
sudo apt-get install libhdf5-dev build-essential
```

### Problem: CUDA/GPU packages

```bash
# Install with GPU support
uv add --optional gpu rapids-singlecell
uv sync --extra gpu
```

### Problem: Slow on first run

This is normal—uv builds a cache. Subsequent runs will be much faster.
