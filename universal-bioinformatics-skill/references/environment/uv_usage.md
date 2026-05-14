# uv Quick Reference Guide

## Installation

```bash
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# Homebrew
brew install uv

# Verify
uv --version
```

## Essential Commands

### Project Setup

```bash
# Create virtual environment
uv venv

# Install all dependencies from pyproject.toml
uv sync

# Install with optional dependency groups
uv sync --extra dev
uv sync --extra gpu
uv sync --all-extras
```

### Running Scripts

```bash
# Run a Python script
uv run python your_analysis.py

# Run Jupyter Lab
uv run jupyter lab

# Run Jupyter Notebook
uv run jupyter notebook

# Run IPython
uv run ipython

# Run any command in the venv
uv run python -c "import scanpy; print(scanpy.__version__)"
```

### Package Management

```bash
# Add a package
uv add scanpy

# Add dev dependency
uv add --dev pytest

# Add with version constraint
uv add "numpy>=1.24,<2.0"

# Remove a package
uv remove package-name

# Update all
uv sync --upgrade

# Update specific package
uv add package-name --upgrade
```

### Environment Info

```bash
# List installed packages
uv pip list

# Show package details
uv pip show scanpy

# Check Python version
uv run python --version

# Export requirements.txt
uv export --format requirements-txt > requirements.txt
```

### Python Version Management

```bash
# Install a specific Python version
uv python install 3.11

# List available versions
uv python list

# Pin Python version for project
uv python pin 3.11
```

## Bioinformatics-Specific Tips

### Running Analysis Scripts

```bash
# Single-cell analysis with scanpy
uv run python analysis/single_cell.py

# Spatial transcriptomics with squidpy
uv run python analysis/spatial.py

# Machine learning biomarker discovery
uv run python analysis/ml_biomarker.py

# Multi-omics network analysis
uv run python analysis/multiomics_network.py
```

### Interactive Analysis

```bash
# Start Jupyter Lab for interactive analysis
uv run jupyter lab

# Register kernel for external Jupyter installations
uv run python -m ipykernel install --user --name bioinfo
```

### Large Dataset Tips

```bash
# For large single-cell datasets, ensure sufficient memory
uv run python -c "
import scanpy as sc
sc.settings.n_jobs = 8  # Use 8 cores
adata = sc.read_h5ad('large_dataset.h5ad')
"
```

## Troubleshooting

| Issue | Solution |
|-------|----------|
| `uv: command not found` | Add `~/.cargo/bin` to PATH |
| Python version mismatch | `uv python install 3.11` |
| Build failures | Install system deps (hdf5, gcc) |
| Slow first install | Normal—cache is built; subsequent runs are fast |
| Permission denied | Don't use sudo with uv |
| Lock file conflict | `uv lock --upgrade` to regenerate |
