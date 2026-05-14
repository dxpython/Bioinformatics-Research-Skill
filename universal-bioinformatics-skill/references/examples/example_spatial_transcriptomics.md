# Case study: Spatial transcriptomics of tumor microenvironment heterogeneity

## 1. Research questions

**Main question:** How does spatial organization of the tumor microenvironment (TME) relate to immune infiltration and progression?

**Sub-questions:**
1. How do transcriptomes differ across tumor core, invasive margin, and stroma?
2. How are immune cells distributed in space—immune excluded vs inflamed?
3. Are there spatial “hotspots” of tumor–immune interaction?
4. How can Visium be integrated with scRNA-seq?
5. Do spatial co-expression modules reveal new functional zones?

## 2. Data and design

### Experiment

| Item | Detail |
|------|--------|
| Platform | 10x Visium (FFPE or fresh-frozen) |
| Tissue | Tumor section spanning core, margin, adjacent normal |
| Replicates | 4–6 sections with biological replication |
| Resolution | ~55 µm spot; roughly 1–10 cells per spot |
| Imaging | Matched H&E |

### File layout (Space Ranger)

```
spaceranger_output/
├── filtered_feature_bc_matrix/
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
├── spatial/
│   ├── tissue_positions_list.csv
│   ├── scalefactors_json.json
│   ├── tissue_hires_image.png
│   └── tissue_lowres_image.png
└── analysis/
```

## 3. Analysis workflow

### Step 1: Load and QC

```python
import scanpy as sc
import squidpy as sq
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

adata = sc.read_visium("spaceranger_output/")

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.spatial(adata, color='n_genes_by_counts', ax=axes[0], title='Genes per spot')
sc.pl.spatial(adata, color='total_counts', ax=axes[1], title='UMI per spot')
sc.pl.spatial(adata, color='pct_counts_mt', ax=axes[2], title='MT %')

adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.pct_counts_mt < 25, :]

print(f"Spots retained: {adata.n_obs}")
```

### Step 2: Normalize, HVGs, cluster

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor='seurat_v3')
adata_hvg = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=50)
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.5)

adata.obs['spatial_cluster'] = adata_hvg.obs['leiden']
sc.pl.spatial(adata, color='spatial_cluster', spot_size=1.2)
```

### Step 3: Region annotation (markers + scores)

```python
tumor_markers = ['EPCAM', 'KRT18', 'MKI67', 'TOP2A']
immune_markers = ['PTPRC', 'CD3D', 'CD68', 'CD19']
stroma_markers = ['COL1A1', 'ACTA2', 'FAP', 'VIM']

sc.pl.spatial(adata, color=tumor_markers + immune_markers + stroma_markers,
              ncols=4, spot_size=1.2)

sc.tl.score_genes(adata, tumor_markers, score_name='tumor_score')
sc.tl.score_genes(adata, immune_markers, score_name='immune_score')
sc.tl.score_genes(adata, stroma_markers, score_name='stroma_score')

scores = adata.obs[['tumor_score', 'immune_score', 'stroma_score']]
adata.obs['region'] = scores.idxmax(axis=1).str.replace('_score', '')
```

### Step 4: Differential expression and SVGs

```python
sc.tl.rank_genes_groups(adata, groupby='region', method='wilcoxon')
sc.pl.rank_genes_groups_dotplot(adata, n_genes=5)

sq.gr.spatial_neighbors(adata, coord_type='generic', n_neighs=6)
sq.gr.spatial_autocorr(adata, mode='moran', genes=adata.var_names[:2000])

svg_df = adata.uns['moranI'].sort_values('I', ascending=False)
top_svg = svg_df.head(20).index.tolist()

sc.pl.spatial(adata, color=top_svg[:6], ncols=3, spot_size=1.2)
```

### Step 5: Spatial communication and co-occurrence

```python
sq.gr.ligrec(
    adata,
    n_perms=1000,
    cluster_key='region',
    use_raw=True
)

sq.pl.ligrec(adata, cluster_key='region',
             source_groups='tumor', target_groups='immune')

sq.gr.nhood_enrichment(adata, cluster_key='spatial_cluster')
sq.pl.nhood_enrichment(adata, cluster_key='spatial_cluster')
```

### Step 6: Spatial modules (Squidpy / SpatialDE)

```python
sq.gr.co_occurrence(
    adata,
    cluster_key='region',
    spatial_key='spatial'
)
sq.pl.co_occurrence(adata, cluster_key='region')

import SpatialDE
counts = pd.DataFrame(adata.X.toarray(), columns=adata.var_names)
coords = adata.obsm['spatial']
sample_info = pd.DataFrame({'x': coords[:, 0], 'y': coords[:, 1]})

results = SpatialDE.run(sample_info, counts)
significant_genes = results[results['qval'] < 0.05]

patterns = SpatialDE.spatial_patterns(
    sample_info, counts,
    significant_genes,
    n_patterns=5, length=1.0
)
```

## 4. Joint analysis with scRNA-seq

### Strategy A: Deconvolution (cell2location)

```python
import cell2location
from cell2location.models import RegressionModel

sc_adata = sc.read_h5ad("scrnaseq_reference.h5ad")

RegressionModel.setup_anndata(sc_adata, labels_key='cell_type')
ref_model = RegressionModel(sc_adata)
ref_model.train(max_epochs=250)

inf_aver = ref_model.export_posterior()

cell2location.models.Cell2location.setup_anndata(adata)
mod = cell2location.models.Cell2location(
    adata, cell_state_df=inf_aver,
    N_cells_per_location=10,
    detection_alpha=20
)
mod.train(max_epochs=30000)

adata = mod.export_posterior(adata, sample_kwargs={'num_samples': 1000})

cell_types = ['T_CD8', 'Macrophage', 'Fibroblast', 'B_cell']
sc.pl.spatial(adata, color=[f'q05_cell_abundance_w_sf_{ct}' for ct in cell_types],
              ncols=2, spot_size=1.2)
```

### Strategy B: Tangram label transfer

```python
import tangram as tg

sc_adata = sc.read_h5ad("scrnaseq_reference.h5ad")
tg.pp_adatas(sc_adata, adata, genes=marker_genes)

ad_map = tg.map_cells_to_space(
    sc_adata, adata,
    mode='cells',
    density_prior='rna_count_based',
    num_epochs=500
)

tg.project_cell_annotations(ad_map, adata, annotation='cell_type')
```

### Strategy C: Multimodal joint embedding

Use frameworks such as MultiVI / MIRA when multi-modal training data are available.

## 5. Figure plan (12 main panels)

| ID | Figure | Type | Purpose |
|----|--------|------|---------|
| Fig 1A | H&E overview | Histology | Anatomy reference |
| Fig 1B | Spatial QC maps | Spatial heatmap | QC on tissue |
| Fig 2A | Spatial clusters on H&E | Spatial | Domains |
| Fig 2B | UMAP of spots | UMAP | Transcriptional space |
| Fig 3A | Marker genes in space | Spatial feature | Zoning |
| Fig 3B | Region DE dotplot | Dot | Molecular identity |
| Fig 4A | Top SVGs | Spatial heatmap | Spatial patterns |
| Fig 4B | SpatialDE patterns | Pattern plot | Metagenes |
| Fig 5A | Deconvolved cell types | Spatial | Composition |
| Fig 5B | Immune gradient | Gradient / summary | Infiltration |
| Fig 6A | Neighborhood enrichment | Heatmap | Co-localization |
| Fig 6B | Spatial LR dotplot | Dot | Communication |

### Supplementary

| ID | Content |
|----|---------|
| Fig S1 | Per-section spot counts |
| Fig S2 | Cross-section cluster concordance |
| Fig S3 | cell2location training curves |
| Fig S4 | Extended SVG gallery |

## 6. Notes

### Visium limitations

1. **Resolution:** Each spot can contain multiple cells—use deconvolution.
2. **Sensitivity:** Whole transcriptome but lower per-gene sensitivity than scRNA-seq.
3. **2D sections:** May not capture full 3D anatomy.
4. **FFPE vs fresh:** FFPE RNA quality is lower; adjust expectations and methods.

### Spatial technology comparison

| Technology | Resolution | Genes | Throughput | Typical use |
|-------------|------------|-------|------------|-------------|
| 10x Visium | 55 µm | Whole transcriptome | Medium | Standard spatial RNA |
| Visium HD | ~2 µm | Whole transcriptome | Medium | Higher resolution |
| MERFISH | Subcellular | 100–1000 targeted | High | Targeted panels |
| Slide-seq | ~10 µm | Whole transcriptome | Medium | High-res whole transcriptome |
| Stereo-seq | ~500 nm | Whole transcriptome | High | Ultra-high resolution |

### Choosing a joint-analysis strategy

| Strategy | Pros | Cons | When to use |
|----------|------|------|---------------|
| cell2location | Probabilistic, quantitative | Heavy compute | Precise cell fractions |
| Tangram | Flexible gene projection | Marker dependent | Gene-level mapping |
| RCTD | Fast, Visium-oriented | Strong assumptions | Quick deconvolution |
| SPOTlight | Lightweight NMF-based | Lower accuracy | Exploratory screens |
