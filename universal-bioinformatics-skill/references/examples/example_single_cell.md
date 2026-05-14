# Case study: scRNA-seq of immune cell states in tumor tissue

## 1. Research questions

**Main question:** How do immune cell compositions and functional states change in the tumor microenvironment?

**Sub-questions:**
1. How do tumor vs adjacent-normal samples differ in cellular composition?
2. Do CD8+ T cells show an exhaustion program?
3. What is the polarization state of tumor-associated macrophages (TAMs)?
4. How does cell–cell communication differ between conditions?
5. What do differentiation trajectories look like for key lineages?

## 2. Data and experimental design

### Design

| Group | Source | Samples | Cells (approx.) |
|-------|--------|---------|-----------------|
| Tumor | Tumor core | 5 | ~5,000 per sample |
| Normal | Adjacent normal | 5 | ~5,000 per sample |

### Sequencing

- 10x Genomics Chromium (3′ v3.1)
- Illumina NovaSeq 6000
- Target ≥20,000 reads per cell

### Inputs

- Cell Ranger output: `filtered_feature_bc_matrix/` (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`)
- Or preprocessed `.h5` / `.h5ad`

## 3. Analysis workflow

### Step 1: Load and QC

```python
import scanpy as sc
import numpy as np
import pandas as pd

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white')

adata_list = []
for sample in ["Tumor_1", "Tumor_2", ..., "Normal_5"]:
    adata_tmp = sc.read_10x_h5(f"data/{sample}/filtered_feature_bc_matrix.h5")
    adata_tmp.obs['sample'] = sample
    adata_tmp.obs['condition'] = 'Tumor' if 'Tumor' in sample else 'Normal'
    adata_tmp.var_names_make_unique()
    adata_list.append(adata_tmp)

adata = sc.concat(adata_list, join='outer')

adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

adata = adata[adata.obs.n_genes_by_counts > 200, :]
adata = adata[adata.obs.n_genes_by_counts < 5000, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

print(f"Cells after filtering: {adata.n_obs}")
```

### Step 2: Normalize and integrate

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata.copy()

sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key='sample')
adata = adata[:, adata.var.highly_variable]

sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)

import harmonypy
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)
sc.external.pp.harmony_integrate(adata, key='sample')

sc.pp.neighbors(adata, use_rep='X_pca_harmony', n_neighbors=30)
```

### Step 3: Embedding and clustering

```python
sc.tl.umap(adata, min_dist=0.3)

for res in [0.3, 0.5, 0.8, 1.0, 1.5]:
    sc.tl.leiden(adata, resolution=res, key_added=f'leiden_{res}')

sc.pl.umap(adata, color=['leiden_0.8', 'condition', 'sample'],
           ncols=3, frameon=False)
```

### Step 4: Cell type annotation

```python
marker_genes = {
    'T_CD8': ['CD8A', 'CD8B', 'GZMB', 'PRF1'],
    'T_CD4': ['CD4', 'IL7R', 'CCR7', 'LEF1'],
    'Treg': ['FOXP3', 'IL2RA', 'CTLA4', 'IKZF2'],
    'NK': ['NKG7', 'GNLY', 'KLRD1', 'NCR1'],
    'B_cell': ['CD19', 'MS4A1', 'CD79A', 'PAX5'],
    'Plasma': ['JCHAIN', 'MZB1', 'SDC1', 'XBP1'],
    'Monocyte': ['CD14', 'LYZ', 'S100A8', 'S100A9'],
    'Macrophage': ['CD68', 'CD163', 'MSR1', 'MARCO'],
    'DC': ['CLEC9A', 'XCR1', 'CD1C', 'FCER1A'],
    'Mast': ['KIT', 'TPSAB1', 'TPSB2', 'CPA3'],
    'Epithelial': ['EPCAM', 'KRT18', 'KRT19', 'CDH1'],
    'Fibroblast': ['COL1A1', 'DCN', 'THY1', 'ACTA2'],
    'Endothelial': ['PECAM1', 'VWF', 'CDH5', 'FLT1'],
}

sc.pl.dotplot(adata, marker_genes, groupby='leiden_0.8',
              standard_scale='var')

cluster_annotation = {
    '0': 'T_CD8',
    '1': 'Macrophage',
    '2': 'T_CD4',
    # ... map each cluster
}
adata.obs['cell_type'] = adata.obs['leiden_0.8'].map(cluster_annotation)

# Optional: CellTypist / SingleR automation
import celltypist
model = celltypist.models.Model.load(model='Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model=model, majority_voting=True)
adata.obs['celltypist'] = predictions.predicted_labels.majority_voting
```

### Step 5: Differential expression (Tumor vs Normal within types)

```python
cell_types = adata.obs['cell_type'].unique()
deg_results = {}

for ct in cell_types:
    adata_ct = adata[adata.obs['cell_type'] == ct].copy()
    if adata_ct.obs['condition'].nunique() < 2:
        continue
    sc.tl.rank_genes_groups(adata_ct, groupby='condition',
                            groups=['Tumor'], reference='Normal',
                            method='wilcoxon')
    deg_results[ct] = sc.get.rank_genes_groups_df(adata_ct, group='Tumor')

cd8_deg = deg_results['T_CD8']
cd8_deg_sig = cd8_deg[(cd8_deg['pvals_adj'] < 0.05) &
                       (abs(cd8_deg['logfoldchanges']) > 0.5)]
```

### Step 6: Functional scoring (exhaustion, M1/M2)

```python
exhaustion_genes = ['PDCD1', 'CTLA4', 'LAG3', 'HAVCR2', 'TIGIT',
                    'TOX', 'ENTPD1', 'LAYN']
cytotoxicity_genes = ['GZMB', 'PRF1', 'GNLY', 'NKG7', 'IFNG',
                      'GZMA', 'GZMK', 'FASLG']

adata_cd8 = adata[adata.obs['cell_type'] == 'T_CD8'].copy()
sc.tl.score_genes(adata_cd8, exhaustion_genes, score_name='exhaustion_score')
sc.tl.score_genes(adata_cd8, cytotoxicity_genes, score_name='cytotoxicity_score')

adata_mac = adata[adata.obs['cell_type'] == 'Macrophage'].copy()
m1_genes = ['NOS2', 'IL1B', 'TNF', 'IL6', 'CD80', 'CD86']
m2_genes = ['CD163', 'MRC1', 'MSR1', 'ARG1', 'TGFB1', 'IL10']

sc.tl.score_genes(adata_mac, m1_genes, score_name='M1_score')
sc.tl.score_genes(adata_mac, m2_genes, score_name='M2_score')
```

### Step 7: CellChat (export to R)

```python
adata.write('adata_for_cellchat.h5ad')
```

```r
library(CellChat)
library(SeuratDisk)

cellchat_tumor <- createCellChat(object = seurat_tumor, group.by = "cell_type")
cellchat_tumor@DB <- CellChatDB.human
cellchat_tumor <- subsetData(cellchat_tumor)
cellchat_tumor <- identifyOverExpressedGenes(cellchat_tumor)
cellchat_tumor <- identifyOverExpressedInteractions(cellchat_tumor)
cellchat_tumor <- computeCommunProb(cellchat_tumor)
cellchat_tumor <- computeCommunProbPathway(cellchat_tumor)
cellchat_tumor <- aggregateNet(cellchat_tumor)

# Repeat for Normal, then compare
cellchat_merge <- mergeCellChat(list(Normal = cellchat_normal,
                                     Tumor = cellchat_tumor),
                                add.names = c("Normal", "Tumor"))

netVisual_diffInteraction(cellchat_merge)
rankNet(cellchat_merge, mode = "comparison")
```

### Step 8: Trajectory (Monocle3 / scVelo)

```python
import monocle3
# Example: CD8+ naive → effector → exhausted ordering
# Build cds, preprocess, reduce_dim, cluster_cells, learn_graph, order_cells

import scvelo as scv
scv.pp.filter_and_normalize(adata_cd8)
scv.pp.moments(adata_cd8)
scv.tl.velocity(adata_cd8)
scv.tl.velocity_graph(adata_cd8)
scv.pl.velocity_embedding_stream(adata_cd8, basis='umap')
```

## 4. Figure plan (14 main panels)

| ID | Figure | Type | Purpose |
|----|--------|------|---------|
| Fig 1A | QC violins | Violin | QC metrics |
| Fig 1B | Cells before/after filter | Bar | QC summary |
| Fig 2A | UMAP by cell type | UMAP | Atlas |
| Fig 2B | UMAP split Tumor/Normal | UMAP | Condition contrast |
| Fig 2C | Cell-type proportions | Stacked bar | Composition |
| Fig 3A | Marker dotplot | Dot | Annotation QC |
| Fig 3B | Feature plots | Feature | Spatial gene signal |
| Fig 4A | CD8 exhaustion score | Violin / box | State |
| Fig 4B | Macrophage M1 vs M2 | Scatter | Polarization |
| Fig 5A | Volcano (per type) | Volcano | DE |
| Fig 5B | DE heatmap | Heatmap | Programs |
| Fig 6A | CellChat network | Network | Communication |
| Fig 6B | LR bubble | Bubble | Interactions |
| Fig 7 | Trajectory | Trajectory | Ordering |

### Supplementary

| ID | Content |
|----|---------|
| Fig S1 | UMAP before/after Harmony |
| Fig S2 | Multi-resolution clustering |
| Fig S3 | CellTypist vs manual labels |
| Fig S4 | RNA velocity stream |
| Fig S5 | Pathway scores per cell type |

## 5. QC thresholds and integration tips

### QC guidelines

| Metric | Suggested | Notes |
|--------|-----------|-------|
| nGenes | 200–5000 | Too low: empty droplets; too high: doublets |
| nUMI | > 500 | Very low UMI cells are often poor quality |
| MT% | < 20% | High MT% suggests lysis / stress |
| Doublets | Scrublet < 0.25 | Remove before clustering |

### Integration methods

| Method | Strength | Use case |
|--------|----------|----------|
| Harmony | Fast, strong default | Routine multi-sample merge |
| Scanorama | Preserves biology | When batch vs biology is subtle |
| BBKNN | Graph-based | Mild batch effects |
| scVI | Deep model | Large atlases |

### Common issues

1. **Doublets:** Remove with Scrublet or DoubletFinder before clustering.
2. **Annotation:** Combine automated labels, manual markers, and literature.
3. **Rare states:** Lower resolution can absorb rare clusters—subcluster if needed.
4. **Pseudoreplication:** For cohort-level DE, consider pseudobulk models.
5. **CellChat:** Distinguish statistical communication scores from true secretion proof.
