# Spatial Transcriptomics Analysis Workflow

## 1. Applicable Scenarios

- Mapping gene expression in tissue spatial context
- Identifying spatial domains and tissue architecture
- Cell type deconvolution of spatial spots
- Spatial co-localization and neighborhood analysis
- Spatial cell-cell communication
- Comparing lesion regions vs. normal regions
- Integration with scRNA-seq for enhanced resolution

## 2. Input Data Requirements

| Input | Description | Format |
|-------|-------------|--------|
| Spatial expression data | Spot/cell × gene matrix with spatial coordinates | h5ad (with obsm['spatial']), 10X Visium (Space Ranger output) |
| Tissue image (optional) | H&E-stained tissue image | PNG/TIFF |
| scRNA-seq reference (optional) | For deconvolution | h5ad / Seurat object |

## 3. Data Validation Checklist

- [ ] Spatial coordinates exist in the data (obsm['spatial'] or tissue_positions)
- [ ] Expression matrix contains counts
- [ ] Tissue image aligns with spot positions
- [ ] Sufficient spots per tissue section (typically 1000–5000 for Visium)

## 4. Recommended Analysis Steps

### Step 1: QC for Spatial Spots
- Total counts per spot
- Number of genes per spot
- Mitochondrial gene percentage
- Spatial distribution of QC metrics (overlay on tissue image)
- Filter low-quality spots

### Step 2: Normalization and Feature Selection
- Normalize (total counts or SCTransform)
- Select highly variable genes
- Same as scRNA-seq but consider spatial variability

### Step 3: Spatial Clustering and Domain Identification
- Standard clustering: PCA → neighbors → Leiden/Louvain
- Spatially-aware clustering: BayesSpace, SpaGCN, STAGATE
- Visualize clusters on tissue coordinates

### Step 4: Spatial Marker Genes
- Identify marker genes for each spatial domain
- Moran's I statistic for spatially variable genes
- spatialDE for spatial differential expression

### Step 5: Spatial Feature Plots
- Plot gene expression on tissue coordinates
- Compare expression patterns across spatial domains
- Overlay on H&E image

### Step 6: Cell Type Deconvolution

| Tool | Method | Features |
|------|--------|----------|
| cell2location | Bayesian model | Gold standard for Visium |
| Tangram | Deep learning mapping | Maps scRNA-seq cells to spots |
| SPOTlight | NMF-based | Fast, lightweight |
| RCTD | Reference-based | Robust to noise |
| stereoscope | Probabilistic | Good for complex tissues |

**Requirement**: All deconvolution methods need a scRNA-seq reference.

### Step 7: Spatial Co-localization
- Neighborhood enrichment analysis (squidpy)
- Co-occurrence analysis: which cell types tend to co-localize?
- Ripley's statistics for spatial point patterns

### Step 8: Spatial Communication
- MISTy: multi-view spatial modeling
- SpaTalk: spatially resolved cell-cell communication
- COMMOT: optimal transport-based communication
- CellChat (adapted for spatial)

### Step 9: Region Comparison
- Define regions of interest (tumor core, tumor margin, normal tissue)
- Differential expression between regions
- Pathway activity differences across regions
- Cell type composition differences

### Step 10: Integration with scRNA-seq
- Map scRNA-seq cell types to spatial spots
- Validate spatial expression patterns with single-cell resolution
- Combine for comprehensive tissue characterization

## 5. Recommended Statistical Methods

| Analysis | Method |
|----------|--------|
| Spatially variable genes | Moran's I, spatialDE |
| Neighborhood enrichment | Permutation test (squidpy) |
| Deconvolution | Bayesian inference (cell2location) |
| Region comparison | Wilcoxon or pseudobulk DE |

## 6. Recommended Figures

| Figure | Description |
|--------|-------------|
| Spatial QC plots | nCounts, nGenes overlaid on tissue |
| Spatial clustering | Cluster assignments on tissue |
| H&E overlay | Clusters on histology image |
| Spatial feature plot | Gene expression heatmap on tissue |
| Deconvolution map | Cell type proportions per spot |
| Neighborhood enrichment | Cell type co-localization heatmap |
| Co-occurrence plot | Spatial co-occurrence curves |
| Region comparison | Volcano plot for region-specific genes |


## 7. Result Interpretation Template

> Spatial transcriptomics analysis of [N] tissue spots identified [K] distinct spatial domains (Figure XA). Domain [X] corresponded to [anatomical region], characterized by high expression of [genes] (Figure XB). Cell type deconvolution using cell2location revealed predominant [cell type] in domain [X], while domain [Y] was enriched for [cell type] (Figure XC). Neighborhood enrichment analysis showed significant co-localization between [cell type A] and [cell type B] (z-score = X.XX, p < 0.05), suggesting spatial interaction in [region] (Figure XD).

## 8. Manuscript Writing Template

### Methods
> Spatial transcriptomics data were generated using the 10x Genomics Visium platform. Data were processed using Scanpy (v1.X) and Squidpy (v1.X). Spots with fewer than [N] total counts were excluded. Cell type deconvolution was performed using cell2location with a scRNA-seq reference from [source]. Spatially variable genes were identified using Moran's I statistic. Neighborhood enrichment analysis was performed with [N] permutations.

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| Low spot resolution (Visium) | Use deconvolution; consider subcellular platforms |
| Tissue damage / holes | Mask damaged regions before analysis |
| Batch effects across sections | Use shared gene sets; batch correction |
| Deconvolution inaccuracy | Validate with IHC/IF; use multiple methods |
| Image alignment issues | Re-run Space Ranger with manual alignment |

## 10. Experimental Validation Suggestions

1. **IHC/IF**: Validate cell type localization predictions
2. **smFISH/RNAscope**: Validate spatial gene expression at single-molecule level
3. **Laser capture microdissection + RNA-seq**: Validate region-specific gene expression
4. **Serial sections**: Compare adjacent sections with different staining
