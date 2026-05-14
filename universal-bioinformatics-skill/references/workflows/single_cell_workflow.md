# Single-Cell RNA-seq Analysis Workflow

## 1. Applicable Scenarios

- Identifying cell types and subtypes in heterogeneous tissues
- Comparing cell composition between conditions (disease vs. normal)
- Discovering novel cell states or rare populations
- Cell-cell communication analysis
- Pseudotime trajectory and lineage inference
- Copy number variation inference from scRNA-seq (inferCNV)
- Gene regulatory network inference (SCENIC)
- Single-cell multiome analysis (joint ATAC + RNA) for chromatin accessibility and gene expression
- CITE-seq / ADT analysis (surface protein + RNA) for protein-level validation

## 2. Input Data Requirements

| Input | Description | Format |
|-------|-------------|--------|
| Expression data | UMI count matrix | h5ad, 10X (barcodes.tsv + features.tsv + matrix.mtx), RDS (Seurat), h5 |
| Multiome data | Paired ATAC + RNA per cell | h5ad (MuData), h5 (10X multiome), fragments.tsv.gz + Seurat RDS |
| CITE-seq/ADT data | RNA + antibody-derived tag counts | h5ad (with protein layer), Seurat RDS (with ADT assay), CITE-seq-Count output |
| Metadata (optional) | Sample origin, condition, batch | CSV/TSV or embedded in h5ad/Seurat |
| Reference markers (optional) | Known cell type markers | CSV/TSV or literature |

## 3. Data Validation Checklist

- [ ] File format is readable (h5ad for scanpy, 10X/RDS for Seurat)
- [ ] Matrix is UMI counts (non-negative integers)
- [ ] Number of cells and genes is reasonable (typically 1K–100K cells, 20K–30K genes)
- [ ] Metadata matches cell barcodes
- [ ] If multiple samples, batch information is available
- [ ] For multiome: fragments file is provided; ATAC peaks are called; paired barcodes match between modalities
- [ ] For CITE-seq: ADT count matrix is available; ADT antibodies are annotated with target protein names

## 4. Recommended Analysis Steps

### Step 1: Quality Control (QC)

**Key QC metrics:**
- `nFeature_RNA` / `n_genes_by_counts`: Number of genes detected per cell
- `nCount_RNA` / `total_counts`: Total UMI counts per cell
- `percent.mt` / `pct_counts_mt`: Percentage of mitochondrial gene counts

**Typical filtering criteria (adjustable by tissue type):**
- 200 < nFeature < 5000 (remove empty droplets and doublets)
- nCount > 500
- percent.mt < 20% (strict: <10%; lenient for certain tissues: <25%)

**Additional QC:**
- Doublet detection: DoubletFinder (Seurat) or Scrublet (scanpy)
- Ambient RNA removal: SoupX or CellBender

### Step 2: Normalization

| Method | Tool | When to Use |
|--------|------|-------------|
| LogNormalize | Seurat | Default for most analyses |
| SCTransform | Seurat | Better for integrating datasets |
| scran pooling | scran | Robust for variable library sizes |
| normalize_total + log1p | scanpy | Default scanpy pipeline |

### Step 3: Feature Selection (HVG)

- Select highly variable genes (HVGs): typically 2000–3000 genes
- Seurat: `FindVariableFeatures(method = "vst")`
- scanpy: `sc.pp.highly_variable_genes(n_top_genes=2000)`

### Step 4: Dimensionality Reduction

1. **PCA**: Compute 50 PCs, use elbow plot to determine significant PCs (typically 15–30)
2. **UMAP**: 2D embedding for visualization (not for quantitative analysis)
3. **t-SNE**: Alternative to UMAP (preserves local structure better, but slower)

### Step 5: Clustering

- Graph-based clustering: Louvain or Leiden algorithm
- Resolution parameter controls granularity (0.1–2.0; default ~0.5–1.0)
- Start with moderate resolution, adjust based on biological knowledge
- Seurat: `FindNeighbors()` + `FindClusters(resolution = 0.8)`
- scanpy: `sc.pp.neighbors()` + `sc.tl.leiden(resolution = 0.8)`

### Step 6: Marker Gene Identification

- Find marker genes per cluster: Wilcoxon rank-sum test (default), t-test, or MAST
- Seurat: `FindAllMarkers(min.pct = 0.25, logfc.threshold = 0.25)`
- scanpy: `sc.tl.rank_genes_groups(method='wilcoxon')`
- Filter: adjusted p < 0.05, log2FC > 0.5, pct.1 > 0.25

### Step 7: Cell Type Annotation

**Manual annotation (recommended):**
1. Examine top marker genes per cluster
2. Compare with known canonical markers
3. Use databases: CellMarker, PanglaoDB, The Human Protein Atlas

**Automated annotation:**
- SingleR (reference-based)
- scType
- CellTypist

**Visualization:**
- DotPlot: marker gene expression and percentage per cluster
- FeaturePlot: individual gene expression on UMAP
- VlnPlot: gene expression distribution per cluster
- Heatmap: top markers across all clusters

### Step 8: Cell Proportion Analysis

- Calculate cell type proportions per sample/condition
- Statistical comparison: chi-square test, Fisher's exact, or permutation test
- Visualization: stacked barplot, proportion comparison plot

### Step 9: Subclustering

- Re-cluster specific cell populations at higher resolution
- Identify subtypes (e.g., CD4+ T cell subtypes: Th1, Th2, Treg, Tfh)
- Repeat marker identification and annotation for subclusters

### Step 10: Pseudotime / Trajectory Analysis

| Tool | Method | Best For |
|------|--------|----------|
| Monocle3 | UMAP + graph | General trajectory inference |
| Slingshot | Minimum spanning tree | Branching lineages |
| scVelo | RNA velocity | Directionality of differentiation |
| PAGA | Graph abstraction (scanpy) | Connectivity between clusters |

### Step 11: Cell-Cell Communication

| Tool | Method | Features |
|------|--------|----------|
| CellChat | Ligand-receptor + signaling | Pathway-level, quantitative |
| CellPhoneDB | Statistical framework | Large L-R database |
| NicheNet | Ligand-target prediction | Links ligands to downstream targets |
| LIANA | Meta-framework | Combines multiple methods |

### Step 12: Gene Regulatory Networks (SCENIC)

- Infer transcription factor (TF) regulons
- Identify cell-type-specific regulatory programs
- Tools: pySCENIC or SCENIC+

### Step 13: Copy Number Inference (inferCNV)

- Infer large-scale chromosomal gains/losses from expression
- Distinguish malignant from non-malignant cells
- Requires a reference cell population (e.g., immune cells)

### Step 14: Single-Cell Multiome Analysis (ATAC + RNA)

**Joint analysis of paired chromatin accessibility (ATAC-seq) and gene expression (RNA-seq) from the same cell.**

**ATAC-specific QC:**
- TSS enrichment score (ENCODE standard: >7 for high-quality cells)
- Nucleosome banding pattern
- Fraction of reads in peaks (FRiP; >40% recommended)
- Fragment size distribution (should show nucleosome periodicity ~147 bp)

**ATAC processing pipeline:**
1. **Peak calling**: MACS2 per sample or pseudo-bulk; recall peaks across all cells
2. **Quantification**: Build gene activity matrix (counts in gene body + 2 kb upstream)
3. **Dimensionality reduction**: TF-IDF normalization + LSI (latent semantic indexing)
4. **Signac (R)** workflow: `FeatureMatrix()` → `RunTFIDF()` → `FindTopFeatures()` → `RunSVD()` → `RunUMAP(dims=2:30, reduction='lsi')`
5. **ArchR (R)** pipeline (scalable alternative): Arrow files → iterative LSI → Harmony integration

**Joint embedding (WNN in Seurat v5, MultiVI in scVI-tools):**
- Seurat: `FindMultiModalNeighbors()` with RNA + ATAC reductions → WNN graph → joint UMAP
- scVI-tools: MultiVI for joint probabilistic modeling of RNA + ATAC

**Linking peaks to genes:**
- Signac: `LinkPeaks()` — correlate ATAC peak accessibility with nearby gene expression
- Cicero (Monocle3): cis-co-accessibility networks
- Identify putative enhancers and distal regulatory elements

**Differential analysis:**
- Differential accessibility (DA) by cluster or condition: Wilcoxon, DESeq2 pseudobulk
- Overlap DA regions with DEGs to find coordinated regulation
- Motif enrichment at DA peaks: Signac `FindMotifs()`, HOMER, MEME-ChIP

**TF activity inference:**
- ChromVAR: infer transcription factor activity from motif accessibility deviation
- Signac: `RunChromVAR()` after `AddMotifs()`

**Key visualizations:**
- Joint UMAP split by modality (RNA, ATAC, WNN)
- Coverage plots at specific loci (e.g., gene body + cis-regulatory elements)
- Motif footprinting plots (TF binding patterns at binding sites)
- Volcano plots of DA peaks vs. DEGs side-by-side

**Caveats:**
- ATAC signal is sparser than RNA; interpret LSI components carefully (first component often captures technical depth)
- Gene activity scores are approximations, not direct expression measurements
- Peak-gene links are correlational, not confirmed regulatory interactions

### Step 15: CITE-seq / ADT Analysis (RNA + Surface Protein)

**Joint analysis of transcriptome and antibody-derived tag (ADT) counts from the same cell, enabling protein-level validation of cell identity and state.**

**ADT-specific QC:**
- ADT total counts per cell (remove cells with <100 ADT UMIs)
- ADT unique features (too few = failed antibody binding)
- Check for potential antibody aggregates (ultra-high ADT counts)
- Antibody isotype controls (if available) to assess non-specific binding

**ADT normalization:**
| Method | Tool | Description |
|--------|------|-------------|
| CLR (centered log-ratio) | Seurat `NormalizeData(method = "CLR")` | Default for compositional data |
| DSB (denoised and scaled by background) | dsb R package | Corrects for ambient and cell-intrinsic background |
| Log-normalize | Seurat | Simpler but ignores compositionality |

**Joint RNA + ADT clustering (WNN):**
- Seurat: `FindMultiModalNeighbors()` with RNA + ADT reductions → WNN graph → joint UMAP
- Protein-based markers can refine cluster separation when RNA signals are ambiguous (e.g., memory vs. naive T cells)

**Differential protein expression:**
- Wilcoxon test on CLR-normalized ADT values per cluster or condition
- Identify proteins that differentiate populations more clearly than their mRNA counterparts

**RNA-protein correlation:**
- Per cluster: correlate ADT signal with corresponding mRNA expression for each marker
- Discordant gene-protein pairs may indicate post-transcriptional regulation
- Visualization: scatter plot (ADT vs. RNA) per feature, biaxial feature plots

**Key applications:**
- **Cell annotation refinement**: Surface proteins (CD4, CD8, CD19, CD14, CD45RA/RO) provide definitive lineage markers
- **T cell differentiation**: CD45RA vs. CD45RO ADT distinguishes naive from memory cells better than RNA
- **Exhaustion markers**: Co-staining of PD-1, TIM-3, LAG-3 ADTs
- **Sample multiplexing**: HTO (hashtag oligo) demultiplexing via CITE-seq-Count or Seurat `HTODemux()`

**Caveats:**
- ADT measurements are limited to the antibody panel (typically 30–200 proteins)
- Some antibodies show non-specific binding at high concentrations
- ADT counts are compositional; CLR or DSB normalization is essential
- Batch effects between panels: include bridge samples or use anchor-based integration

## 5. Recommended Statistical Methods

| Analysis | Method |
|----------|--------|
| DE between clusters | Wilcoxon rank-sum test |
| DE between conditions | MAST, pseudobulk DESeq2 |
| Proportion comparison | Chi-square, Fisher's exact |
| Trajectory inference | Monocle3, scVelo |
| Communication significance | Permutation test (CellChat) |
| Multiome DA | Wilcoxon, pseudobulk DESeq2 on peak counts |
| Multiome TF activity | chromVAR deviation scores |
| CITE-seq ADT DE | Wilcoxon on CLR/DSB-normalized values |

**Important**: For between-condition comparisons, **pseudobulk analysis** (aggregating cells per sample, then using DESeq2/edgeR) is statistically more appropriate than cell-level tests.

## 6. Recommended Figures

| Figure | Description |
|--------|-------------|
| QC violin plots | nFeature, nCount, percent.mt before/after filtering |
| UMAP by cluster | Cell clusters colored by identity |
| UMAP by condition | Cells colored by sample/condition |
| DotPlot | Marker gene expression per cluster |
| FeaturePlot | Key gene expression on UMAP |
| Heatmap | Top 5-10 markers per cluster |
| Stacked barplot | Cell proportion per condition |
| Pseudotime trajectory | Cells colored by pseudotime |
| CellChat circle plot | Communication network strength |
| CellChat bubble plot | Ligand-receptor pairs |
| inferCNV heatmap | Chromosomal gains/losses |
| Multiome joint UMAP | RNA, ATAC, WNN embeddings side-by-side |
| Coverage plot | ATAC signal at specific loci |
| Motif footprinting | TF binding patterns |
| ADT feature plot | Surface protein expression on UMAP |
| RNA vs. ADT biaxial plot | Gene-protein correlation per cluster |
| ADT ridge/heatmap | Protein marker expression per cluster |


## 7. Result Interpretation Template

### Clustering Results
> Unsupervised clustering of [N] cells identified [K] distinct clusters (Figure XA). Based on canonical marker gene expression, these clusters were annotated as [list cell types] (Figure XB). [Cell type 1] was characterized by high expression of [Marker1] and [Marker2], while [Cell type 2] showed specific expression of [Marker3] (Figure XC).

### Proportion Changes
> Comparison of cell type proportions between [Condition1] and [Condition2] revealed a significant increase in [Cell type X] (XX% vs. YY%, p = X.XX) and a decrease in [Cell type Y] (XX% vs. YY%, p = X.XX) (Figure XD).

### Communication Results
> Cell-cell communication analysis identified [N] significant ligand-receptor interactions. [Cell type A] showed enhanced signaling to [Cell type B] via the [pathway] pathway in [Condition1] compared to [Condition2] (Figure XE).

## 8. Manuscript Writing Template

### Methods
> Single-cell RNA sequencing data were processed using Seurat (v5.X.X) / scanpy (v1.X.X). Cells with fewer than 200 or more than 5000 detected genes, or with mitochondrial gene percentage exceeding 20%, were excluded. After quality filtering, [N] cells remained for downstream analysis. Data were normalized using [method], and the top 2000 highly variable genes were selected. Principal component analysis was performed, and the first [N] PCs were used for UMAP visualization and graph-based clustering (Leiden algorithm, resolution = X.X). Cluster-specific marker genes were identified using the Wilcoxon rank-sum test (adjusted p < 0.05, log2FC > 0.5). Cell types were annotated based on known canonical markers.

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| Batch effects between samples | Use Harmony, BBKNN, or Seurat integration (CCA/RPCA) |
| Overclustering | Reduce resolution; merge similar clusters with biological justification |
| Underclustering | Increase resolution; subcluster specific populations |
| Doublets | Run DoubletFinder/Scrublet before analysis |
| Ambient RNA contamination | Run SoupX or CellBender |
| Rare cell types missed | Targeted subclustering; lower min.pct threshold |
| Annotation uncertainty | Use multiple references; validate with independent markers |
| Multiome: low TSS enrichment | Relax per-cell TSS cutoff; use pseudo-bulk aggregation |
| Multiome: LSI captures depth, not biology | Exclude dimension 1 from downstream analysis; verify with FRiP |
| Multiome: few peak-gene links | Use broader genomic window; increase cells per cluster for correlation |
| CITE-seq: ADT background noise | Apply DSB normalization; use isotype controls to set thresholds |
| CITE-seq: antibody cross-reactivity | Validate with independent protein methods (flow cytometry); flag ambiguous markers |

## 10. Experimental Validation Suggestions

1. **Flow cytometry / FACS**: Validate cell proportions
2. **Immunofluorescence (IF)**: Validate co-expression of marker genes at protein level
3. **IHC**: Validate spatial distribution in tissue sections
4. **Functional assays**: Validate cell state predictions (proliferation, activation)
5. **Spatial transcriptomics**: Validate spatial relationships between cell types
6. **CRISPR perturbation**: Validate predicted regulators from SCENIC
7. **ATAC-qPCR / reporter assay**: Validate predicted enhancer-promoter interactions from multiome analysis
8. **Flow cytometry**: Validate CITE-seq surface protein abundance differences; an independent antibody panel
