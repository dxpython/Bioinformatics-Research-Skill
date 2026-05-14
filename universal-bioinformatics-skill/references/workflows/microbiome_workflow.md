# Microbiome Analysis Workflow

## 1. Applicable Scenarios

- Gut, oral, skin, or environmental microbiome profiling
- Comparing microbial communities between conditions
- Identifying differentially abundant taxa
- Functional prediction from 16S data
- Metagenomics analysis (species-level resolution)
- Microbiome-metabolome-host integration

## 2. Input Data Requirements

| Input | Description | Format |
|-------|-------------|--------|
| OTU/ASV table | Taxa (rows) × Samples (columns) abundance | CSV/TSV/BIOM |
| Taxonomy table | OTU/ASV to taxonomy mapping (Kingdom to Species) | CSV/TSV |
| Sample metadata | Sample IDs, group labels, covariates | CSV/TSV |
| Phylogenetic tree (optional) | For UniFrac distances | Newick format |
| Representative sequences (optional) | For functional prediction | FASTA |

## 3. Data Validation Checklist

- [ ] OTU/ASV IDs match between abundance table and taxonomy
- [ ] Sample IDs match between abundance table and metadata
- [ ] Sequencing depth is sufficient (rarefaction curve plateaus)
- [ ] Remove mitochondrial and chloroplast sequences
- [ ] Filter low-abundance OTUs (e.g., present in <2 samples or <0.01% total)

## 4. Recommended Analysis Steps

### Step 1: Data Preprocessing
1. Remove singletons and rare OTUs/ASVs
2. Check sequencing depth across samples
3. Rarefaction curve analysis
4. Optional: rarefy to minimum depth (controversial but standard)

### Step 2: Alpha Diversity

| Index | Measures | Interpretation |
|-------|----------|----------------|
| Observed species | Richness | Number of species present |
| Chao1 | Richness (with rare species estimate) | Estimated total richness |
| Shannon | Richness + Evenness | Higher = more diverse and even |
| Simpson | Dominance | Higher = more even distribution |
| Faith's PD | Phylogenetic diversity | Incorporates evolutionary distances |

- Compare between groups: Wilcoxon/Kruskal-Wallis test
- Visualization: boxplots with statistical comparisons

### Step 3: Beta Diversity

| Distance | Type | Features |
|----------|------|----------|
| Bray-Curtis | Abundance-based | Most common; considers abundance |
| Jaccard | Presence/absence | Ignores abundance |
| Weighted UniFrac | Phylogenetic + abundance | Incorporates phylogenetic tree |
| Unweighted UniFrac | Phylogenetic + presence | Phylogenetic tree required |

- Ordination: PCoA (Principal Coordinates Analysis), NMDS
- Statistical testing: PERMANOVA (adonis2), ANOSIM
- Visualization: PCoA/NMDS plots colored by group

### Step 4: Community Composition
- Stacked barplot at Phylum/Family/Genus level
- Relative abundance comparison
- Top N taxa across groups

### Step 5: Differential Abundance Analysis

| Tool | Method | Features |
|------|--------|----------|
| LEfSe | LDA Effect Size | Classic; identifies biomarker taxa |
| DESeq2 | Negative binomial | Borrowed from RNA-seq |
| ANCOM-BC | Bias-corrected | Handles compositionality |
| MaAsLin2 | General linear model | Handles covariates well |
| ALDEx2 | CLR transformation | Compositional data approach |

### Step 6: Functional Prediction (16S)
- PICRUSt2: Predict metagenome from 16S
- Tax4Fun2: KEGG pathway prediction
- BugBase: Predict phenotypes (pathogenicity, biofilm, etc.)

### Step 7: Metagenomics (if available)
- Taxonomic profiling: MetaPhlAn, Kraken2
- Functional profiling: HUMAnN3
- Genome assembly: MEGAHIT, metaSPAdes
- MAG recovery: MetaBAT2, CONCOCT

### Step 8: Microbiome-Metabolome Integration
- Correlate taxa abundance with metabolite levels
- Spearman correlation with FDR correction
- Network visualization of significant correlations
- Procrustes analysis of ordination results
- Mantel test for distance matrix correlation

## 5. Recommended Statistical Methods

| Analysis | Recommended Method | Alternative |
|----------|-------------------|-------------|
| Alpha diversity | Wilcoxon / Kruskal-Wallis | Linear model with covariates |
| Beta diversity | PERMANOVA (adonis2, vegan) | ANOSIM, MRPP |
| Dispersion | PERMDISP2 (betadisper) | — (required before PERMANOVA) |
| Differential abundance | ANCOM-BC2, MaAsLin2 | DESeq2, ALDEx2, LEfSe |
| Multiple testing | Benjamini-Hochberg FDR | Storey's q-value |
| Correlation (microbe-metabolite) | Spearman correlation | SparCC, SPIEC-EASI (compositional-aware) |
| Functional prediction | PICRUSt2 + Welch's t-test | Tax4Fun2, Piphillin |
| Ordination | PCoA (Bray-Curtis, UniFrac) | NMDS, DPCoA |

- Compositional data requires special handling: use centered log-ratio (CLR) transformation or compositional-aware methods
- PERMANOVA assumes equal group dispersions—always test with PERMDISP2 first
- For differential abundance: ANCOM-BC2 controls FDR and handles compositionality better than LEfSe
- Report effect sizes (e.g., R² from PERMANOVA, LDA scores from LEfSe) alongside p-values

## 6. Recommended Figures

| Figure | Description |
|--------|-------------|
| Rarefaction curves | Sequencing depth adequacy |
| Alpha diversity boxplots | Shannon, Chao1 per group |
| PCoA/NMDS plot | Beta diversity ordination |
| Stacked barplot | Community composition at phylum/genus |
| LEfSe cladogram | Differentially abundant taxa |
| LEfSe barplot | LDA scores |
| Heatmap | Top differential taxa |
| Correlation network | Microbiome-metabolome |

## 7. Result Interpretation Template

> 16S rRNA amplicon sequencing of [N] samples yielded [M] ASVs after quality filtering. Alpha diversity analysis revealed significantly higher Shannon diversity in [Group1] compared to [Group2] (p = X.XX, Wilcoxon test) (Figure XA). PCoA based on Bray-Curtis distance showed significant separation between groups (PERMANOVA R² = X.XX, p = 0.001) (Figure XB). At the phylum level, [Group1] showed higher relative abundance of [Phylum] (Figure XC). LEfSe analysis identified [N] differentially abundant taxa, with [Genus1] and [Genus2] enriched in [Group1] (LDA > 3.0, p < 0.05) (Figure XD).

## 8. Manuscript Writing Template

### Methods

> **16S rRNA gene sequencing and microbiome analysis**
>
> [Fecal / mucosal / environmental] samples were collected from [N] subjects and DNA was extracted using the [kit name] according to the manufacturer's protocol. The V[3-4] hypervariable region of the 16S rRNA gene was amplified using primers [primer sequences] and sequenced on an [Illumina MiSeq / NovaSeq] platform (2 × 250 bp paired-end). Raw sequencing reads were processed using QIIME2 (vX.X) with the DADA2 plugin for amplicon sequence variant (ASV) inference, quality filtering, and chimera removal. Taxonomy was assigned using the SILVA (v138) / Greengenes2 reference database. ASVs classified as mitochondria and chloroplasts were removed. Samples were rarefied to [X] reads for diversity analyses. Alpha diversity indices (Shannon, Chao1, and observed ASVs) were calculated and compared between groups using the [Wilcoxon rank-sum test / Kruskal-Wallis test]. Beta diversity was assessed using Bray-Curtis, weighted UniFrac, and unweighted UniFrac distances, visualized by principal coordinate analysis (PCoA), and statistically tested by PERMANOVA (adonis2, 999 permutations, vegan package). Differential abundance analysis was performed using [ANCOM-BC2 / MaAsLin2 (vX.X) / LEfSe], with [age, sex, BMI, and medication] included as covariates. Functional potential was predicted using PICRUSt2 with the MetaCyc database. Microbiome-metabolome associations were assessed by Spearman's rank correlation with BH adjustment.

## 9. Common Issues and Risks

| Issue | Solution |
|-------|----------|
| Uneven sequencing depth | Rarefy or use rarefaction-free methods (DESeq2) |
| Compositional data bias | Use CLR transformation (ALDEx2, ANCOM-BC) |
| Low taxonomic resolution (16S) | Use ASVs instead of OTUs; consider metagenomics |
| Batch effects | Include in statistical model; use same protocols |
| Contamination | Include negative controls; use decontam package |
| Multiple testing | Apply FDR correction for all comparisons |

## 10. Experimental Validation Suggestions

1. **qPCR**: Validate abundance of key taxa with species-specific primers
2. **Cultivation**: Isolate and culture key species
3. **Germ-free mouse models**: Test causal relationships
4. **Fecal microbiota transplant (FMT)**: Validate functional effects
5. **Metabolite measurements**: Validate predicted metabolic functions
