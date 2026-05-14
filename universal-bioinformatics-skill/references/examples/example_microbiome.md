# Case study: Gut microbiome and metabolic disease (T2DM)

## 1. Research questions

**Main question:** How do gut microbial community structure and predicted functions associate with metabolic disease (e.g. type 2 diabetes)?

**Sub-questions:**
1. How does microbiome composition differ between healthy controls and T2DM?
2. Which genera/species correlate with glycemic and lipid traits?
3. What functional shifts does PICRUSt2 / Tax4Fun2 predict?
4. How do fecal/serum metabolites change?
5. Can integrated microbiome–metabolite analysis highlight host–microbe axes?

## 2. Data and design

### Cohort

| Group | Criteria | n | Biospecimens |
|-------|----------|---|--------------|
| Control | Healthy adults, fasting glucose < 5.6 mmol/L | 50 | Stool + serum |
| T2DM | ADA criteria | 50 | Stool + serum |

### Omics layers

| Layer | Platform | Output |
|-------|----------|--------|
| 16S rRNA | Illumina MiSeq, V3–V4 | OTU/ASV table |
| Untargeted metabolomics | LC-MS/MS (+/− modes) | Feature intensity matrix |
| Clinical labs | Routine biochemistry | FBG, HbA1c, TG, TC, HDL, LDL, … |

## 3. Analysis workflow

### Step 1: 16S processing (QIIME2)

```bash
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left-f 0 --p-trim-left-r 0 \
  --p-trunc-len-f 280 --p-trunc-len-r 240 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --p-n-threads 8

qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza
```

### Step 2: Alpha diversity

```r
library(phyloseq)
library(microbiome)
library(ggpubr)

ps <- qza_to_phyloseq(
  features = "table.qza",
  tree = "rooted-tree.qza",
  taxonomy = "taxonomy.qza",
  metadata = "metadata.tsv"
)

ps_rare <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 42)

alpha_div <- estimate_richness(ps_rare, measures = c("Observed", "Shannon", "Simpson"))
alpha_div$group <- sample_data(ps_rare)$group

wilcox.test(Shannon ~ group, data = alpha_div)

ggboxplot(alpha_div, x = "group", y = "Shannon",
          color = "group", add = "jitter",
          palette = c("#00AFBB", "#FC4E07")) +
  stat_compare_means(method = "wilcox.test")
```

### Step 3: Beta diversity

```r
library(vegan)

dist_bray <- phyloseq::distance(ps_rare, method = "bray")
dist_unifrac <- phyloseq::distance(ps_rare, method = "unifrac")
dist_wunifrac <- phyloseq::distance(ps_rare, method = "wunifrac")

ordination <- ordinate(ps_rare, method = "PCoA", distance = dist_bray)
plot_ordination(ps_rare, ordination, color = "group") +
  stat_ellipse(level = 0.95) +
  theme_bw()

adonis2(dist_bray ~ group, data = as(sample_data(ps_rare), "data.frame"),
        permutations = 999)

anosim(dist_bray, sample_data(ps_rare)$group)

bd <- betadisper(dist_bray, sample_data(ps_rare)$group)
permutest(bd, pairwise = TRUE)
```

### Step 4: Differential abundance

```r
# LEfSe export + run locally or Galaxy

library(DESeq2)
dds <- phyloseq_to_deseq2(ps, ~ group)
dds <- DESeq(dds, test = "Wald", fitType = "parametric")
res <- results(dds, contrast = c("group", "T2DM", "Control"))
res_sig <- subset(res, padj < 0.05)

library(ANCOMBC)
ancom_res <- ancombc2(
  data = ps,
  fix_formula = "group",
  p_adj_method = "BH",
  alpha = 0.05,
  tax_level = "Genus"
)

library(Maaslin2)
maaslin_res <- Maaslin2(
  input_data = as.data.frame(otu_table(ps)),
  input_metadata = as.data.frame(sample_data(ps)),
  output = "maaslin2_output",
  fixed_effects = c("group"),
  normalization = "TSS",
  transform = "LOG"
)
```

### Step 5: Functional prediction

```r
# PICRUSt2 CLI example:
# picrust2_pipeline.py -s rep-seqs.fasta -i table.biom -o picrust2_output -p 8

library(Tax4Fun2)
# Map 16S to KEGG KO / pathway abundances

# Differential testing on pathway matrices mirrors taxon-level models
```

### Step 6: Metabolomics

```r
library(mixOmics)

metab <- read.csv("metabolite_abundance.csv", row.names = 1)

pca_res <- pca(metab, ncomp = 5)
plotIndiv(pca_res, group = metadata$group, legend = TRUE)

plsda_res <- plsda(metab, metadata$group, ncomp = 3)
plotIndiv(plsda_res, group = metadata$group, legend = TRUE,
          ellipse = TRUE, title = "PLS-DA")

vip_scores <- vip(plsda_res)
key_metabolites <- names(which(vip_scores[, 1] > 1.0))

diff_metab <- apply(metab, 2, function(x) {
  wilcox.test(x ~ metadata$group)$p.value
})
diff_metab_adj <- p.adjust(diff_metab, method = "BH")
sig_metabolites <- names(which(diff_metab_adj < 0.05))
```

### Step 7: Microbe–metabolite associations

```r
library(psych)

cor_res <- corr.test(
  diff_genus_abundance,
  diff_metabolite_abundance,
  method = "spearman",
  adjust = "BH"
)

library(pheatmap)
pheatmap(
  cor_res$r,
  display_numbers = matrix(ifelse(cor_res$p < 0.05, "*", ""),
                           nrow = nrow(cor_res$r)),
  color = colorRampPalette(c("blue", "white", "red"))(100),
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

library(igraph)
cor_matrix <- cor_res$r
cor_matrix[cor_res$p >= 0.05] <- 0

g <- graph_from_adjacency_matrix(
  abs(cor_matrix) > 0.5,
  mode = "undirected",
  diag = FALSE
)
```

### Step 8: Link to clinical variables

```r
envfit_res <- envfit(ordination$vectors, clinical_data, permutations = 999)

mantel(dist_bray, dist(clinical_data$HbA1c), method = "spearman", permutations = 999)

dbrda_res <- dbrda(dist_bray ~ HbA1c + BMI + FBG + TG, data = clinical_data)
anova(dbrda_res, by = "terms")
```

## 4. Figure plan (12 main panels)

| ID | Figure | Type | Purpose |
|----|--------|------|---------|
| Fig 1A | Alpha diversity | Box | Richness |
| Fig 1B | PCoA (Bray–Curtis) | Scatter | Beta diversity |
| Fig 2A | Phylum stacked bar | Stacked bar | Community structure |
| Fig 2B | Top 20 genera heatmap | Heatmap | Abundance patterns |
| Fig 3A | LEfSe bar | Horizontal bar | Biomarker taxa |
| Fig 3B | DA volcano | Volcano | Effect sizes |
| Fig 4A | PLS-DA scores | Scatter | Metabolome separation |
| Fig 4B | VIP top 20 | Bar | Key metabolites |
| Fig 5A | Correlation heatmap | Heatmap | Microbe–metabolite |
| Fig 5B | Association network | Network | Multi-omic links |
| Fig 6A | dbRDA triplot | Triplot | Clinical drivers |
| Fig 6B | Genus vs HbA1c | Scatter | Example association |

### Supplementary

| ID | Content |
|----|---------|
| Fig S1 | Rarefaction curves |
| Fig S2 | PERMANOVA tables |
| Fig S3 | PICRUSt2 pathway contrasts |
| Fig S4 | Representative MS2 spectra |
| Fig S5 | Mediation / SEM schematics |

## 5. Multi-layer association framework

```
Layer 1: Microbial structure
├── Alpha diversity vs metabolic traits
├── Enterotype / community types
└── Core network modules

Layer 2: Microbe ↔ metabolite
├── Spearman correlation matrices
├── CCA / RDA joint ordination
├── MOFA across omics
└── Procrustes concordance

Layer 3: Causal inference (when justified)
├── Mediation (microbe → metabolite → clinical outcome)
├── SEM
└── Mendelian randomization (if GWAS instruments exist)
```

### Mediation sketch

```r
library(mediation)

model_m <- lm(butyrate ~ Bacteroides + age + sex + BMI, data = combined_data)
model_y <- lm(FBG ~ Bacteroides + butyrate + age + sex + BMI, data = combined_data)

med_res <- mediate(model_m, model_y,
                   treat = "Bacteroides", mediator = "butyrate",
                   boot = TRUE, sims = 1000)
summary(med_res)
```

## 6. Practical notes

### 16S considerations

1. **Compositional data:** Use CLR transforms or ANCOM-BC rather than raw relative abundances alone.
2. **Rarefaction:** Required for some alpha metrics; optional for robust beta distances.
3. **Reference DB:** Silva vs Greengenes vs GTDB can change taxonomy.
4. **PCR bias:** Primer pair affects detected taxa.

### Metabolomics considerations

1. **Batch:** Track QC samples for drift correction.
2. **ID confidence:** Distinguish Level 1 (standards) vs Level 2 (database) annotations.
3. **Ion modes:** Integrate positive/negative carefully to avoid duplicate features.
4. **Validation:** Targeted panels strengthen untargeted hits.
