# Case study: Transcriptomic mechanisms of high-fat-diet–induced hepatic lipid dysregulation in mice

## 1. Research questions

**Main question:** How does a high-fat diet (HFD) reshape transcriptional networks to drive hepatic lipid dysregulation in mice?

**Sub-questions:**
1. What are the major transcriptional changes in liver after HFD?
2. Which lipid-metabolism pathways are activated or repressed?
3. How do target programs of key TFs (e.g. SREBP1, PPARγ) shift?
4. Do differential genes overlap known NAFLD/NASH disease genes?
5. Can we nominate therapeutic intervention targets?

## 2. Data sources and experimental design

### Design

| Group | Diet | Duration | n | Tissue |
|-------|------|----------|---|--------|
| Control | Chow (10% kcal from fat) | 16 weeks | 4 | Liver |
| HFD | High-fat diet (60% kcal from fat) | 16 weeks | 4 | Liver |

### Sequencing

- Platform: Illumina NovaSeq 6000
- Library: poly(A) mRNA-seq
- Read length: PE150
- Depth: ≥30M clean reads per sample

### Preprocessing (shell sketch)

```bash
# Step 1: QC
fastp -i sample_R1.fq.gz -I sample_R2.fq.gz \
      -o sample_clean_R1.fq.gz -O sample_clean_R2.fq.gz \
      -h sample_fastp.html -j sample_fastp.json

# Step 2: Align to mouse mm10/GRCm38
STAR --runThreadN 8 \
     --genomeDir /ref/mm10_STAR_index \
     --readFilesIn sample_clean_R1.fq.gz sample_clean_R2.fq.gz \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix sample_

# Step 3: Quantify
featureCounts -a /ref/gencode.vM25.annotation.gtf \
              -o counts.txt \
              -T 8 -p --countReadPairs \
              sample1.bam sample2.bam ...
```

## 3. Analysis workflow

### Step 1: QC and exploration

```r
library(DESeq2)

count_data <- read.table("counts.txt", header = TRUE, row.names = 1)
col_data <- data.frame(
  condition = factor(c(rep("Control", 4), rep("HFD", 4))),
  row.names = colnames(count_data)
)

dds <- DESeqDataSetFromMatrix(count_data, col_data, design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 4
dds <- dds[keep, ]

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = "condition")

library(pheatmap)
sample_cor <- cor(assay(vsd))
pheatmap(sample_cor, annotation_col = col_data)
```

### Step 2: Differential expression

```r
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "HFD", "Control"))
res <- res[order(res$padj), ]

deg <- subset(res, abs(log2FoldChange) > 1 & padj < 0.05)
up_genes <- rownames(subset(deg, log2FoldChange > 0))
down_genes <- rownames(subset(deg, log2FoldChange < 0))

cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")
```

### Step 3: GO/KEGG enrichment

```r
library(clusterProfiler)
library(org.Mm.eg.db)

ego <- enrichGO(
  gene = up_genes,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

gene_ids <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

ekegg <- enrichKEGG(
  gene = gene_ids$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

# Run separate enrichments for up/down and compare
```

### Step 4: GSEA

```r
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)
gene_list <- sort(gene_list, decreasing = TRUE)

library(msigdbr)
hallmark_mm <- msigdbr(species = "Mus musculus", category = "H")

gsea_res <- GSEA(
  gene_list,
  TERM2GENE = hallmark_mm[, c("gs_name", "gene_symbol")],
  pvalueCutoff = 0.05
)

# Lipid-related Hallmark sets of interest:
# HALLMARK_FATTY_ACID_METABOLISM
# HALLMARK_ADIPOGENESIS
# HALLMARK_CHOLESTEROL_HOMEOSTASIS
```

### Step 5: Lipid pathway focus

```r
lipid_pathways <- list(
  fatty_acid_synthesis = c("Fasn", "Acaca", "Acacb", "Scd1", "Elovl6"),
  fatty_acid_oxidation = c("Cpt1a", "Cpt2", "Acadm", "Acadl", "Hadha"),
  cholesterol_synthesis = c("Hmgcr", "Hmgcs1", "Sqle", "Fdft1", "Cyp51"),
  triglyceride_synthesis = c("Dgat1", "Dgat2", "Gpam", "Agpat2", "Lpin1"),
  lipid_droplet = c("Plin2", "Plin5", "Cidea", "Cidec", "Fitm2"),
  transcription_factors = c("Srebf1", "Srebf2", "Pparg", "Ppara", "Nr1h3")
)

lipid_deg <- res[unlist(lipid_pathways), ]

library(ComplexHeatmap)
lipid_expr <- assay(vsd)[unlist(lipid_pathways), ]
# ... grouped heatmap
```

### Step 6: TF regulatory context

```r
# Option A: SCENIC (resource-intensive)

# Option B: curated TF targets
srebp1_targets <- c("Fasn", "Scd1", "Acaca", "Gpam", "Elovl6")
ppara_targets <- c("Cpt1a", "Acox1", "Ehhadh", "Acadm", "Cyp4a10")

# Query ChEA3 / DoRothEA for TF–target enrichment in DEGs

# Option C: WGCNA
library(WGCNA)
# Build signed network, relate modules to HFD phenotype, extract hub genes
```

### Step 7: Overlap with disease gene sets

```r
nafld_genes <- c("PNPLA3", "TM6SF2", "MBOAT7", "HSD17B13", ...)

overlap <- intersect(toupper(rownames(deg)), nafld_genes)

fisher.test(matrix(c(
  length(overlap),
  length(setdiff(rownames(deg), nafld_genes)),
  length(setdiff(nafld_genes, rownames(deg))),
  total_genes - length(union(rownames(deg), nafld_genes))
), nrow = 2))
```

## 4. Figure plan (10 main panels)

| ID | Name | Type | Purpose |
|----|------|------|---------|
| Fig 1A | PCA | Scatter | QC / separation |
| Fig 1B | Sample correlation heatmap | Heatmap | Consistency |
| Fig 2A | Volcano | Volcano | Global DE |
| Fig 2B | Top 30 DEG heatmap | Heatmap | Key genes |
| Fig 3A | GO-BP bubble | Bubble | Processes |
| Fig 3B | KEGG bar | Bar | Pathways |
| Fig 4A | GSEA ridge (lipid) | Ridge | Directional enrichment |
| Fig 4B | GSEA running score | Enrichment plot | Single pathway |
| Fig 5 | Lipid gene heatmap | Heatmap | Lipid focus |
| Fig 6 | TF–target network | Network | Mechanism cartoon |

### Supplementary

| ID | Content |
|----|---------|
| Fig S1 | fastp QC summary |
| Fig S2 | Expression density after normalization |
| Fig S3 | WGCNA module–trait relationships |
| Fig S4 | Extended GO dotplot |

## 5. Mechanism narrative

### Model

```
HFD (16 weeks)
      │
      ▼
Hepatic lipid overload
      │
      ├──► SREBP1c hyperactivation
      │         │
      │         ▼
      │    De novo lipogenesis ↑↑
      │    (Fasn↑, Scd1↑, Acaca↑)
      │
      ├──► PPARα downregulation
      │         │
      │         ▼
      │    β-oxidation ↓↓
      │    (Cpt1a↓, Acadm↓)
      │
      ├──► Inflammatory programs
      │    (NF-κB, TNF-α, IL-6)
      │         │
      │         ▼
      │    Macrophage infiltration → NASH progression
      │
      └──► ER stress
           (Xbp1s↑, Atf6↑, Ddit3↑)
                │
                ▼
           Steatosis → steatohepatitis
```

### Discussion paragraph template

> We profiled livers from HFD-fed mice and identified [X] up- and [Y] downregulated genes. Enrichment highlighted dysregulated lipid metabolism, especially de novo lipogenesis and fatty-acid oxidation.
>
> SREBP1c targets (Fasn, Scd1, Acaca) rose while PPARα-driven β-oxidation genes (Cpt1a, Acadm) fell—a “lipogenic switch” linked to steatosis [citation].
>
> GSEA further indicated inflammatory and ER-stress programs, consistent with progression from simple steatosis toward steatohepatitis. [Gene] emerged as a hub in the co-expression network and merits functional follow-up.

## 6. Practical notes

### Mouse RNA-seq specifics

1. **Genome build:** Use GRCm38/mm10 or GRCm39/mm39 consistently with annotations.
2. **Gene IDs:** Mouse symbols differ from human (e.g. Fasn vs FASN).
3. **Ortholog mapping:** Use `biomaRt` when comparing to human cohorts.
4. **Replicates:** Aim for n ≥ 4 per arm in animal studies.
5. **Batch:** If sequenced in batches, include batch as a covariate.

### Common reviewer comments

| Comment | Response |
|---------|----------|
| “Low n” | Report power analysis |
| “Need protein validation” | Western blot / IHC for key genes |
| “Other tissues?” | Add adipose or serum if available |
| “Clinical translation” | Discuss validation in human NAFLD cohorts |
| “Single time point” | Acknowledge limitation; propose time course |
