# Epigenomics Analysis Workflow

## 1. Applicable Scenarios

1. **DNA methylation (450K/850K array or WGBS/RRBS)**: Identify differentially methylated positions (DMPs) and regions (DMRs).
2. **ATAC-seq**: Chromatin accessibility, regulatory regions, and transcription factor activity.
3. **ChIP-seq**: Histone marks (H3K4me3, H3K27ac, H3K27me3, etc.) or transcription factor binding.
4. **CUT&Tag / CUT&RUN**: Low-input epigenomic profiling.
5. **Epigenome–transcriptome integration**: Link methylation/chromatin state to gene expression.
6. **Cancer epigenetic biomarkers**: Methylation-based biomarker discovery.

## 2. Input Data Requirements

### DNA Methylation (Array)
| Data type | Format | Description |
|-----------|--------|-------------|
| Beta matrix | CSV / RDS | Rows = CpG probes (cg IDs), columns = samples, values 0–1 |
| Raw IDAT | .idat | Illumina scanner raw signal files |
| Sample metadata | CSV | sample_id, group, age, sex, tissue, Sentrix_ID, etc. |

### ATAC-seq / ChIP-seq
| Data type | Format | Description |
|-----------|--------|-------------|
| BAM | .bam | Aligned reads, sorted + indexed |
| Peaks | BED / narrowPeak | Peak calling output |
| Fragment file (ATAC) | .tsv.gz | Fragment coordinates per read pair |
| BigWig | .bw | Signal tracks for visualization |
| Sample metadata | CSV | sample_id, group, IP_target (ChIP) |

**Key requirements**:
- Methylation array: ≥3 samples per group; ≥5 recommended.
- ATAC-seq: ≥2 biological replicates per group; ≥3 recommended.
- ChIP-seq: Input / IgG controls; ≥2 replicates per group.
- Consistent genome build (hg19/hg38, mm9/mm10).
- ATAC-seq: check fragment size distribution (nucleosome periodicity expected).

## 3. Data Validation Checklist

### 3.1 DNA Methylation Array QC
```
1. Signal quality:
   - Fraction of probes with detection p-value < 0.01 (sample-level > 95%)
   - Samples with many failed probes may be low quality

2. Probe filtering:
   - Remove probes with detection p > 0.01 in > 5% of samples
   - Remove SNP-associated probes (Zhou et al. annotation)
   - Remove cross-reactive probes
   - Remove sex-chromosome probes if sex differences are not of interest

3. Normalization:
   - Functional normalization (funnorm) recommended
   - Or BMIQ / SWAN (Infinium I/II type bias)

4. Batch effects:
   - Inspect technical covariates (Sentrix_ID, Sentrix_Position)
   - ComBat / SVA correction

5. QC plots:
   - Beta density plots
   - MDS plot (sample clustering)
   - Sex prediction check
```

### 3.2 ATAC-seq QC
```
1. Mapping:
   - Mapping rate > 80%
   - Remove mitochondrial reads (typically < 20% is acceptable)
   - Remove PCR duplicates

2. Fragment size distribution:
   - Nucleosome periodicity (~200 bp spacing)
   - Sub-nucleosomal fragments (< 150 bp) reflect open chromatin

3. Library complexity:
   - NRF (Non-Redundant Fraction) > 0.7
   - PBC1 (PCR Bottleneck Coefficient 1) > 0.7

4. FRiP (Fraction of Reads in Peaks):
   - FRiP > 0.2 (prefer > 0.3)

5. TSS enrichment:
   - TSS enrichment score > 6 indicates high quality
   - Reflects enrichment at transcription start sites
```

### 3.3 ChIP-seq QC
```
1. Mapping rate > 80%
2. Remove duplicates and low-quality reads (MAPQ >= 10)
3. FRiP > 0.01 (broad marks) or > 0.1 (narrow marks)
4. NSC (Normalized Strand Coefficient) > 1.05
5. RSC (Relative Strand Coefficient) > 0.8
6. IDR (Irreproducible Discovery Rate) for replicate concordance
7. Fingerprint plot to assess ChIP enrichment
```

## 4. Recommended Analysis Steps

### Step 1: DNA Methylation

#### 1A: DMP analysis
```
- Tools: minfi + limma (recommended), ChAMP, methylKit
- Use M-values for statistics; Beta for biological interpretation
- Model: M ~ group + age + sex + batch (SVA surrogate variables)
- Thresholds: |ΔBeta| > 0.1 (or 0.2) AND adjusted p < 0.05
- Output: ΔBeta, p-value, adj.p-value, gene annotation per CpG
```

#### 1B: DMR analysis
```
- Tools: DMRcate (recommended), bumphunter, dmrseq, comb-p
- DMRs are often more robust than single DMPs (spatial correlation of CpGs)
- DMRcate: lambda = 1000, C = 2
- Filters: ≥3 CpGs per DMR, mean |ΔBeta| > 0.1
- Annotate DMRs to genes, promoters, enhancers, etc.
```

#### 1C: Functional annotation of methylation
```
- CpG context: CpG island, shore, shelf, open sea
- Gene-centric: TSS200, TSS1500, 5'UTR, 1stExon, Body, 3'UTR, Intergenic
- Link promoter methylation to gene silencing
- Enhancer methylation (with H3K4me1/H3K27ac when available)
```

### Step 2: ATAC-seq

#### 2A: Peak calling
```
- Tool: MACS2 (recommended)
  macs2 callpeak -t sample.bam -f BAMPE --nomodel --shift -100 --extsize 200
  --keep-dup all -g hs -q 0.05 --broad (for broad peaks)
- Or Genrich (ATAC-aware alternative)
- Consensus peaks: intersect or merge across samples/replicates (DiffBind/bedtools)
```

#### 2B: Differential accessibility
```
- Tools: DiffBind (recommended), DESeq2, csaw
- Consensus peak set → read counts → differential analysis
- DiffBind: dba() → dba.count() → dba.contrast() → dba.analyze()
- Thresholds: |log2FC| >= 1.0, FDR < 0.05
```

#### 2C: Peak annotation and functional analysis
```
- ChIPseeker: annotatePeak() for genes and genomic features
- GREAT: region-based functional enrichment
- Homer: GO/pathway enrichment for peaks
```

### Step 3: Motif enrichment
```
- Homer findMotifsGenome.pl: de novo + known motifs
- MEME-ChIP / MEME suite
- Input: sequences of differential peaks (FASTA)
- Background: genome or non-differential peaks
- Identify enriched TF binding motifs
- Homer: fast for large data; MEME: sensitive for smaller sets
```

### Step 4: TF activity inference
```
- chromVAR (ATAC-seq): motif accessibility-based TF activity
  - Input: fragment count matrix + motif DB (JASPAR/CIS-BP)
  - Output: TF deviation score (z-score)
  - Differential TF activity

- TOBIAS (ATAC-seq): TF footprinting
  - Finer bound vs unbound states

- diffTF (ChIP-seq/ATAC-seq): compare TF activity between conditions
```

### Step 5: Epigenome–transcriptome integration
```
5.1 Methylation vs expression
    - Promoter methylation vs expression: expect negative correlation
    - Gene body methylation vs expression: may be positive
    - Scatter plots + Spearman correlation

5.2 Accessibility vs expression
    - Promoter accessibility vs expression
    - Enhancer accessibility vs target genes
    - GREAT or activity-by-contact models for enhancer–target links

5.3 Multi-layer epigenomics
    - ChromHMM / Segway: chromatin states
    - Combine histone marks to define functional states
    - e.g. Active promoter = H3K4me3 + H3K27ac
          Poised enhancer = H3K4me1 + H3K27me3
```

## 5. Recommended Statistical Methods

| Goal | Method | Tool |
|------|---------|------|
| DMP | limma on M-values | minfi + limma |
| DMR | Gaussian kernel smoothing | DMRcate |
| ATAC/ChIP diff peaks | Negative binomial | DiffBind (DESeq2 mode) |
| Peak calling | Model-based | MACS2 |
| Peak annotation | Nearest gene + features | ChIPseeker |
| Motif enrichment | De novo + known | Homer, MEME |
| TF activity | Deviation / footprinting | chromVAR, TOBIAS |
| Methylation–expression | Spearman | R (Hmisc) |
| Chromatin state | HMM segmentation | ChromHMM |
| Batch correction | SVA / ComBat | sva |
| CpG feature enrichment | Fisher's exact | missMethyl |
| Pathway (methylation) | gometh (probe-count bias) | missMethyl |

## 6. Recommended Figures

| Figure | Purpose | Tool |
|--------|---------|------|
| **Beta density** | Global methylation distribution | ggplot2, minfi |
| **Volcano** | DMPs / differential peaks | EnhancedVolcano |
| **Heatmap** | DMP/DMR/peak patterns | ComplexHeatmap |
| **Manhattan** | Genome-wide DMP view | qqman, CMplot |
| **Annotation pie** | Genomic distribution of peaks/DMPs | ChIPseeker |
| **TSS enrichment** | ATAC QC | deepTools, ggplot2 |
| **Fragment size** | Nucleosome pattern (ATAC) | ggplot2 |
| **Signal heatmap** | Peak-centered signals | deepTools computeMatrix + plotHeatmap |
| **Profile plot** | TSS/TES signal | deepTools |
| **Motif logo** | Enriched motifs | Homer, ggseqlogo |
| **TF deviation heatmap** | chromVAR TF scores | ComplexHeatmap |
| **Scatter** | Methylation vs expression | ggplot2 |
| **Genome browser** | Locus-specific tracks | IGV, Gviz, trackViewer |
| **Circos** | Genome-wide methylation/accessibility | circlize |
| **ChromHMM emissions** | State definitions | ChromHMM |

---


## 7. Result Interpretation Template

### DMP results
```
Using Illumina [450K/EPIC] methylation arrays, after QC and functional normalization,
[N] CpG sites were analyzed. With |ΔBeta| > [0.1/0.2] and adjusted p < 0.05,
[X] differentially methylated positions (DMPs) were identified, including
[hyper count] hypermethylated and [hypo count] hypomethylated sites. DMPs were enriched in
[genomic regions, e.g. promoter, gene body, intergenic]; [X%] fell in CpG islands and shores.
```

### ATAC-seq results
```
ATAC-seq identified [N] consensus peaks, of which [X] showed differential accessibility
(|log2FC| >= 1.0, FDR < 0.05). Differential peaks were enriched in [promoter/distal intergenic/intronic]
regions. Motif enrichment showed that [condition]-specific open chromatin was enriched for
[TF] motifs (p = [X], enrichment = [Y]-fold), suggesting increased [TF] activity in [disease/condition].
chromVAR confirmed elevated deviation scores for [TF] in [group] (p = [Z]).
```

### Integration results
```
Integrating methylation and transcriptomics, promoter methylation of [X] genes was negatively
correlated with expression (Spearman r < -0.5, p < 0.05); genes such as [gene list] showed
hypermethylation with downregulated expression, consistent with epigenetic silencing.
```

## 8. Manuscript Writing Template

### Methods
```
DNA methylation analysis:
DNA methylation was profiled using Illumina Infinium MethylationEPIC BeadChip. Raw IDAT
files were processed with minfi (v[version]). Probes with detection p-value > 0.01 in > 5%
of samples, SNP-associated probes, cross-reactive probes, and sex-chromosome probes were
removed. Data were normalized with functional normalization. Batch effects were corrected
with SVA. DMPs were identified with limma on M-values, adjusting for [covariates]. DMPs with
|ΔBeta| > [threshold] and adjusted p < 0.05 were significant. DMRs were called with DMRcate
(lambda=1000, C=2). GO enrichment used gometh (missMethyl) to account for probe-number bias.

ATAC-seq analysis:
ATAC-seq reads were aligned to [hg38/mm10] with Bowtie2. Mitochondrial reads, PCR duplicates,
and MAPQ < 10 reads were removed. Peaks were called with MACS2 [parameters]. Differential
accessibility used DiffBind (DESeq2 mode). Motif enrichment used Homer findMotifsGenome.pl.
TF activity was inferred with chromVAR and the JASPAR 2022 motif database.
```

## 9. Common Issues and Risks

| Issue | Risk | Mitigation |
|-------|------|------------|
| Array batch effects | High | SVA/ComBat; balance chip layout at design |
| Cell-composition confounding (blood) | High | EpiDISH / RefFreeEWAS cell proportions |
| Too many DMPs | Medium | Raise |ΔBeta| to 0.2; focus on DMRs |
| High mtDNA in ATAC | High | Optimize library prep; filter mt reads |
| Low FRiP | High | Check library quality; may need rerun |
| Missing ChIP Input | High | Input/IgG controls are essential |
| Inconsistent peaks across replicates | Medium | IDR framework; more replicates |
| Array probe coverage bias | Medium | gometh; WGBS for genome-wide view |
| Non-specific histone antibodies | High | Validated antibodies; fingerprint plots |

## 10. Experimental Validation Suggestions

### Methylation validation
```
1. Bisulfite sequencing (BSP/MSP)
   - Validate 3–5 key DMRs
   - BSP cloning for single-CpG quantification

2. Pyrosequencing
   - Precise quantification of selected CpGs

3. 5-Aza demethylation treatment
   - Test whether silenced genes re-express after demethylation
```

### Chromatin validation
```
4. ChIP-qPCR
   - Validate histone marks or TF binding at loci

5. ATAC-qPCR
   - Validate accessibility at selected regions

6. Luciferase reporter
   - Test enhancer/promoter activity; mutate motifs and compare
```

### Functional validation
```
7. CRISPR-Cas9 editing
   - Edit regulatory regions
   - dCas9-DNMT3A / dCas9-TET1 for targeted methylation editing
   - dCas9-p300 / dCas9-KRAB for targeted chromatin modulation

8. Epigenetic drugs
   - DNMT inhibitors (5-Aza, Decitabine)
   - HDAC inhibitors (SAHA, TSA)
   - BET inhibitors (JQ1)
   - Assess target gene expression and phenotypes

9. 3D genome (advanced)
   - Hi-C / 4C-seq for enhancer–promoter loops
   - CRISPRi at enhancers and measure target gene changes
```
