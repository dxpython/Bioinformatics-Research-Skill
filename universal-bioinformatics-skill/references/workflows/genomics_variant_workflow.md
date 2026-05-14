# Genomic Variant Analysis Workflow

## 1. Applicable Scenarios

1. **Somatic mutation analysis**: SNVs, indels, CNVs from WES/WGS.
2. **TMB / MSI**: Tumor mutational burden and microsatellite instability for immunotherapy stratification.
3. **Mutation signatures**: Deconvolve mutational processes and etiology.
4. **maftools visualization**: Oncoplots, waterfall plots, and related summaries.
5. **GWAS**: Genome-wide association for disease/traits.
6. **Population genetics**: Structure, selection, linkage disequilibrium.
7. **CNV analysis**: Gains and losses genome-wide.
8. **Driver gene discovery**: Distinguish drivers from passengers.

## 2. Input Data Requirements

### Somatic mutation analysis
| Data type | Format | Description |
|-----------|--------|-------------|
| MAF | .maf | TCGA-style mutation table (maftools) |
| VCF | .vcf / .vcf.gz | Variant calls (GATK / Mutect2 / Strelka2) |
| Clinical | TSV | sample_id, cancer_type, stage, survival, treatment |
| CNV segments | .seg | Copy-number segments (GISTIC2) |

### GWAS
| Data type | Format | Description |
|-----------|--------|-------------|
| Genotypes | PLINK (.bed/.bim/.fam) | SNP matrix |
| Phenotype | TSV | Binary or quantitative traits |
| Covariates | TSV | PCs, age, sex, etc. |

**Key requirements**:
- MAF must include: Hugo_Symbol, Chromosome, Start_Position, Variant_Classification, Variant_Type, Tumor_Sample_Barcode.
- VCF should follow GATK Best Practices.
- Tumor–normal pairing recommended to separate somatic vs germline.
- GWAS: suggest ≥1000 samples (cases + controls) combined.
- Consistent reference (hg19/hg38).

## 3. Data Validation Checklist

### 3.1 Mutation QC
```
1. MAF completeness:
   - Required fields present
   - Summarize Variant_Classification categories
   - Mutations per sample

2. Variant filters:
   - VAF >= 0.05 (somatic)
   - Alt depth >= 3 supporting reads
   - Total depth >= 20x
   - Remove PoN artifacts
   - Remove gnomAD/ExAC AF > 0.01

3. Sample quality:
   - Hypermutators or contamination (very high counts)
   - Low mutation burden (low coverage)
   - Tumor purity (ABSOLUTE / ASCAT)
```

### 3.2 GWAS QC
```
1. SNP-level:
   - Call rate > 95% per SNP
   - MAF > 0.01 (or 0.05)
   - HWE p > 1e-6 in controls

2. Sample-level:
   - Call rate > 95% per sample
   - Sex check (X heterozygosity)
   - Relatedness (IBD, pi_hat < 0.2)
   - Population stratification (PCA, remove outliers)

3. Imputation (if used):
   - INFO score > 0.3 (or stricter 0.8)
```

## 4. Recommended Analysis Steps

### Step 1: Somatic mutation landscape
```
1.1 maftools summary
    - read.maf()
    - mafSummary(), plotmafSummary()

1.2 Oncoplot (waterfall)
    - oncoplot() for top mutated genes (top = 20–30)
    - Optional clinical annotation bars
    - Pathway-grouped views

1.3 Mutation frequencies
    - Per-gene frequencies vs known driver lists
    - Variant classes: Missense, Nonsense, Frameshift, Splice_Site, etc.

1.4 Lollipop plots
    - lollipopPlot() for domain hotspots (e.g. TP53, PIK3CA, KRAS)
```

### Step 2: TMB and MSI
```
2.1 TMB
    - Nonsynonymous mutations / targeted region size (Mb)
    - Typical WES capture ~30–50 Mb
    - TMB-high commonly >= 10 mut/Mb (FDA pembrolizumab context)
    - maftools::tmb()

2.2 MSI
    - MSIsensor / MSIsensor2 / MSIsensor-pro; MANTIS
    - MSI-H thresholds depend on tool (e.g. MSIsensor score >= 3.5)
    - MSI-H associated with immunotherapy benefit

2.3 Joint TMB–MSI
    - Partially overlapping information; combined view for IO eligibility
```

### Step 3: Mutation signatures
```
3.1 Trinucleotide context matrix
    - 96-type spectrum; maftools::trinucleotideMatrix()

3.2 Decomposition
    - NMF for de novo signatures
    - MutationalPatterns, SigProfiler, deconstructSigs
    - Match to COSMIC v3.3 reference

3.3 Interpretation (examples)
    - SBS1/5: clock/aging
    - SBS2/13: APOBEC
    - SBS4: tobacco
    - SBS6/15/20/26: MMR deficiency (dMMR)
    - SBS7: UV
    - SBS10a/b: POLE proofreading defect
    - SBS22: aristolochic acid

3.4 Clinical associations
    - Signatures vs stage, prognosis, treatment response
    - Signature exposures vs TMB/MSI
```

### Step 4: CNV analysis
```
4.1 CNV calling
    - WES: GATK CNV / CNVkit / FACETS
    - SNP array: ASCAT / PennCNV
    - Low-pass WGS: QDNAseq / ichorCNA

4.2 GISTIC2
    - Recurrent focal/arm-level events
    - Amplification (score > 0.1, q < 0.25); deletion (score < -0.1, q < 0.25)

4.3 Visualization
    - Genome-wide heatmaps
    - GISTIC plots
    - Single-sample profiles

4.4 Integrate CNV with mutations
    - Amplification + activating mutation
    - Deletion + mutation (two-hit biallelic loss)
```

### Step 5: Driver genes
```
5.1 Statistical tools
    - MutSigCV (background mutation rate)
    - dNdScv (dN/dS)
    - OncodriveFML (functional impact bias)
    - OncodriveCLUSTL (mutation clustering)

5.2 Clinical annotation
    - OncoKB, CIViC, COSMIC Cancer Gene Census

5.3 Pathway view
    - Summarize by oncogenic pathways (RTK/RAS, PI3K/AKT, WNT, MYC, NOTCH, Hippo, NRF2, TGFβ, TP53, cell cycle)
    - maftools::OncogenicPathways()
```

### Step 6: GWAS (when applicable)
```
6.1 Association testing
    - Logistic (case–control) or linear (QTL) models
    - Y ~ SNP + PC1 + PC2 + ... + age + sex
    - PLINK2, REGENIE, SAIGE (large biobanks), BOLT-LMM

6.2 Reporting
    - Genome-wide p < 5e-8; suggestive p < 1e-5
    - Manhattan + QQ plots

6.3 Follow-up
    - Fine-mapping: FINEMAP, SuSiE, PAINTOR
    - Colocalization: coloc, eCAVIAR
    - Mendelian randomization; PRS (PRSice2, LDpred2)
```

## 5. Recommended Statistical Methods

| Goal | Method | Tool |
|------|---------|------|
| Mutation landscape | Frequencies + plots | maftools |
| Differential mutation | Fisher's exact | maftools::mafCompare |
| TMB | Nonsynonymous / Mb | maftools, custom scripts |
| MSI | Microsatellite instability score | MSIsensor2 |
| Signatures | NMF | MutationalPatterns, SigProfiler |
| CNV | Depth + BAF | CNVkit, FACETS |
| Recurrent CNV | GISTIC2 | GISTIC2.0 |
| Drivers | Mutation-rate models | MutSigCV, dNdScv |
| GWAS | Logistic/linear regression | PLINK2, REGENIE |
| Fine-mapping | Credible sets | SuSiE, FINEMAP |
| PRS | C+T / LDpred2 / PRScs | PRSice2, LDpred2 |

## 6. Recommended Figures

| Figure | Purpose | Tool |
|--------|---------|------|
| **Oncoplot / waterfall** | Mutation overview | maftools |
| **Lollipop** | Gene-level hotspots | maftools, trackViewer |
| **TMB barplot** | Compare TMB across samples | ggplot2 |
| **Signature profile** | 96-context spectrum | MutationalPatterns |
| **Signature contributions** | Exposure barplots | ggplot2 |
| **CNV heatmap** | Genome-wide copy number | ComplexHeatmap, CopyNumberPlots |
| **GISTIC** | Recurrent CNV regions | GISTIC2 |
| **Manhattan** | GWAS associations | qqman, CMplot |
| **QQ plot** | Inflation check | qqman |
| **LocusZoom** | Locus detail | locuszoom |
| **Co-mutation** | Mutual exclusivity / co-occurrence | maftools::somaticInteractions |
| **Rainfall** | Mutation density along chromosomes | maftools |
| **VAF density** | VAF distributions | ggplot2 |
| **Forest** | Mutation vs prognosis | forestplot, ggplot2 |

---


## 7. Result Interpretation Template

### Mutation landscape
```
Among [N] [cancer type] samples, [X] nonsynonymous somatic mutations were detected (median [Y]
per sample, range [min]–[max]). Missense mutations were most frequent ([X%]), followed by
Nonsense ([Y%]) and Frameshift ([Z%]). Top mutated genes included [TP53] ([X%]), [PIK3CA] ([Y%]),
and [KRAS] ([Z%]), consistent with known [cancer type] genomics.
```

### TMB / MSI
```
Median TMB was [X] mutations/Mb; [Y] patients ([Z%]) were TMB-high (>= 10 mut/Mb). MSI analysis
identified [N] MSI-high cases. TMB-high and/or MSI-high patients may benefit from immune checkpoint blockade.
```

### Mutation signatures
```
NMF extracted [N] major signatures matching COSMIC SBS references including [SBS1 (aging),
SBS2/13 (APOBEC), SBS[X]]. [X%] of samples were dominated by [signature], suggesting [etiology].
Signature exposure correlated with [clinical feature] (Spearman r = [X], p = [Y]).
```

## 8. Manuscript Writing Template

### Methods
```
Somatic mutation analysis:
Somatic mutations were called from tumor–normal WES with Mutect2 (GATK v[version]) following
Best Practices. Variants were filtered with FilterMutectCalls and additional rules: VAF >= 0.05,
alt depth >= 3, total depth >= 20, gnomAD population AF < 0.01. MAF files were built with vcf2maf
and analyzed in maftools (v[version]).

TMB was defined as nonsynonymous mutations per megabase of capture. MSI was assessed with MSIsensor2.
Signatures were extracted with MutationalPatterns (v[version]) and compared to COSMIC v3.3.
Drivers were evaluated with dNdScv.

GWAS (if applicable):
QC used PLINK2 excluding SNPs with call rate < 95%, MAF < 0.01, or HWE p < 1e-6, and samples with
call rate < 95% or cryptic relatedness (pi_hat > 0.2). PCA assessed stratification. Association used
logistic regression adjusting for top [N] PCs, age, and sex.
```

## 9. Common Issues and Risks

| Issue | Risk | Mitigation |
|-------|------|------------|
| No matched normal | High | PoN filtering; interpret cautiously |
| Low tumor purity | High | ABSOLUTE/ASCAT; affects VAF/CNV |
| Inconsistent TMB definition | Medium | Report panel/WES assumptions clearly |
| Signature overfitting | Medium | Cross-validate number of signatures |
| GWAS population stratification | High | PCA covariates; genomic inflation ~1 |
| GWAS multiple testing | High | Genome-wide 5e-8 threshold |
| Discordant CNV tools | Medium | Multi-tool consensus; check purity/ploidy |
| Passenger mis-called as driver | Medium | Combine statistical and functional evidence |
| Hypermutator outliers | Medium | Flag and analyze separately |

## 10. Experimental Validation Suggestions

### Mutation validation
```
1. Sanger sequencing for key somatic calls (especially low VAF)

2. ddPCR for precise VAF (liquid biopsy / ctDNA)

3. Targeted panel sequencing in independent cohorts for recurrent genes
```

### Functional validation
```
4. Site-directed mutants; compare WT vs mutant protein activity and cellular phenotypes

5. CRISPR/base/prime editing to introduce or correct mutations

6. Drug sensitivity assays (e.g. EGFR mutations vs TKIs)
```

### Clinical validation
```
7. Cross-cohort validation (TCGA/ICGC); compare mutation prevalence across ancestries

8. Prognostic/predictive value vs outcomes and therapies; prospective multi-center studies when possible
```
