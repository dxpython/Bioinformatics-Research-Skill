# General bioinformatics study design template

> Use this template for studies based on public repositories or your own data. Replace every `[placeholder]` with project-specific information.

---

## 1. Background

### How to write

Provide three layers: broad context → current knowledge in the subfield → the gap your study addresses.

### Template

> [Disease or biological process] is characterized by [incidence/prevalence/severity]. Recent work in [subfield] shows [what is known], but [key gap] remains unclear.
>
> In [tissue/species/cell type], [molecule/pathway/phenotype] has been reported (citations), yet systematic study of its [mechanism/regulatory network/clinical value] is lacking.
>
> With rapid growth in [sequencing/multi-omics/public databases], [RNA-seq / scRNA-seq / ATAC-seq / proteomics / metabolomics] data now enable [integrated analysis/systematic mining].

### Example

> Hepatocellular carcinoma (HCC) is the sixth most common cancer and the third leading cause of cancer death. Immune checkpoint inhibitors (ICIs) have transformed HCC care, yet only ~20% of patients benefit, and predictive biomarkers and underlying mechanisms are poorly defined.
>
> Heterogeneous immune infiltration in the tumor microenvironment (TME) shapes immunotherapy response. Prior work often focused on single subsets (e.g. CD8+ T cells), while multi-cell interaction networks remain incompletely understood.
>
> Large HCC transcriptomic cohorts in TCGA and GEO, together with single-cell atlases, offer an opportunity to dissect TME heterogeneity systematically.

---

## 2. Core scientific question

### How to write

State one **empirically testable** question—not a vague theme. Prefer a single-sentence interrogative form.

### How to derive the question

1. **Phenomenon-driven**: What mechanism explains an observation (e.g. a gene highly expressed in tumors)?
2. **Clinic-driven**: Can we discover biomarkers for [prognosis/response/diagnosis]?
3. **Technology-driven**: Can [new method/data type] reveal [patterns/subpopulations/pathways] we could not see before?
4. **Conflict-driven**: Prior studies disagree—what is the consistent rule?

### Template

> In [disease/tissue/species], how does [molecular process/cell type/pathway] [regulate/influence/predict] [phenotype/prognosis/treatment response]?

### Example

> In HCC, how do cuproptosis-related gene expression programs associate with TME remodeling and patient prognosis?

---

## 3. Hypotheses

### How to write

Explicit hypotheses guide analysis and interpretation. State null (H0) and alternative (H1) hypotheses.

### Template

| Type | Statement |
|------|-----------|
| **H0** | [Target molecule/gene set/feature] has no statistically significant association with [prognosis/phenotype/treatment response]. |
| **H1** | [Target molecule/gene set/feature] significantly [separates/predicts/associates with] [prognosis/phenotype/treatment response], potentially via [pathway/cell type/regulatory mechanism]. |

### Example

| Type | Statement |
|------|-----------|
| **H0** | Cuproptosis-related gene expression is not associated with overall survival (OS) in HCC. |
| **H1** | A cuproptosis-based risk score stratifies HCC patients into high- and low-risk groups that differ in survival, immune infiltration, and drug sensitivity. |

---

## 4. Data sources

### How to write

List every source: repository, data type, sample size, and access path. Transparency supports reproducibility.

### Template

| Source | Repository / cohort | Data type | Sample size | Role | Access |
|--------|---------------------|-----------|-------------|------|--------|
| Public | TCGA-[project] | RNA-seq (TPM/FPKM/counts) | n = [N] | Training | GDC / UCSC Xena |
| Public | GEO [GSE ID] | Microarray / RNA-seq | n = [N] | Validation | GEO website |
| Public | [ICGC, CGGA, ArrayExpress, …] | [Type] | n = [N] | External validation | [How obtained] |
| Internal | [Lab / hospital] | [Type] | n = [N] | Experimental validation | [IRB / ID] |
| Collaborative | [Partner] | [Type] | n = [N] | [Role] | Data agreement |

### Notes

- RNA-seq: raw counts for differential analysis; TPM/FPKM for cross-sample display comparisons.
- Microarray: watch platform differences and batch effects.
- Clinical tables should include survival time and status, TNM stage, age, sex, etc.
- scRNA-seq: record 10x chemistry version and approximate cell numbers.

---

## 5. Analysis workflow

### How to write

Describe steps with goals, methods, and tools. A flowchart often helps.

### Template

```
Step 1: Data acquisition and preprocessing
  ├── Goal: high-quality expression matrix + clinical metadata
  ├── Approach: [download → format conversion → QC → normalization]
  ├── Tools: R (TCGAbiolinks / GEOquery) or Python (pandas)
  └── Outputs: normalized matrix + phenotype table

Step 2: [Core analysis 1, e.g. differential expression]
  ├── Goal: [detect DE genes]
  ├── Model: [e.g. DESeq2 negative binomial]
  ├── Filters: |log2FC| > [cutoff] and adjusted p-value < [cutoff]
  ├── Tools: [R/Python packages]
  └── Outputs: DE table + volcano + heatmap

Step 3: [Core analysis 2, e.g. functional enrichment]
  ├── Goal: [interpret DE genes]
  ├── Methods: [GO/KEGG/GSEA/Reactome]
  ├── Tools: clusterProfiler / GSEA / Metascape
  └── Outputs: enrichment tables + dot/bar plots

Step 4: [Core analysis 3, e.g. prognostic modeling]
  ├── Goal: [build and validate a risk score]
  ├── Methods: [LASSO-Cox / stepwise Cox / random survival forest]
  ├── Validation: [internal CV + external cohorts]
  ├── Tools: glmnet / survival / survminer
  └── Outputs: risk formula + KM + time-dependent ROC + nomogram

Step 5: [Extension, e.g. immune deconvolution]
  ├── Goal: [link risk score to immune context]
  ├── Methods: [CIBERSORT / ssGSEA / ESTIMATE / xCell]
  ├── Tools: [specific packages]
  └── Outputs: cell abundance contrasts + correlation analyses

Step 6: [Experimental validation plan]
  ├── Goal: [validate key findings]
  ├── Methods: [qPCR / Western / IHC, …]
  └── Outputs: wet-lab readouts + statistics
```

---

## 6. Key figures

### How to write

Figures tell the story. Below is a typical checklist for omics papers.

### Suggested figure plan

| Panel | Type | Purpose | Step | Common tools |
|-------|------|---------|------|----------------|
| Fig 1A | Flowchart | Study design overview | Intro | BioRender / draw.io |
| Fig 1B | Venn | Dataset overlap | Steps 1–2 | VennDiagram / ggVennDiagram |
| Fig 2A | Volcano | DE overview | Step 2 | ggplot2 / EnhancedVolcano |
| Fig 2B | Heatmap | Top DE patterns | Step 2 | pheatmap / ComplexHeatmap |
| Fig 3 | Dot/bar enrichment | Functional enrichment | Step 3 | clusterProfiler |
| Fig 4A | Forest plot | Uni/multivariable Cox | Step 4 | forestplot / survminer |
| Fig 4B | KM curves | Risk stratification | Step 4 | survminer |
| Fig 4C | Time-dependent ROC | Discrimination | Step 4 | timeROC / survivalROC |
| Fig 4D | Nomogram | Individualized risk | Step 4 | rms |
| Fig 5 | Box/violin | Immune score contrasts | Step 5 | ggplot2 / ggpubr |
| Fig 6 | Scatter | Risk vs immune metrics | Step 5 | ggplot2 |
| Fig 7 | Wet-lab panels | Experimental validation | Step 6 | GraphPad Prism |

### Supplementary figures (examples)

- PCA / t-SNE / UMAP for QC and batch structure
- Calibration curves for nomogram accuracy
- DCA for clinical net benefit
- GSEA running enrichment score plots

---

## 7. Expected outcomes

### How to write

Anticipate results from literature and hypotheses—this guides interpretation when patterns differ.

### Template

| Step | Expected pattern | If consistent | If not |
|------|------------------|---------------|--------|
| DE | ~[range] DEGs; [signature] mostly up/down | Proceed downstream | Revisit thresholds or QC |
| Enrichment | Enrichment in [pathways, e.g. immune/metabolic] | Supports mechanism story | Explore novel pathways with literature support |
| Prognostic model | AUC > 0.65; KM separation p < 0.05 | Useful model | Tune features, increase N, or change algorithm |
| Immune context | High risk = “cold” immune desert | Explains therapy discordance | Could be “hot” with immune escape—dig deeper |

### Example

> We expect a cuproptosis-based risk score to separate HCC survival. High-risk tumors may show immunosuppressive infiltration (e.g. Tregs, M2 macrophages) and lower predicted benefit from ICIs.

---

## 8. Novelty

### How to write

Novelty need not be revolutionary. Typical angles:

1. **First systematic report** of [gene set/pathway/cell state] in [disease/tissue].
2. **Method combination** integrating [method A] with [method B].
3. **New angle** on a classic problem (e.g. a regulated cell death program).
4. **Translational tool** such as a prognostic nomogram or web resource.
5. **Multi-omics** validation (transcriptome + genome + proteome).
6. **Multi-cohort** replication for robustness.

### Template

> Novel aspects of this work:
> 1. First systematic view of [target] in [disease] from the perspective of [angle].
> 2. A [prognostic/diagnostic] model based on [features], validated in [K] independent cohorts.
> 3. Links between [target] and [immune context/drug response/pathway], informing [clinical implication].

---

## 9. Risks and contingencies

### How to write

List risks and backups before analysis stalls.

| Risk | Description | Mitigation |
|------|-------------|------------|
| Small N | Rare subtype with <50 samples | Merge GEO sets with batch correction; meta-analytic models |
| Data quality | Strong batch effects or heavy missingness | ComBat / SVA; imputation where justified |
| Overfitting | Train performs, test fails | Stronger regularization; fewer features; bootstrap optimism correction |
| No significant DE | Empty DE list | Relax thresholds; nonparametric tests; revisit contrasts |
| Unexpected direction | Gene predicted up is down | Check sample labels; tissue-specific literature |
| Failed external validation | Model does not transfer | Platform effects; simpler models; report as limitation |

---

## 10. Experimental validation

### How to write

Dry-lab claims gain weight with targeted wet-lab follow-up.

| Goal | Method | When to use | Notes |
|------|--------|-------------|-------|
| mRNA validation | qPCR | Key DE genes | Tissues or cell lines; ≥3 biological replicates |
| Protein validation | Western blot | Protein-level change | Antibody specificity + loading controls |
| Tissue localization | IHC | Spatial expression + intensity | Human Protein Atlas can support planning |
| Subcellular localization | IF | Compartmentalization | Colocalization (e.g. Pearson) |
| Functional perturbation | siRNA/shRNA/CRISPR | Causal tests | Measure proliferation, migration, apoptosis, etc. |
| Large-scale protein | TMA | High-N protein–outcome links | Useful for prognostic claims |
| Pathway pharmacology | Inhibitors/activators | Pathway necessity/sufficiency | Include positive/negative controls |

### Minimum expectations from reviewers

1. **Expression**: qPCR for 3–5 core genes (tumor vs normal).
2. **Protein**: WB or IHC for ≥1–2 core targets.
3. **Clinical biospecimens**: ideally ≥20 paired tumor/adjacent samples.
4. **Statistics**: paired t-test for paired samples; Mann–Whitney for unpaired/non-Gaussian data.

### Wet-lab design sketch

```
Validation 1: mRNA of core genes
  ├── Method: qRT-PCR
  ├── Samples: [n] paired tumor/adjacent
  ├── Targets: [Gene1], [Gene2], [Gene3]
  ├── Reference genes: GAPDH and ACTB
  ├── Statistics: paired t-test / Wilcoxon signed-rank
  └── Expectation: [Gene1] and [Gene2] higher in tumor

Validation 2: Protein-level confirmation
  ├── Method: IHC
  ├── Samples: [n] tumor + [n] normal
  ├── Target: [Protein1]
  ├── Scoring: H-score or IRS
  └── Expectation: higher staining in tumor, linked to stage/outcome
```

---

## How to use this document

1. Replace every `[placeholder]` with study-specific text.
2. Skip sections that do not apply to your design.
3. Extend or shorten the workflow and figure plan as needed.
4. Complete a first draft before heavy computation—it is your roadmap.
5. Revise as results arrive.

---

*Template v1.0 | Suited to transcriptomics, proteomics, metabolomics, single-cell omics, and related bioinformatics studies.*
