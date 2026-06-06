# From NCBI GEO to KEGG Pathway Analysis
## Investigating Bufalin's Potential to Mitigate Doxorubicin-Induced Cardiotoxicity

> A complete, reproducible R bioinformatics pipeline — from public gene expression data (NCBI GEO) to differential expression analysis and KEGG pathway enrichment — focusing on **doxorubicin (DOX) cardiotoxicity** and the cardiac glycoside **bufalin**.

---

## Table of Contents

- [Overview](#overview)
- [Scientific Background](#scientific-background)
- [Pipeline Architecture](#pipeline-architecture)
- [Prerequisites](#prerequisites)
- [Repository Structure](#repository-structure)
- [How to Run](#how-to-run)
- [Analysis Steps](#analysis-steps)
- [Key Results & Interpretation](#key-results--interpretation)
- [Adapting This Pipeline](#adapting-this-pipeline)
- [References](#references)
- [License](#license)

---

## Overview

This repository provides a fully reproducible, end-to-end bioinformatics workflow to:

- Download microarray/RNA-seq data from **NCBI GEO**
- Perform quality control and normalization
- Run differential expression analysis with **limma** (or DESeq2 for RNA-seq)
- Conduct **KEGG Over-Representation Analysis (ORA)** and **Gene Set Enrichment Analysis (GSEA)**
- Visualize perturbed pathways with **pathview**
- Explore mechanistic overlaps between DOX cardiotoxicity and **bufalin** (Na⁺/K⁺-ATPase inhibitor)

**Main Scientific Question:**
> Can bufalin reduce doxorubicin-induced cardiotoxicity? Through which molecular pathways (calcium signaling, cardiac muscle contraction, apoptosis, oxidative stress) might it exert a protective — or additive toxic — effect?

---

## Scientific Background

| Agent | Mechanism | Cardiac Effect |
|---|---|---|
| **Doxorubicin (DOX)** | Topoisomerase II inhibitor, ROS generation | Oxidative stress, calcium overload, mitochondrial dysfunction, apoptosis |
| **Bufalin** | Na⁺/K⁺-ATPase inhibitor (cardiac glycoside) | Modulates intracellular Ca²⁺ via Na⁺/Ca²⁺ exchanger (NCX1/SLC8A1) |

**Hypothesis:** At sub-toxic concentrations, bufalin's modulation of calcium homeostasis may counteract DOX-induced calcium overload. Alternatively, combined inhibition of Na⁺/K⁺-ATPase and ROS production could be additive — or even synergistically toxic.

This pipeline provides the computational framework to test these hypotheses using publicly available transcriptomic data.

---

## Pipeline Architecture

```
NCBI GEO (GSE59672)
        │
        ▼
01_download_and_qc.R
  ├── GEOquery::getGEO()
  ├── Boxplots, density plots, PCA
  └── Normalization (quantile / RMA)
        │
        ▼
02_differential_expression.R
  ├── limma: design matrix, contrast matrix
  ├── eBayes, topTable
  ├── Volcano plot, heatmap
  └── Gene list export (ENTREZ IDs)
        │
        ▼
03_kegg_analysis.R
  ├── ORA: enrichKEGG (over-representation)
  ├── GSEA: gseKEGG (ranked gene list)
  ├── Dot plots, enrichment plots
  └── pathview pathway maps
        │
        ▼
04_bufalin_mechanism.R
  ├── Candidate gene overlap (ATP1A1, SLC8A1, RYR2, ...)
  ├── Sub-network analysis
  └── Hypothesis summary table
```

---

## Prerequisites

### System Requirements

- **R** ≥ 4.2
- **RStudio** (recommended)
- ~2 GB disk space for GEO data and pathway images

### R Packages

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GEOquery",        # GEO data download
  "limma",           # Differential expression (microarray)
  "DESeq2",          # Differential expression (RNA-seq, optional)
  "clusterProfiler", # ORA and GSEA
  "pathview",        # KEGG pathway visualization
  "org.Mm.eg.db",   # Mouse gene annotation
  "org.Hs.eg.db",   # Human gene annotation (for adaptation)
  "AnnotationDbi",   # Gene ID conversion
  "enrichplot",      # Enrichment visualization
  "DOSE",            # Disease enrichment
  "KEGGREST",        # KEGG REST API access
  "ComplexHeatmap",  # Advanced heatmaps (optional, improves figures)
  "fgsea"            # Fast GSEA (optional alternative)
))

# CRAN packages
install.packages(c(
  "dplyr", "ggplot2", "pheatmap", "ggrepel",
  "tidyr", "stringr", "patchwork", "RColorBrewer"
))
```

---

## Repository Structure

```
bufalin-dox-cardiotoxicity-pipeline/
│
├── README.md
│
├── scripts/
│   ├── 00_install_packages.R          # One-time package installation
│   ├── 01_download_and_qc.R           # GEO download + quality control
│   ├── 02_differential_expression.R   # limma DEA + volcano + heatmap
│   ├── 03_kegg_analysis.R             # ORA + GSEA + pathview
│   └── 04_bufalin_mechanism.R         # Bufalin-specific exploration
│
├── reports/
│   └── Cardiotoxicity_Bufalin_Report.Rmd  # Full reproducible report (HTML/PDF)
│
├── results/
│   ├── tables/
│   │   ├── DEG_DOX_vs_CTL.csv         # Differentially expressed genes
│   │   ├── KEGG_ORA_results.csv       # ORA enrichment table
│   │   └── KEGG_GSEA_results.csv      # GSEA enrichment table
│   ├── figures/
│   │   ├── PCA_plot.png
│   │   ├── volcano_plot.png
│   │   ├── heatmap_top50.png
│   │   ├── dotplot_ORA.png
│   │   └── gsea_enrichment_plots/
│   └── pathview/
│       ├── mmu04260.png               # Cardiac muscle contraction
│       ├── mmu04020.png               # Calcium signaling
│       └── ...
│
└── data/                              # Auto-generated; not tracked by git
    └── GSE59672/
```

> **Note:** The `data/` directory is excluded from version control (add to `.gitignore`). Raw GEO data is downloaded automatically by `01_download_and_qc.R`.

---

## How to Run

### Option A — Script by script

```bash
# 1. Clone the repository
git clone https://github.com/yourusername/bufalin-dox-cardiotoxicity-pipeline.git
cd bufalin-dox-cardiotoxicity-pipeline

# 2. Open in RStudio and run scripts in order
```

```r
source("scripts/00_install_packages.R")   # First time only
source("scripts/01_download_and_qc.R")
source("scripts/02_differential_expression.R")
source("scripts/03_kegg_analysis.R")
source("scripts/04_bufalin_mechanism.R")
```

### Option B — Full R Markdown report

```r
rmarkdown::render("reports/Cardiotoxicity_Bufalin_Report.Rmd",
                  output_format = "html_document",
                  output_file   = "results/Cardiotoxicity_Bufalin_Report.html")
```

This generates a single self-contained HTML report with all figures, tables, and interpretation.

---

## Analysis Steps

### 1. Data Download & Quality Control (`01_download_and_qc.R`)

- **Dataset:** [GSE59672](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59672) — mouse model of acute DOX-induced cardiotoxicity
- Download via `GEOquery::getGEO("GSE59672", GSEMatrix = TRUE)`
- Normalization: RMA (microarray) or TMM/DESeq2 normalization (RNA-seq)
- QC outputs: boxplots of expression distribution, density plots, PCA to confirm group separation

### 2. Differential Expression Analysis (`02_differential_expression.R`)

- **Framework:** `limma` with `voom` (if count data) or directly on normalized intensities
- **Model:** `~0 + group` (no intercept, facilitates contrasts)
- **Contrast:** `DOX - CTL`
- **Thresholds:** `adj.P.Val < 0.05` AND `|logFC| > 1`
- **Outputs:** volcano plot, heatmap of top 50 DEGs, full DEG table with ENTREZ IDs

### 3. KEGG Pathway Analysis (`03_kegg_analysis.R`)

#### Over-Representation Analysis (ORA)
Uses the significant DEG list (binary input).

```r
ego <- enrichKEGG(gene         = sig_entrez,
                  organism     = "mmu",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05)
```

#### Gene Set Enrichment Analysis (GSEA)
Uses the full ranked gene list (logFC-ranked), capturing both up and down regulation.

```r
gse <- gseKEGG(geneList     = ranked_list,  # named numeric vector, sorted by logFC
               organism     = "mmu",
               minGSSize    = 15,
               pAdjustMethod = "BH",
               seed         = 42)
```

#### Pathways of Interest

| KEGG ID | Pathway | Relevance |
|---|---|---|
| mmu04260 | Cardiac muscle contraction | DOX direct target |
| mmu04020 | Calcium signaling pathway | Bufalin mechanism |
| mmu04261 | Adrenergic signaling in cardiomyocytes | Cardiac stress response |
| mmu04210 | Apoptosis | DOX-induced cell death |
| mmu04932 | Non-alcoholic fatty liver disease | Mitochondrial dysfunction overlap |
| mmu05410 | Hypertrophic cardiomyopathy | Chronic cardiotoxicity |
| mmu05414 | Dilated cardiomyopathy | Chronic cardiotoxicity |
| mmu04370 | VEGF signaling pathway | Cardioprotective signaling |

#### Pathway Visualization

```r
pathview(gene.data  = logFC_vector,
         pathway.id = "mmu04260",
         species    = "mmu",
         out.suffix = "DOX_cardiotoxicity")
```

### 4. Bufalin Mechanism Exploration (`04_bufalin_mechanism.R`)

- Extract expression values for **bufalin target genes** from the DOX dataset:

| Gene | Protein | Role |
|---|---|---|
| `Atp1a1` | Na⁺/K⁺-ATPase α1 | Direct bufalin target |
| `Slc8a1` | NCX1 | Na⁺/Ca²⁺ exchanger |
| `Ryr2` | RyR2 | Ryanodine receptor (SR Ca²⁺ release) |
| `Cacna1c` | Cav1.2 | L-type Ca²⁺ channel |
| `Serca2a (Atp2a2)` | SERCA2a | SR Ca²⁺ reuptake |
| `Bcl2`, `Bax` | BCL-2/BAX | Apoptosis regulation |
| `Nfe2l2` | NRF2 | Oxidative stress response |

- Overlap analysis: which of these are differentially expressed under DOX?
- Interpret directionality: would bufalin's effect on these genes oppose or amplify DOX-induced dysregulation?

---

## Key Results & Interpretation

### Expected Findings

Under DOX treatment, you should observe significant dysregulation of:
- **Calcium handling genes** (downregulation of `Atp2a2`/SERCA2a, `Ryr2`)
- **Mitochondrial and oxidative stress genes** (upregulation of `Nox4`, downregulation of `Sod2`)
- **Sarcomeric proteins** (`Myh6`, `Myh7`, `Tnni3`)
- **Apoptosis markers** (upregulation of `Bax`, downregulation of `Bcl2`)

### Bufalin Hypothesis Assessment

| Scenario | Evidence Pattern | Interpretation |
|---|---|---|
| **Protective** | Bufalin reverses Ca²⁺ overload genes; NRF2 pathway activated | Low-dose bufalin may be cardioprotective |
| **Neutral** | No overlap between DOX DEGs and bufalin targets | Different mechanisms, no interaction |
| **Additive toxicity** | Bufalin targets are already dysregulated by DOX in the same direction | Risk of combined cardiotoxicity |

> **Important:** Computational results must be validated with in vitro (cardiomyocyte cell lines) and in vivo experiments before drawing clinical conclusions.

---

## Adapting This Pipeline

| Goal | Change |
|---|---|
| Different GEO dataset | Change `GEO_ID <- "GSE59672"` in script 01 |
| Human data | Replace `org.Mm.eg.db` → `org.Hs.eg.db` and `organism = "hsa"` |
| RNA-seq input | Use `DESeq2` in script 02, `voom()` in limma, or switch to edgeR |
| Different drug | Update candidate gene list in script 04 |
| Add bufalin-treated dataset | Add a second GEO accession and compare contrasts in script 02 |
| Reactome instead of KEGG | Replace `enrichKEGG` → `enrichPathway` (ReactomePA package) |

---

## References & Resources

- **Dataset:** [GSE59672 — NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59672)
- **KEGG Database:** https://www.kegg.jp
- **clusterProfiler:** [Bioconductor vignette](https://bioconductor.org/packages/clusterProfiler)
- **limma:** [Bioconductor user guide](https://bioconductor.org/packages/limma)
- **pathview:** [Bioconductor vignette](https://bioconductor.org/packages/pathview)
- Lipshultz et al. (2010). Doxorubicin cardiotoxicity — mechanisms and prevention. *Lancet Oncology*
- Zhang et al. (2022). Bufalin as a cardioprotective agent — evidence and mechanisms. *Front. Pharmacol.*

---

## License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for details.

---

## Contributing

Contributions are welcome! Feel free to open an issue or pull request to:

- Add support for additional datasets (bufalin-treated, combination DOX+bufalin)
- Implement an RNA-seq branch (DESeq2-based)
- Add network/PPI analysis (STRINGdb, igraph)
- Add a Shiny app for interactive pathway exploration

---

*Made with ❤️ for toxicogenomics and drug discovery research.*  
*Questions or need help adapting the code? [Open an issue](../../issues)!*
