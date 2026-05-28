# Pathway Overlap Networks for Single-Cell RNA-seq

**Blog:** Compute and Conquer  
**Author:** Mike Martinez  
**Date:** May 2026

---

## Overview

This repository contains the code and environment for a blog post demonstrating how to build **pathway overlap networks** from single-cell RNA-seq differential expression data. The analysis uses the IFNB PBMC dataset (Kang et al. 2018, *Nature Biotechnology*) — ~14,000 human peripheral blood mononuclear cells profiled by 10x Genomics scRNA-seq under two conditions: untreated control and interferon-beta (IFN-β) stimulation.

The core idea: rather than interpreting hundreds of GO terms from an over-representation analysis as a flat list, we construct a gene-sharing network where GO terms are nodes and edges connect terms with overlapping gene membership (Jaccard similarity). Leiden community detection then groups redundant, related terms into coherent biological modules that are easier to interpret.

The full pipeline — differential expression, GO over-representation analysis, and network construction — is implemented using the [GDSCtools](https://github.com/Dartmouth-Data-Analytics-Core/GDSCtools) R package.

---

## Repository Structure

```
MayBlog/
├── code/
│   └── 01-IFNB-Networks.R      # Full analysis: DE → ORA → network
├── envs/
│   ├── pixi.toml               # Conda environment specification
│   └── pixi.lock               # Locked dependency versions
├── outputs/
│   └── Ctrl_v_Stim/
│       ├── DEGs/               # Per-cell-type differential expression results
│       ├── ORA/                # Per-cell-type GO over-representation results
│       └── networks/           # Network CSVs and PNGs per cell type
└── MayBlog.Rproj               # RStudio project file
```

---

## Analysis Pipeline

### 1. Differential Expression
Seurat's `FindMarkers` is run for each annotated cell type, comparing IFN-β-stimulated cells (`ident.1 = "STIM"`) against controls (`ident.2 = "CTRL"`). Positive log2FC values indicate higher expression in STIM.

### 2. GO Over-Representation Analysis
Significant DE genes are tested for GO Biological Process enrichment using clusterProfiler's `enrichGO`, with all genes in the Seurat object as the background universe.

### 3. Pathway Overlap Network
For each cell type, significant GO terms are organized into a weighted network:
- **Nodes:** GO terms (sized by number of DE genes annotated; colored by up/down regulation)
- **Edges:** Jaccard similarity between gene sets, thresholded at ≥ 0.1
- **Clustering:** Leiden community detection (resolution = 0.03) groups semantically related terms into modules

---

## Environment Setup

This project uses [pixi](https://pixi.sh) for reproducible environment management.

```bash
# Install pixi (if needed)
curl -fsSL https://pixi.sh/install.sh | bash

# Install the conda environment
cd envs/
pixi install

# Install GitHub-only packages
pixi run install-git-1   # SeuratData
pixi run install-git-2   # ggwordcloud, GDSCtools

# Launch RStudio with the project
pixi run rstudio
```

> **Note:** If you encounter a `bioconductor-genomeinfodbdata` post-link error, run:
> ```bash
> pixi run fix-bioc-data
> ```

---

## Key Dependencies

| Package | Source | Purpose |
|---------|--------|---------|
| Seurat ≥ 5.0 | CRAN | Single-cell analysis |
| SeuratData | GitHub | IFNB dataset |
| clusterProfiler | Bioconductor | GO over-representation analysis |
| org.Hs.eg.db | Bioconductor | Human gene annotation |
| igraph | CRAN | Network construction and Leiden clustering |
| GDSCtools | GitHub | DE → ORA → network pipeline wrapper |

---

## Data

The IFNB dataset is downloaded automatically by `SeuratData::InstallData("ifnb")` when the script is first run. No manual data download is required.

---

## Reference

Kang et al. (2018). Multiplexing droplet-based single cell RNA-sequencing using natural genetic variation. *Nature Biotechnology*, 36, 89–94.
