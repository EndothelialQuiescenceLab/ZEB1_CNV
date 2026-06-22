# ZEB1_CNV

Analysis code accompanying the manuscript:

> **Loss of endothelial ZEB1 enhances neovascularisation during choroidal neovascularisation.**
> Beazley-Long N\*, Horder JL\*, Green KR\*, Bourne JH, Lynch AP, Ahmed NS, Tabrizi ZB, Diez-Pinel G, Carroll CP, Allen CL, Cresswell CT, Harris PE, Ferreira I, Mosqueira D, Rayes J, Mongan NP, King J, Denning C, McIntyre AP, Stemmler MP, Brabletz S, Brabletz T, Bates DO, Benest AV.
> *Endothelial Quiescence Group, University of Nottingham.* Corresponding author: andrew.benest@nottingham.ac.uk
> \* Equal contribution.

This repository contains all bioinformatics code used to generate the scRNA-seq and bulk RNA-seq analyses in the paper. It does **not** contain raw or processed data — see [Data availability](#data-availability) below.

## Overview of the analyses

The manuscript investigates the role of the transcription factor **ZEB1** in endothelial cells (ECs) during choroidal neovascularisation (CNV), a hallmark of wet age-related macular degeneration. The analyses fall into three blocks:

1. **scRNA-seq re-analysis** of a publicly available murine choroidal EC dataset (Rohlenova et al., 2020) to characterise Zeb1 expression across EC phenotypes and identify Zeb1 enrichment in lesion ECs (manuscript Fig. 1).
2. **In-silico knockout** of Zeb1 in lesion ECs via gene-regulatory-network perturbation (`scTenifoldKnk`), with downstream functional enrichment to predict inflammatory and angiogenic consequences (manuscript Fig. 2).
3. **Bulk RNA-seq** of siRNA-mediated ZEB1 knockdown in HUVECs, with DESeq2 differential expression, functional enrichment, and overlap with the in-silico predictions (manuscript Fig. 5).

## Repository contents

| File | Purpose |
| --- | --- |
| `00_scrnaseq_preprocessing.Rmd` | Builds a Seurat object from the raw, normalised and batch-corrected count matrices, metadata, and t-SNE coordinates downloaded from the [Murine ECTax shiny app](https://endotheliomics.shinyapps.io/murine_ectax/). Writes `00_CNV_raw_norm_batch_object.rds`. |
| `01_scrnaseq_analyses.Rmd` | All downstream scRNA-seq analyses: t-SNE plots, Zeb1 feature/density/violin plots, per-cluster and lesion-vs-non-lesion DE, dot/violin/heatmap visualisations of EndoMT and canonical EC markers, `scTenifoldKnk` in-silico Zeb1 KO on the lesion subset, GO/KEGG/Hallmark enrichment of predicted DR genes, inflammation- and angiogenesis-focused heatmaps and dotplots, and Zeb1+/- stratified DE in tip, lesion, and non-lesion subsets. Produces Figs 1, 2 and supplementary tables 1–2. |
| `02_bulkrnaseq_preprocessing.sh` | Bash pipeline that processes the HUVEC ZEB1 KD / non-silenced FASTQs through FastQC → cutadapt → STAR (GRCh38, Ensembl release 107) → featureCounts. Outputs a gene-level count matrix consumed by `03_bulkrnaseq_analyses.Rmd`. |
| `03_bulkrnaseq_analyses.Rmd` | Full DESeq2 workflow for the HUVEC ZEB1 KD experiment: size factor estimation, BioMart symbol annotation, rlog transformation, PCA, `apeglm`-shrunk DE, GO/KEGG/Hallmark enrichment with a custom non-zero universe, inflammation/angiogenesis-focused dotplots and heatmaps, and orthologue-based overlap with the scTenifoldKnk predictions from `01_*`. Produces Fig. 5 and the relevant supplementary tables. |
| `bulkRNA_preprocessing_env.yml` | Conda environment spec used to run `02_bulkrnaseq_preprocessing.sh` (`fastqc 0.11.9`, `STAR 2.7.9`, `subread 2.0.0`, `cutadapt 4.4`). |
| `bulkRNA_preprocessing_locked.yml` | Fully resolved/locked version of the conda environment for exact reproduction. |
| `renv.lock`, `renv/`, `.Rprofile` | [`renv`](https://rstudio.github.io/renv/) project lockfile and activation scaffolding pinning the R/Bioconductor toolchain (R 4.5.1, Bioconductor 3.21; key packages include `Seurat 5.3.0`, `DESeq2 1.48.2`, `scTenifoldKnk 1.0.2`, `clusterProfiler`, `ComplexHeatmap 2.24.1`, `scCustomize 3.2.0`, `msigdbr`, `biomaRt`). |
| `code.Rproj` | RStudio project file. |
| `.gitignore` | Excludes RStudio user state (`.Rproj.user`, `.Rhistory`, `.RData`, `.Ruserdata`). |

## Data availability

| Dataset | Source |
| --- | --- |
| HUVEC ZEB1 KD bulk RNA-seq (this study) | GEO accession **[GSE333449](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE333449)** |
| Murine choroidal EC scRNA-seq (re-analysed) | [Murine ECTax shiny app](https://endotheliomics.shinyapps.io/murine_ectax/) (Rohlenova et al., *Cell Metab.* 2020) — raw / normalised / batch-corrected matrices, metadata, and t-SNE coordinates |
| Reference genome (bulk RNA-seq) | Ensembl GRCh38, release 107 |

Sample → condition mapping for the HUVEC bulk experiment is documented inline at the top of `02_bulkrnaseq_preprocessing.sh` (3 NS controls, 3 ZEB1 KD, paired-end).

## Expected directory layout

All `*.Rmd` files use relative paths of the form `../data/...` and `../output/...`. Clone this repo into a `code/` subdirectory and create sibling `data/` and `output/` directories alongside it:

```
project_root/
├── code/                              <-- this repo
│   ├── 00_scrnaseq_preprocessing.Rmd
│   ├── 01_scrnaseq_analyses.Rmd
│   ├── 02_bulkrnaseq_preprocessing.sh
│   ├── 03_bulkrnaseq_analyses.Rmd
│   └── ...
├── data/
│   ├── raw_counts/{Data,Metadata}.csv         # from Murine ECTax
│   ├── normalized_data/{Data,Metadata}.csv    # from Murine ECTax
│   ├── batch_corrected_data/{Data,Metadata}.csv
│   ├── Mouse_eye_tSNE_table.csv               # t-SNE coords from Rohlenova et al.
│   └── bulkrnaseq/
│       ├── raw_data/*.fastq                   # HUVEC FASTQs (GSE333449)
│       ├── groups.csv                         # sample-to-condition table
│       └── counts_processed.csv               # post-featureCounts matrix
└── output/
    ├── Robjects/    # serialised Seurat objects and in-silico KO RData
    ├── data/        # exported CSVs (DE tables, enrichment results)
    └── plots/       # SVG figures
```

`02_bulkrnaseq_preprocessing.sh` additionally expects a STAR index for GRCh38 at `~/star_index/human/starindex107HS/` (path can be adjusted at the top of the script).

## Reproducing the analyses

### R environment (scRNA-seq + bulk DE)

The R toolchain is managed with `renv`. From an R session opened at the repo root:

```r
# install.packages("renv")
renv::restore()
```

This installs the exact package versions recorded in `renv.lock` (R 4.5.1, Bioconductor 3.21). `scTenifoldKnk` is RAM-hungry — the in-silico KO step in `01_scrnaseq_analyses.Rmd` was run on a 128 GB machine.

### Conda environment (bulk RNA-seq preprocessing)

```bash
conda env create -f bulkRNA_preprocessing_env.yml      # human-readable spec
# or, for exact reproduction:
conda env create -f bulkRNA_preprocessing_locked.yml
conda activate bulkRNA_preprocessing
```

### Run order

Scripts are numbered and intended to be run in sequence:

1. `00_scrnaseq_preprocessing.Rmd` — knit to build the Seurat object from the Murine ECTax matrices.
2. `01_scrnaseq_analyses.Rmd` — knit to reproduce the scRNA-seq figures and the scTenifoldKnk in-silico KO outputs.
3. `02_bulkrnaseq_preprocessing.sh` — quality control, trimming, alignment, quantification of the HUVEC FASTQs. Edit `main_path` and the STAR index path at the top, then run.
4. `03_bulkrnaseq_analyses.Rmd` — knit to reproduce the bulk DE/enrichment analyses; also consumes `06_cnv_lesion_zeb1KO_DRGs.csv` from step 2 for the overlap with the in-silico KO predictions.

The bulk preprocessing workflow follows the protocol described in [Díez Pinel et al. (2022) *Methods Mol. Biol.* 2441, 369–426](https://pubmed.ncbi.nlm.nih.gov/35099752/).

## Citation

If you use code from this repository, please cite the manuscript above and the original [Rohlenova et al. choroidal scRNA-seq](https://endotheliomics.shinyapps.io/murine_ectax/) dataset for the scRNA-seq analyses.

## Contact

Questions about the code: open an issue, or contact Joseph L. Horder.
Questions about the study: andrew.benest@nottingham.ac.uk
