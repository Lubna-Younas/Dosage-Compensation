# Dosage Compensation

Analysis scripts supporting **"Dosage compensation is regulated by segregated chromatin modules during development in Drosophila"**.

Maintained by [Lubna Younas](https://github.com/Lubna-Younas).

## Overview

This repository contains analysis code for RNA-seq, ChIP-seq, Hi-C, and single-cell RNA-seq in *Drosophila melanogaster*, used to study dosage compensation and chromatin organization.

## Repository contents

| File | Contents |
|------|----------|
| [ChIP-seq](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/ChIP-seq) | Read quality control and trimming, Bowtie2 alignment, ChIP/input comparisons, MACS3 peak calling, deepTools profiles and heatmaps, and Intervene overlaps. |
| [Hi-C/Alignment](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/Hi-C/Alignment) | Restriction-fragment preparation and HiC-Pro processing commands. |
| [Hi-C/Interactions matrix](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/Hi-C/Interactions%20matrix) | HiC-Pro configuration and processing, contact-matrix analysis, compartments and saddle plots, contact aggregation, distance-decay plots, and TAD calling. |
| [RNA-seq](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/RNA-seq) | Trimmomatic processing, RSEM quantification, STAR alignment, and featureCounts gene counts. |
| [scRNA-seq data processing](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/scRNA-seq%20data%20processing) | StringTie transcriptome preparation, Cell Ranger reference construction and counting, and Seurat/R visualization code. |

## Get the repository

```bash
git clone https://github.com/Lubna-Younas/Dosage-Compensation.git
cd Dosage-Compensation
```

## Using the code

The files collect analysis commands and code sections, including SLURM job examples, shell commands, Python code, and R code. Select the sections needed for an analysis and adapt input files, reference files, software paths, and scheduler settings to the computing environment.

Software used in the code includes:

- **ChIP-seq:** fastp, Bowtie2, SAMtools, bedtools, MACS3, deepTools, and Intervene.
- **Hi-C:** HiC-Pro, NumPy, Cooler, cooltools, FAN-C, and HiCExplorer.
- **RNA-seq:** Trimmomatic, RSEM, Bowtie2, STAR, SAMtools, and Subread/featureCounts.
- **Single-cell RNA-seq:** StringTie, Cell Ranger, Seurat, SeuratObject, ggplot2, dplyr, tidyverse, DT, ggpubr, and RColorBrewer.

## Citation

If you use this code or data, please cite:

> Lubna Younas, Mujahid Ali, Xinpei Zhang, Huangyi He, Catherine Regnard, and Qi Zhou.
> Dosage compensation is regulated by segregated chromatin modules during development in Drosophila.
> *Molecular Systems Biology*. [In press]

## Contact and support

**Lubna Younas — First Author**

- Email: [lubna.ma528@gmail.com](mailto:lubna.ma528@gmail.com)
- GitHub: [Lubna Younas](https://github.com/Lubna-Younas)
- Repository: [Dosage-Compensation](https://github.com/Lubna-Younas/Dosage-Compensation)
- Questions and bug reports: [Repository issues](https://github.com/Lubna-Younas/Dosage-Compensation/issues)

## Acknowledgements

This work was supported by the European Research Council Starting Grant.
