# Dosage Compensation — Lubna Younas

Analysis code for **"Dosage compensation is regulated by segregated chromatin modules during development in Drosophila"**.

**Repository owner and maintainer:** [Lubna Younas](https://github.com/Lubna-Younas)  
**Repository:** [Lubna-Younas/Dosage-Compensation](https://github.com/Lubna-Younas/Dosage-Compensation)

---

## Overview

This repository contains analysis commands and code for integrating **RNA-seq**, **ChIP-seq**, **Hi-C**, and **single-cell RNA-seq** data to study dosage compensation and chromatin organization in *Drosophila melanogaster*.

The uploaded code covers:

- Read quality control, trimming, alignment, and gene-level quantification.
- ChIP/input enrichment, peak calling, signal profiles, and overlap analysis.
- Hi-C reference preparation, contact-matrix processing, compartments, saddle plots, contact aggregation, distance-decay plots, and TAD calling.
- Single-cell reference preparation, Cell Ranger counting, and Seurat-based visualization.
- Gene-expression, chromatin-signal, and contact-map outputs for downstream comparison.

The analysis files collect multiple command blocks and, in some cases, different programming languages. Use the relevant sections for each analysis after adapting paths, inputs, and computing settings.

---

## Repository Structure

```text
Dosage-Compensation/
├── README.md                     # Project overview and analysis guide
├── ChIP-seq                      # ChIP-seq processing and enrichment analysis
├── RNA-seq                       # RNA-seq processing and gene quantification
├── scRNA-seq data processing      # Single-cell preparation and visualization
└── Hi-C/
    ├── Alignment                 # Reference preparation and HiC-Pro commands
    └── Interactions matrix       # Hi-C configuration and downstream analysis
```

| Analysis file | Main contents |
|---------------|---------------|
| [ChIP-seq](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/ChIP-seq) | fastp quality control; Bowtie2/SAMtools alignment processing; bedtools coverage and gene-level ChIP/input ratios; MACS3 peak calling; deepTools profiles and heatmaps; Intervene overlaps. |
| [RNA-seq](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/RNA-seq) | Trimmomatic processing; RSEM/Bowtie2 quantification; STAR alignment; SAMtools indexing and alignment summaries; featureCounts gene counts. |
| [Hi-C/Alignment](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/Hi-C/Alignment) | Restriction-fragment preparation, chromosome sizes, Bowtie2 reference indexing, and HiC-Pro processing commands. |
| [Hi-C/Interactions matrix](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/Hi-C/Interactions%20matrix) | HiC-Pro configuration; contact matrices; Cooler/cooltools and FAN-C analyses; HiCExplorer aggregation, distance-decay plots, and TAD calling. |
| [scRNA-seq data processing](https://github.com/Lubna-Younas/Dosage-Compensation/blob/Dosage-Compensation/scRNA-seq%20data%20processing) | StringTie transcriptome preparation, Cell Ranger reference construction and counting, Seurat plotting functions, and testis-expression visualizations. |

---

## System Requirements

### Hardware and Computing Environment

The code includes **SLURM** job examples and paths from the computing environment used for the analyses. CPU, memory, storage, and runtime requirements depend on the dataset and selected analysis.

Before submitting a job:

- Update input, output, reference, and software paths.
- Match the requested CPUs and memory to the tool settings in that command block.
- Set the SLURM partition, runtime, and other scheduler options for the available cluster.
- Provide enough storage for sequencing inputs, intermediate alignments, and matrix outputs.

### Software Dependencies

Install or load the tools required for the analysis being used.

| Analysis | Software used in the uploaded code |
|----------|------------------------------------|
| Read quality control and trimming | fastp, Trimmomatic |
| Alignment and reference indexing | Bowtie2, STAR, SAMtools |
| RNA-seq quantification | RSEM, Subread/featureCounts |
| ChIP-seq analysis | bedtools, MACS3, deepTools, Intervene |
| Hi-C processing | HiC-Pro |
| Hi-C matrix and chromatin-organization analysis | Cooler, cooltools, FAN-C, HiCExplorer |
| Single-cell reference preparation and counting | StringTie, Cell Ranger |
| Restriction-fragment Python code | Python, NumPy |
| Single-cell visualization | R, Seurat, SeuratObject, ggplot2, dplyr, tidyverse, DT, ggpubr, RColorBrewer |

Executable paths and module names in the files identify the original local installations. Update them for the environment where the analysis will run.

---

## Installation

### Step 1: Clone the Repository

```bash
git clone https://github.com/Lubna-Younas/Dosage-Compensation.git
cd Dosage-Compensation
```

### Step 2: Prepare the Software Environment

Use the software available on the cluster or prepare an environment containing the dependencies needed for the selected analysis.

For example, on a system providing these environment modules:

```bash
module load bowtie2
module load samtools
module load bedtools
```

Load the additional tools listed in the relevant analysis file. The repository currently provides analysis command collections; environment files and installation scripts are not included.

### Step 3: Prepare Analysis Jobs

Copy the required command block into a separate shell, SLURM, Python, or R script as appropriate. Supply the input files and update paths and parameters before running that script.

Some repository files contain several job examples or combine shell commands with Python or R code, so select the relevant sections separately.

---

## Data Requirements

### Input Data

| Data type | Typical format | Used for |
|-----------|----------------|----------|
| RNA-seq reads | Paired-end FASTQ | Trimming, alignment, expression quantification, and gene counting |
| ChIP-seq and matching input reads | Paired-end FASTQ | Alignment, enrichment comparisons, and peak calling |
| Hi-C reads | Paired-end FASTQ | HiC-Pro processing and contact-matrix construction |
| Single-cell reads | Assay-compatible FASTQ | Cell Ranger counting |
| Reference genome | FASTA | Alignment indices, restriction fragments, and single-cell reference construction |
| Gene annotation | GTF/GFF, as required by the selected tool | Gene quantification, transcriptome preparation, and feature-based analysis |
| Genomic regions and gene lists | BED or text | Signal aggregation, overlaps, and gene selection |
| Processed expression objects | RDS, CSV, or count tables | Seurat and expression visualizations |

### Reference and Sample Preparation

Use compatible genome and annotation versions throughout each analysis, including consistent chromosome names. Prepare the required reference indices and organize samples with clear identifiers for assay, sex, tissue, stage, and replicate where applicable.

The sequencing data, reference files, processed matrices, and sample-specific inputs are supplied separately from this code repository.

---

## Usage

### Step 1: Pre-process Raw Sequencing Data

#### ChIP-seq Quality Control and Alignment

Use the quality-control and alignment sections in **ChIP-seq** to:

1. Process paired-end reads with fastp.
2. Align reads with Bowtie2.
3. Filter and sort alignments with SAMtools.
4. Generate duplicate-removed BAM files and BAM indices.
5. Calculate per-base coverage with bedtools.

**Outputs shown in the code:** trimmed FASTQ files, sorted and duplicate-removed BAM files, BAM indices, and coverage-count files.

#### RNA-seq Read Preparation

The **RNA-seq** file contains a Trimmomatic command example for paired-end read preparation. Adapt the sample filenames, adapter file, quality encoding, and output settings for the reads being processed.

---

### Step 2: Quantify RNA-seq Expression

The **RNA-seq** file provides two analysis routes:

| Route | Main steps | Outputs shown in the code |
|-------|------------|---------------------------|
| RSEM with Bowtie2 | Quantify expression against an RSEM reference | RSEM expression outputs and alignment files requested by the command |
| STAR with featureCounts | Align reads, index the BAM, summarize alignments, and count annotated exons by gene | Sorted BAM and BAM index, alignment report, featureCounts table, and extracted gene-count table |

The featureCounts examples write files using the following suffixes:

```text
.exon.gene_name.txt
.exon.gene_name.counts.txt
```

Use sample-specific output prefixes and the annotation appropriate for the reference genome.

---

### Step 3: Process Hi-C Data and Analyze Chromatin Organization

#### Prepare the Reference

The **Hi-C/Alignment** file includes:

- Restriction-fragment generation using a `^GATC` recognition-site example.
- FASTA indexing and chromosome-size extraction.
- Bowtie2 reference-index construction.
- HiC-Pro processing and cluster-submission commands.

The restriction digest example produces:

```text
Dmel_digested.genome.bed
```

Use the restriction sites and ligation settings appropriate for the Hi-C experiment.

#### Configure and Run HiC-Pro

The **Hi-C/Interactions matrix** file contains the configuration text and processing examples. Prepare a separate configuration file, such as the `config-hicpro_Dmel.txt` referenced by the commands, using the settings for the experiment.

Review read-pair naming, reference files, restriction fragments, ligation site, matrix resolutions, filtering, and normalization settings.

#### Analyze Contact Matrices

| Analysis | Tools represented in the code |
|----------|-------------------------------|
| Matrix conversion, chromosome selection, and balancing | Cooler, HiCExplorer, FAN-C |
| Expected contacts, compartments, and saddle plots | cooltools, FAN-C |
| Male/female contact-map comparison | FAN-C |
| Aggregate contact analysis | HiCExplorer |
| Contact frequency versus genomic distance | HiCExplorer |
| TAD calling | HiCExplorer |

**Outputs shown in the code:** contact matrices, `.cool` files, expected-contact tables, compartment vectors and strength outputs, aggregate-contact plots, distance-decay plots/data, and TAD output files.

---

### Step 4: Analyze ChIP-seq Enrichment and Genomic Overlaps

The downstream sections in **ChIP-seq** cover:

1. Summarizing ChIP and input coverage across annotated regions with bedtools.
2. Calculating gene-level ChIP/input ratios.
3. Calling peaks with MACS3 using a matching input BAM.
4. Creating coverage and comparison bigWig tracks with deepTools.
5. Building scaled-region matrices and plotting heatmaps.
6. Examining overlaps between genomic-region lists with Intervene.

**Outputs shown in the code:** gene-level coverage and ratio tables, peak-call outputs, bigWig files, compressed signal matrices, heatmap PDFs, and overlap-analysis outputs.

---

### Step 5: Prepare Single-cell Data and Generate Visualizations

Use the relevant sections in **scRNA-seq data processing** for:

- StringTie transcript assembly and merging.
- Selecting transcript features and preparing an augmented annotation.
- Cell Ranger reference construction and counting.
- Gene-name table preparation.
- Seurat-based UMAP, violin, feature, heatmap, and dot-plot visualizations.
- Testis cell-type and expression summaries.

The code uses sample-specific BAMs, gene lists, annotation files, FASTA files, and expression objects. It also references local helper scripts such as `merged_and_ordered.R` and `seurmo.r`; these must be supplied separately if those command blocks are used.

The plotting functions can be adapted to the cell annotations, gene sets, and expression objects for the analysis.

---

## Statistical Analysis

The Hi-C TAD-calling example in **Hi-C/Interactions matrix** specifies `--correctForMultipleTesting fdr`.

Additional manuscript-level statistical analyses should be documented alongside their corresponding data and analysis code. The current repository files provide processing, quantification, chromatin-analysis, and visualization examples; they do not include a separate manuscript-wide statistical-testing script.

---

## Troubleshooting

### Common Issues

| Issue | Check |
|-------|-------|
| File or directory not found | Replace the original cluster paths and sample filenames with existing inputs and output locations. |
| Software command not found | Load the required module or activate the environment containing that tool; update absolute executable paths. |
| A file fails when run as one script | Select the relevant command blocks and separate shell, SLURM, Python, and R sections into appropriate scripts. |
| SLURM job fails or runs out of resources | Inspect the job logs and match CPU, memory, and runtime requests to the command and dataset. |
| Reference or annotation mismatch | Confirm the genome version, chromosome names, indices, and annotation are compatible. |
| HiC-Pro configuration or input error | Check paired-read naming, reference paths, restriction fragments, ligation site, and matrix settings. |
| Single-cell helper script or input is missing | Supply the local helper scripts and sample-specific files referenced by the selected section. |
| Seurat plot uses an unavailable object or metadata field | Load the intended expression object and adapt object names, features, and metadata columns. |

### Getting Help

For questions or bug reports, open an issue in [Lubna Younas's Dosage-Compensation repository](https://github.com/Lubna-Younas/Dosage-Compensation/issues).

---

## Citation

If you use this code or data, please cite:

> Lubna Younas et al.  
> Dosage compensation is regulated by segregated chromatin modules during development in Drosophila.  
> *Molecular Systems Biology*. [In press]

**Repository:** Lubna Younas, *Dosage-Compensation*.  
[GitHub repository](https://github.com/Lubna-Younas/Dosage-Compensation)

---

## Contact

**Lubna Younas — First Author and Repository Maintainer**

- **Email:** [lubna.ma528@gmail.com](mailto:lubna.ma528@gmail.com)
- **GitHub:** [Lubna Younas](https://github.com/Lubna-Younas)
- **Repository:** [Dosage-Compensation](https://github.com/Lubna-Younas/Dosage-Compensation)
- **Support:** [Repository issues](https://github.com/Lubna-Younas/Dosage-Compensation/issues)

---

## Acknowledgements

This work was supported by the European Research Council Starting Grant.
