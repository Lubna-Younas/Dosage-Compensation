# chromatin_3D

Analysis pipeline for **"Developmental and Epigenetic Variations of Drosophila Dosage Compensated Genes with Their Distinct Chromatin Architecture and Segregated Modules"**

---

## Overview

This repository contains custom scripts and analysis pipelines for integrating **RNA-seq**, **ChIP-seq**, and **Hi-C** data to study dosage compensation (DC) and 3D chromatin organization in *Drosophila melanogaster*. The pipeline enables:

- Classification of dosage-compensated genes (CDC, FDC, SDC)
- Hi-C matrix processing, normalization, and compartment/TAD calling
- Long-range contact (LRC) and short-range contact (SRC) quantification
- ChIP-seq peak calling, aggregation, and correlation with expression
- Integration of sequence elements (PionX, HAS, 1.688³F) with chromatin architecture
- Visualization (heatmaps, browser tracks, box plots, network plots)

---

## Repository Structure

```
chromatin_3D/
|-- README.md                          # This file
|-- requirements.txt                   # Python dependencies
|-- environment.yml                    # Conda environment (optional)
|-- scripts/
|   |-- 01_preprocess/
|   |   |-- hicpro_config.conf        # HiC-Pro configuration template
|   |   |-- run_hicpro.sh             # Hi-C processing wrapper
|   |   |-- chipseq_pipeline.sh       # ChIP-seq alignment and peak calling
|   |   |-- rnaseq_pipeline.sh        # RNA-seq quantification pipeline
|   |
|   |-- 02_dosage_compensation/
|   |   |-- classify_dc_genes.R       # Classify CDC/FDC/SDC from ChIP data
|   |   |-- compare_male_female.R     # Male vs. female expression analysis
|   |   |-- tau_score.R               # Tissue specificity (TAU) calculation
|   |
|   |-- 03_hic_analysis/
|   |   |-- call_compartments.py      # PC1 eigenvector calculation
|   |   |-- call_tads.py              # TAD boundary calling (HiCExplorer)
|   |   |-- quantify_lrcs.py          # LRC/SRC frequency quantification
|   |   |-- call_loops_mustache.py    # Loop calling with Mustache
|   |   |-- insulation_score.py       # Boundary insulation profiles
|   |
|   |-- 04_chipseq_analysis/
|   |   |-- deeptools_aggregation.sh  # Anchor-centered ChIP aggregation
|   |   |-- correlate_chip_expression.R  # ChIP-RNA correlation analysis
|   |   |-- plot_enrichment_profiles.R   # Gene body enrichment plots
|   |
|   |-- 05_figures/
|   |   |-- generate_figures.R        # Main and supplementary figure generation
|   |   |-- plot_browser_tracks.py    # pyGenomeTracks locus-level views
|   |   |-- go_enrichment.R           # GO term enrichment and network plots
|   |
|   |-- 06_statistics/
|       |-- wilcoxon_fdr.R            # Wilcoxon tests with BH-FDR correction
|       |-- permutation_test.py       # Permutation-based significance testing
|
|-- data/
|   |-- reference/
|   |   |-- dm6.fa                    # Drosophila genome (dm6/BDGP6)
|   |   |-- dm6.gtf                   # Gene annotations
|   |   |-- pionx_sites.bed           # PionX element coordinates
|   |   |-- has_sites.bed             # HAS coordinates
|   |   |-- satellite_1688.bed        # 1.688³F satellite coordinates
|   |
|   |-- raw/                          # Raw sequencing data (user-provided)
|   |-- processed/                    # Processed outputs (generated)
|
|-- config/
|   |-- sample_metadata.csv           # Sample information and file paths
|   |-- color_scheme.json             # Figure color palette
```

---

## System Requirements

### Hardware
- **RAM**: Minimum 32 GB (64 GB recommended for Hi-C analysis)
- **Storage**: ~200 GB free space for intermediate files
- **CPU**: Multi-core processor recommended (HiC-Pro supports parallelization)

### Software Dependencies

#### Core Tools (must be installed separately)
| Tool | Version | Purpose | Installation |
|------|---------|---------|-------------|
| HiC-Pro | 2.11.1 | Hi-C processing | [GitHub](https://github.com/nservant/HiC-Pro) |
| HiCExplorer | 3.4.3 | TAD/compartment calling | `conda install hicexplorer` |
| Mustache | latest | Loop calling | `pip install mustache-dnase` |
| MACS2 | 2.2.7.1 | ChIP-seq peak calling | `pip install MACS2` |
| deepTools | 3.5.0 | ChIP aggregation | `pip install deeptools` |
| Bowtie2 | 2.4.4 | Read alignment | `conda install bowtie2` |
| SAMtools | 1.15 | BAM manipulation | `conda install samtools` |
| bedtools | 2.30.0 | BED file operations | `conda install bedtools` |
| pyGenomeTracks | 3.6 | Browser track plotting | `pip install pyGenomeTracks` |

#### R Packages
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "reshape2", "pheatmap", 
                   "VennDiagram", "igraph", "ggraph", "clusterProfiler", 
                   "org.Dm.eg.db", "DESeq2", "GenomicRanges", "rtracklayer"))
```

#### Python Packages
```bash
pip install -r requirements.txt
```

`requirements.txt` contents:
```
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0
seaborn>=0.11.0
h5py>=3.0.0
cooler>=0.8.0
biopython>=1.79
pybedtools>=0.8.0
```

---

## Installation

### Option 1: Conda Environment (Recommended)

```bash
# Clone repository
git clone https://github.com/ankushsawant/chromatin_3D.git
cd chromatin_3D

# Create conda environment
conda env create -f environment.yml
conda activate chromatin_3d

# Install additional R packages (run in R)
Rscript scripts/install_r_packages.R
```

### Option 2: Manual Installation

```bash
# Clone repository
git clone https://github.com/ankushsawant/chromatin_3D.git
cd chromatin_3D

# Install Python dependencies
pip install -r requirements.txt

# Install R dependencies
Rscript scripts/install_r_packages.R
```

---

## Data Requirements

### Input Data (User Must Provide)

| Data Type | Format | Source |
|-----------|--------|--------|
| **RNA-seq** | Paired-end FASTQ | Generated in this study or GEO |
| **ChIP-seq** | Paired-end FASTQ | Generated in this study or GEO |
| **Hi-C** | Paired-end FASTQ | Generated in this study or GEO |
| **Genome** | FASTA | dm6/BDGP6 from FlyBase |
| **Annotation** | GTF/GFF3 | FlyBase release 6 |

### Published Datasets Used

| Dataset | GEO Accession | Description |
|---------|--------------|-------------|
| Embryo Hi-C | GSE15292 | Stage 5, 8, 16 embryos (Hou et al., 2012) |
| Larva Hi-C | GSE220639 | Third instar larvae (Szabo et al., 2023) |
| Head/Ovary ChIP | GSE44210 | Adult head and ovary (Brown et al., 2014) |
| Testis ChIP | GSE15292 | Adult testis (Gan et al., 2010) |

New data generated in this study: **GEO GSE260170**

---

## Usage

### Step 1: Pre-process Raw Sequencing Data

#### Hi-C Processing
```bash
# Configure HiC-Pro
# Edit scripts/01_preprocess/hicpro_config.conf with your paths

# Run HiC-Pro
bash scripts/01_preprocess/run_hicpro.sh \
    --input data/raw/hic/ \
    --output data/processed/hic/ \
    --config scripts/01_preprocess/hicpro_config.conf
```

**Output**: `.matrix`, `.bed`, and `.cool` files for each sample

#### ChIP-seq Processing
```bash
bash scripts/01_preprocess/chipseq_pipeline.sh \
    --input data/raw/chipseq/ \
    --output data/processed/chipseq/ \
    --genome data/reference/dm6.fa \
    --annotation data/reference/dm6.gtf
```

**Output**: Aligned BAMs, peak files (`.narrowPeak`), and bigWig coverage files

#### RNA-seq Processing
```bash
bash scripts/01_preprocess/rnaseq_pipeline.sh \
    --input data/raw/rnaseq/ \
    --output data/processed/rnaseq/ \
    --genome data/reference/dm6.fa \
    --annotation data/reference/dm6.gtf
```

**Output**: Gene-level TPM/FPKM counts and DESeq2 normalized counts

---

### Step 2: Classify Dosage-Compensated Genes

```bash
Rscript scripts/02_dosage_compensation/classify_dc_genes.R \
    --chip data/processed/chipseq/ \
    --rnaseq data/processed/rnaseq/ \
    --output data/processed/dc_classification/ \
    --cutoff auto
```

**Parameters:**
- `--chip`: Directory containing ChIP-seq bigWig or bedGraph files
- `--rnaseq`: Directory containing RNA-seq count matrices
- `--output`: Output directory for classification results
- `--cutoff`: ChIP enrichment cutoff for X vs. autosome separation (default: auto)

**Output files:**
- `cdc_genes.txt` — Constitutively dosage-compensated genes
- `fdc_genes.txt` — Facultatively dosage-compensated genes
- `sdc_genes.txt` — Specifically dosage-compensated genes
- `non_dc_genes.txt` — Non-dosage-compensated genes

---

### Step 3: Hi-C Analysis — Compartments and TADs

#### Call A/B Compartments
```bash
python scripts/03_hic_analysis/call_compartments.py \
    --cool data/processed/hic/sample.cool \
    --resolution 25000 \
    --output data/processed/hic/compartments/
```

#### Call TAD Boundaries
```bash
python scripts/03_hic_analysis/call_tads.py \
    --cool data/processed/hic/sample.cool \
    --resolution 5000 \
    --output data/processed/hic/tads/
```

#### Quantify LRCs and SRCs
```bash
python scripts/03_hic_analysis/quantify_lrcs.py \
    --cool data/processed/hic/sample.cool \
    --resolution 25000 \
    --lrc-threshold 1000000 \
    --output data/processed/hic/lrcs/
```

**Parameters:**
- `--lrc-threshold`: Distance threshold (bp) defining long-range contacts (default: 1,000,000 = 1 Mb)
- Contacts ≤ 1 Mb = SRCs; Contacts > 1 Mb = LRCs

#### Call Loops with Mustache
```bash
python scripts/03_hic_analysis/call_loops_mustache.py \
    --cool data/processed/hic/sample.cool \
    --resolution 5000 \
    --output data/processed/hic/loops/
```

---

### Step 4: ChIP-seq Aggregation and Correlation

#### Anchor-Centered ChIP Aggregation
```bash
bash scripts/04_chipseq_analysis/deeptools_aggregation.sh \
    --peaks data/processed/hic/loops/loop_anchors.bed \
    --signal data/processed/chipseq/sample.bw \
    --output data/processed/chipseq/aggregation/ \
    --before 5000 \
    --after 5000
```

#### ChIP-RNA Correlation
```bash
Rscript scripts/04_chipseq_analysis/correlate_chip_expression.R \
    --chip data/processed/chipseq/genebody_enrichment.txt \
    --rna data/processed/rnaseq/tpm_counts.txt \
    --output data/processed/correlation/
```

---

### Step 5: Generate Figures

```bash
# Main figures
Rscript scripts/05_figures/generate_figures.R \
    --input data/processed/ \
    --output figures/ \
    --figure-type main

# Supplementary figures
Rscript scripts/05_figures/generate_figures.R \
    --input data/processed/ \
    --output figures/supplementary/ \
    --figure-type supplementary

# Browser tracks (locus-level views)
python scripts/05_figures/plot_browser_tracks.py \
    --config config/browser_tracks.ini \
    --region chrX:12000000-16600000 \
    --output figures/browser_chrX_locus.png
```

---

## Statistical Analysis

All two sided Wilcoxon rank-sum tests reported in the manuscript use **Benjamini-Hochberg false discovery rate (FDR) correction** applied within each figure/analysis class.

```bash
# Example: FDR correction for Figure 2 LRC comparisons
Rscript scripts/06_statistics/wilcoxon_fdr.R \
    --input data/processed/hic/lrcs/lrc_values.txt \
    --group-by figure \
    --method BH \
    --output data/processed/statistics/
```

---

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| HiC-Pro fails with memory error | Increase `--cpu` and split input files; or use `hicpro_latest` with cluster mode |
| Mustache loop calling returns no loops | Ensure sufficient sequencing depth (>100M valid pairs) and check resolution parameter |
| deepTools `computeMatrix` runs out of memory | Use `--numberOfProcessors` and `--binSize` parameters to reduce memory load |
| R package `org.Dm.eg.db` installation fails | Install via Bioconductor: `BiocManager::install("org.Dm.eg.db")` |
| pyGenomeTracks throws font errors | Install MS fonts: `sudo apt-get install msttcorefonts` (Linux) or use `--font` parameter |

### Getting Help

For bugs or questions, please open an issue on GitHub:  
https://github.com/ankushsawant/chromatin_3D/issues

---

## Citation

If you use this code or data, please cite:

> Lubna Younas1,2,3, Mujahid Ali2,4, Xinpei Zhang5, Huangyi He1,5, Catherine Regnard6, Qi Zhou1,2,5,7,8*. 
> Dosage compensation is regulated by segregated chromatin modules during development in Drosophila
> *Molecular Systems Biology*. [In press]

---

## Contact

Lubna Younas — First Author  
Email: [lubna.ma528@gmail.com]  
GitHub: [LubnaYounas](https://github.com/Lubna-Younas/Dosage-Compensation/edit/Dosage-Compensation)

---

## Acknowledgements

This work was supported by the European Research Council Starting Grant.
