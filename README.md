# SNP Frequency Annotation & Visualization Pipeline

## 📊 Overview

This repository provides an R-based pipeline to analyze and visualize **SNP allele frequency data** from a VCF file, with a focus on **chromosome 22**. It compares **global allele frequencies (AF)** with those in the **South Asian (SAS)** population, annotates nearby genes using **Ensembl**, and generates plots including:

- Allele frequency comparison (Global vs SAS)
- SNP density
- Gene-wise SNP counts
- Manhattan-style −log10(AF) plots
- Delta AF visualization
- Heatmaps of AF and SAS_AF

The dataset used is from the **1000 Genomes Project** (Phase 3, chromosome 22).

  

## 🧬 Input (User-provided)

- VCF file:  
  Example file path in the script:
  ```
  ~/my_1000g_data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
  ```

> **Note**: Only the first 500 SNPs are processed for demonstration.

  

## ⚙️ Requirements

Install the following R packages:

```r
install.packages("vcfR")
install.packages("ggplot2")
install.packages("dplyr")
install.packages("reshape2")
install.packages("GenomicRanges")
BiocManager::install("biomaRt")
```

  

## 🚀 How to Run

1. Clone the repository:

```bash
git clone https://github.com/Dhrupad210/SNP-frequency-annotation-visualization.git
cd SNP-frequency-annotation-visualization
```

2. Open and edit the script `snp_frequency_pipeline.R`:
   - Update the path to your VCF file (`vcf_path` variable).

3. Run the script in **R or RStudio**:
   ```r
   source("snp_frequency_pipeline.R")
   ```

4. The script will:
   - Extract SNPs and allele frequencies
   - Filter SNPs with `SAS_AF > 0.05`
   - Annotate nearby genes using Ensembl (±100 kb)
   - Generate plots as `.png` files
   - Save a CSV file with merged SNP-gene data

  

## 📁 Repository Structure

```
SNP-frequency-annotation-visualization/
├── snp_frequency_pipeline.R   # Main R script
└── README.md                  # Documentation
```

> Note: Output CSV and PNG files are generated at runtime but not included in this repository.

  

## 🧠 Key Biological Concepts

- **AF**: Global allele frequency across all populations.
- **SAS_AF**: Frequency specific to the South Asian population.
- **Delta AF (AF - SAS_AF)**: Highlights population-specific allele shifts.
- **Nearby Genes**: Identified using Ensembl within ±100,000 bp of each SNP using `biomaRt`.

## 📈 Output Visualizations

| Plot | Description |
|------|-------------|
| `RawAF_vs_SAS_AF_scatter.png` | Global vs South Asian Allele Frequency |
| `SNP_density_chr22.png`       | SNP density along chromosome 22 |
| `Gene_SNP_count_barplot.png`  | Genes with ≥2 nearby SNPs |
| `Manhattan_style_AF_plot.png` | −log10(AF) Manhattan-style plot |
| `delta_AF_barplot.png`        | Difference between Global AF and SAS_AF |
| `heatmap_AF.png`              | Heatmap of AF and SAS_AF across SNPs |


## 👤 Author

**Dhrupad Banerjee**  
GitHub: [@Dhrupad210](https://github.com/Dhrupad210)

  

## 📜 License

This project is open for academic and educational purposes.
