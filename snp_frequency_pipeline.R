# snp_frequency_pipeline.R

# Load required libraries
library(vcfR)
library(biomaRt)
library(ggplot2)
library(dplyr)
library(GenomicRanges)
library(reshape2)

# Load VCF data (500 SNPs from chr22)
vcf_path <- "~/my_1000g_data/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
vcf <- read.vcfR(vcf_path, nrows = 500)

# Extract allele frequencies (total and South Asian)
af_total <- sapply(strsplit(extract.info(vcf, "AF"), ","), function(x) as.numeric(x[1]))
sas_af <- sapply(strsplit(extract.info(vcf, "SAS_AF"), ","), function(x) as.numeric(x[1]))

# Create SNP info table
snp_info <- data.frame(
  rsid = vcf@fix[, "ID"],
  chr = vcf@fix[, "CHROM"],
  pos = as.numeric(vcf@fix[, "POS"]),
  AF = af_total,
  SAS_AF = sas_af,
  stringsAsFactors = FALSE
)

# Filter SNPs common in South Asians (SAS_AF > 0.05)
common_SAS <- snp_info %>% filter(SAS_AF > 0.05)

# Annotate nearby genes (±100 kb) using Ensembl
ensembl_gene <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

genes_nearby_all <- do.call(rbind, lapply(1:nrow(common_SAS), function(i) {
  chr <- common_SAS$chr[i]
  pos <- common_SAS$pos[i]
  rsid <- common_SAS$rsid[i]
  
  genes <- getBM(
    attributes = c("external_gene_name", "ensembl_gene_id", 
                   "start_position", "end_position", "chromosome_name"),
    filters = c("chromosome_name", "start", "end"),
    values = list(chr, pos - 100000, pos + 100000),
    mart = ensembl_gene
  )
  
  if (nrow(genes) > 0) {
    genes$rsid <- rsid
    return(genes)
  } else {
    return(NULL)
  }
}))

# Merge SNP and gene info
merged_info <- merge(common_SAS, genes_nearby_all, by = "rsid")
write.csv(merged_info, "SNPs_with_nearby_genes_SAS_AF_gt_0.05.csv", row.names = FALSE)

# Visualization 1: Global AF vs SAS_AF scatter plot
png("RawAF_vs_SAS_AF_scatter.png", width = 800, height = 600)
ggplot(merged_info, aes(x = AF, y = SAS_AF)) +
  geom_point(alpha = 0.7, color = "darkred") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue") +
  labs(title = "Global vs South Asian Allele Frequency", x = "Global AF", y = "SAS AF") +
  theme_minimal()
dev.off()

# Visualization 2: SNP density plot
png("SNP_density_chr22.png", width = 800, height = 600)
ggplot(merged_info, aes(x = pos)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  labs(title = "SNP Density on Chromosome 22", x = "Genomic Position", y = "Density") +
  theme_minimal()
dev.off()

# Visualization 3: Barplot of genes with ≥2 nearby SNPs
top_genes <- merged_info %>%
  count(external_gene_name, sort = TRUE) %>%
  filter(!is.na(external_gene_name) & external_gene_name != "", n >= 2)

png("Gene_SNP_count_barplot.png", width = 800, height = 600)
ggplot(top_genes, aes(x = reorder(external_gene_name, n), y = n)) +
  geom_bar(stat = "identity", fill = "forestgreen") +
  coord_flip() +
  labs(title = "Genes with ≥2 SNPs (SAS_AF > 0.05)", x = "Gene", y = "SNP Count") +
  theme_minimal()
dev.off()

# Visualization 4: Manhattan-style plot of −log10(AF)
merged_info$log_AF <- -log10(merged_info$AF)

png("Manhattan_style_AF_plot.png", width = 800, height = 600)
ggplot(merged_info, aes(x = pos, y = log_AF)) +
  geom_point(color = "purple", alpha = 0.6) +
  labs(title = "SNP Manhattan-style Plot (−log10(AF))", 
       x = "Position on chr22", y = "−log10(AF)") +
  theme_minimal()
dev.off()

# Visualization 5: Delta AF (AF - SAS_AF) bar plot
merged_info$delta_AF <- merged_info$AF - merged_info$SAS_AF

ggplot(merged_info, aes(x = pos, y = delta_AF)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Global vs SAS Delta AF", x = "Position", y = "ΔAF (AF − SAS_AF)") +
  theme_classic()
ggsave("delta_AF_barplot.png", width = 8, height = 5, dpi = 300)

# Visualization 6: Heatmap of AF and SAS_AF
melted <- melt(merged_info[, c("pos", "AF", "SAS_AF")], id.vars = "pos")

ggplot(melted, aes(x = variable, y = as.factor(pos), fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Heatmap of Allele Frequencies", x = "Frequency Type", y = "SNP Position") +
  theme_minimal()
ggsave("heatmap_AF.png", width = 6, height = 7, dpi = 300)

# Save session history
savehistory("my_r_session_commands.Rhistory")
