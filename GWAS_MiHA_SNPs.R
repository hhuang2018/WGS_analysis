source("util.R")
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(ggplot2)
GWAS_samp_table_fp <- "../ClinVar/GWASH/Metadata/newfnldataib1203.csv"
GWAS_file_fp <- "../ClinVar/GWASH/MIHA/data/"

GWAS_sample_table <-read.csv(GWAS_samp_table_fp, header = T)


MiHA_genotypes <- read.table(paste0(GWAS_file_fp, "MiHAgeno.ped"))

num_SNPs <- dim(MiHA_genotypes)[2]
num_samples <- dim(MiHA_genotypes)[1]

