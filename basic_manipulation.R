source('util.R', echo = FALSE)

# read GRCh38 Gene list
# GRCh38_gene_list <- read.table(file = "../UCSCgenome/database/refGene.txt", header = FALSE, sep = "\t")
# colnames(GRCh38_gene_list) <- c("bin", "RefSeqName", "chrom", "strand", "txStart", "txEnd",
#                                 "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
#                                 "id", "GeneName", "cdsStartStat", "cdsEndStat", "exonFrames")
# save(GRCh38_gene_list, file = "../Data/GRCh38_gene_list.RData")
load("../Data/GRCh38_gene_list.RData")

library(vcfR)

Chrom <- "chr6"

VCF_file_dir <- "../Analysis/"

# donor's VCF - Chromosome 
vcf_file_D <- paste0(VCF_file_dir, "n_197_D_annotated.vcf.gz")
vcf_D <- read.vcfR(vcf_file_D, verbose = FALSE)

D_chr <- as.data.frame(vcf_D@fix[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
D_chr$POS <- as.integer(D_chr$POS)

rm(vcf_D)

# recipient's VCF - Chromosome 
vcf_file_R <- paste0(VCF_file_dir, "n_197_R_annotated.vcf.gz")
vcf_R <- read.vcfR(vcf_file_R, verbose = FALSE)

R_chr <- as.data.frame(vcf_R@fix[vcf_R@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
R_chr$POS <- as.integer(R_chr$POS)

rm(vcf_R)

# Check HLA genes
GeneName <- "HLA-DQB1"
geneInfo <- GRCh38_gene_list[which(GRCh38_gene_list$GeneName == GeneName)[1], ]

D_stats <- gene_variant_stats(geneInfo, D_chr)
R_stats <- gene_variant_stats(geneInfo, R_chr)
