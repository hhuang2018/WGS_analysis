source('util.R', echo = FALSE)
library(vcfR)

Chrom <- "chr6"

#VCF_file_dir <- "/mnt/scratch/hhuang/hli_vcf_renamed/"
VCF_file_dir <- "../Analysis/"
# donor's VCF - Chromosome 
vcf_file_D <- paste0(VCF_file_dir, "n_197_D_annotated.vcf.gz")
vcf_D <- read.vcfR(vcf_file_D, verbose = FALSE)

D_chr_variants <- as.data.frame(vcf_D@fix[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
D_chr_variants$POS <- as.integer(D_chr_variants$POS)

D_chr_gt <- as.data.frame(vcf_D@gt[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
D_chr_gt <- D_chr_gt[(D_chr_variants$FILTER == "PASS"), ]
D_chr_variants <- D_chr_variants[(D_chr_variants$FILTER == "PASS"), ] # only look at high-quality variants

D_chr_meta <- vcf_D@meta

num_meta <- length(D_chr_meta)
format_id <- which(grepl("FORMAT=", D_chr_meta))
info_id <- which(grepl("INFO=", D_chr_meta))

# Format
strsplit(gsub("##FORMAT=", "",D_chr_meta[format_id]), ",")

# Annotation

annotation_info <- parse_meta_info(D_chr_meta[info_id], D_chr_variants$INFO)


# #D_chr_meta[info_id]
# ANN_id <- which(grepl("ANN", D_chr_meta[info_id]))
# ann_colnames <- gsub(" ", "", unlist(strsplit(unlist(strsplit(D_chr_meta[info_id[ANN_id]],"'"))[2], "\\|")))
# 
# aa <- sapply(1:length(D_chr_variants$INFO), 
#              function(x) unlist(strsplit(unlist(strsplit(D_chr_variants$INFO[x], split = "ANN="))[2], "\\|")))
# ab <- as.data.frame(aa)
# #
