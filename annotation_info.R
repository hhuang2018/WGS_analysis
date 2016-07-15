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

D_chr_meta <- vcf_D@meta
D_chr_gt <- as.data.frame(vcf_D@gt[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
D_chr_gt <- D_chr_gt[(D_chr_variants$FILTER == "PASS"), ]
D_chr_variants <- D_chr_variants[(D_chr_variants$FILTER == "PASS"), ] # only look at high-quality variants

# rm(vcf_D)