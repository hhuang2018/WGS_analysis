source('util.R', echo = FALSE)

load("../Data/GRCh38_gene_list.RData")

library(vcfR)

VCF_file_dir <- "/mnt/scratch/hhuang/hli_vcf_annotated/"
Donor_file <- "n_197_D_62284237_annotated.vcf.gz"
Recipient_file <- "n_197_R_1005928_annotated.vcf.gz"
output_dir <- "/mnt/scratch/hhuang/hli_vcf_annotation_stats/"
# donor's VCF - Chromosome 
vcf_file_D <- paste0(VCF_file_dir, Donor_file)
vcf_D <- read.vcfR(vcf_file_D, verbose = FALSE)

# recipient's VCF - Chromosome 
vcf_file_R <- paste0(VCF_file_dir, Recipient_file)
vcf_R <- read.vcfR(vcf_file_R, verbose = FALSE)

output_filename <- paste0(output_dir, unlist(strsplit(Donor_file, "_"))[1], "_",
                          unlist(strsplit(Donor_file, "_"))[2], "_",
                          "R_D_annotated_variants_")
# pdf(paste0("../Output/",output_filename,".pdf"))
D_chr_meta <- vcf_D@meta
R_chr_meta <- vcf_R@meta

for(chr in 1:22){
  cat("Chromosome ", chr, ": ")
  
  Chrom <- paste0("chr", chr)
  
  # donor's VCF - Chromosome 
  D_chr <- as.data.frame(vcf_D@fix[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  D_chr$POS <- as.integer(D_chr$POS)
  
  D_chr_gt <- as.data.frame(vcf_D@gt[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  D_chr_gt <- D_chr_gt[(D_chr_variants$FILTER == "PASS"), ]
  D_chr_variants <- D_chr_variants[(D_chr_variants$FILTER == "PASS"), ] # only look at high-quality variants
  
  D_info_id <- which(grepl("INFO=", D_chr_meta))
  
  D_annotation_info <- parse_meta_info(D_chr_meta[D_info_id], D_chr_variants$INFO)
  # rm(vcf_D)
  
  # recipient's VCF - Chromosome 
  R_chr <- as.data.frame(vcf_R@fix[vcf_R@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  R_chr$POS <- as.integer(R_chr$POS)
  
  R_chr_gt <- as.data.frame(vcf_R@gt[vcf_R@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  R_chr_gt <- R_chr_gt[(R_chr_variants$FILTER == "PASS"), ]
  R_chr_variants <- R_chr_variants[(R_chr_variants$FILTER == "PASS"), ] # only look at high-quality variants
  
  R_info_id <- which(grepl("INFO=", R_chr_meta))
  
  R_annotation_info <- parse_meta_info(R_chr_meta[R_info_id], R_chr_variants$INFO)
  # rm(vcf_R)
  
  #

  save(D_annotation_info, R_annotation_info, file = paste0(output_filename, "_chr", chr, ".RData"))
  
#   plot(diff_num$total, 
#        main = paste0("Full Gene region (Chormosome", chr,")"),
#        xlab = "Gene Index",
#        ylab = "Number of different variants")
#   plot(diff_num$intron, 
#        main = paste0("Gene Intron regions (Chormosome", chr,")"),
#        xlab = "Gene Index",
#        ylab = "Number of different variants")
#   plot(diff_num$exon, 
#        main = paste0("Gene Exon regions (Chormosome", chr,")"),
#        xlab = "Gene Index",
#        ylab = "Number of different variants")
#   
#   save(diff_num, file = paste0("../Output/",output_filename, "_chr", chr,".RData"))
#   
  cat("Done!\n")
}

###########
rm(vcf_D)
rm(vcf_R)
VCF_file_dir <- "/mnt/scratch/hhuang/hli_vcf_annotated/"
Donor_file <- "a_208_D_528085277_annotated.vcf.gz"
Recipient_file <- "a_208_R_1018111_annotated.vcf.gz"
output_dir <- "/mnt/scratch/hhuang/hli_vcf_annotation_stats/"
# donor's VCF - Chromosome 
vcf_file_D <- paste0(VCF_file_dir, Donor_file)
vcf_D <- read.vcfR(vcf_file_D, verbose = FALSE)

# recipient's VCF - Chromosome 
vcf_file_R <- paste0(VCF_file_dir, Recipient_file)
vcf_R <- read.vcfR(vcf_file_R, verbose = FALSE)

output_filename <- paste0(output_dir, unlist(strsplit(Donor_file, "_"))[1], "_",
                          unlist(strsplit(Donor_file, "_"))[2], "_",
                          "R_D_annotated_variants_")
# pdf(paste0("../Output/",output_filename,".pdf"))
D_chr_meta <- vcf_D@meta
R_chr_meta <- vcf_R@meta

for(chr in 1:22){
  cat("Chromosome ", chr, ": ")
  
  Chrom <- paste0("chr", chr)
  
  # donor's VCF - Chromosome 
  D_chr <- as.data.frame(vcf_D@fix[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  D_chr$POS <- as.integer(D_chr$POS)
  
  D_chr_gt <- as.data.frame(vcf_D@gt[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  D_chr_gt <- D_chr_gt[(D_chr_variants$FILTER == "PASS"), ]
  D_chr_variants <- D_chr_variants[(D_chr_variants$FILTER == "PASS"), ] # only look at high-quality variants
  
  D_info_id <- which(grepl("INFO=", D_chr_meta))
  
  D_annotation_info <- parse_meta_info(D_chr_meta[D_info_id], D_chr_variants$INFO)
  # rm(vcf_D)
  
  # recipient's VCF - Chromosome 
  R_chr <- as.data.frame(vcf_R@fix[vcf_R@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  R_chr$POS <- as.integer(R_chr$POS)
  
  R_chr_gt <- as.data.frame(vcf_R@gt[vcf_R@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  R_chr_gt <- R_chr_gt[(R_chr_variants$FILTER == "PASS"), ]
  R_chr_variants <- R_chr_variants[(R_chr_variants$FILTER == "PASS"), ] # only look at high-quality variants
  
  R_info_id <- which(grepl("INFO=", R_chr_meta))
  
  R_annotation_info <- parse_meta_info(R_chr_meta[R_info_id], R_chr_variants$INFO)
  # rm(vcf_R)
  
  #
  
  save(D_annotation_info, R_annotation_info, file = paste0(output_filename, "_chr", chr, ".RData"))
  
  #   plot(diff_num$total, 
  #        main = paste0("Full Gene region (Chormosome", chr,")"),
  #        xlab = "Gene Index",
  #        ylab = "Number of different variants")
  #   plot(diff_num$intron, 
  #        main = paste0("Gene Intron regions (Chormosome", chr,")"),
  #        xlab = "Gene Index",
  #        ylab = "Number of different variants")
  #   plot(diff_num$exon, 
  #        main = paste0("Gene Exon regions (Chormosome", chr,")"),
  #        xlab = "Gene Index",
  #        ylab = "Number of different variants")
  #   
  #   save(diff_num, file = paste0("../Output/",output_filename, "_chr", chr,".RData"))
  #   
  cat("Done!\n")
}
