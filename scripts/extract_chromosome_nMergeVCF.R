#!/usr/bin/env Rscript
library("optparse")
# paired_vcf_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"
#paired_vcf_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_annotated_RefSeq_canonical_paired_noPadding/"

#chr <- 2
#output_dir <- paste0("/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_annotated_RefSeq_canonical_paired_noPadding/vcf_chr", chr, "/")

##
# Rscript scripts/extract_chromosome_nMergeVCF.R -c 1 -i /mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_annotated_RefSeq_canonical_paired_noPadding/ -o /mnt/cloudbiodata_nfs_2/users/hhuang/ > /mnt/cloudbiodata_nfs_2/users/hhuang/README_chr1
#

option_list = list(
  make_option(c("-c", "--chr"), type="numeric", default=NULL, 
              help="chromosome number", metavar="numeric"),
  make_option(c("-i", "--input_dir"), type="character", default=NULL, 
              help="chromosome number", metavar="character"),
  make_option(c("-o", "--output_dir"), type="character", default="out.txt", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

paired_vcf_dir <- opt$input_dir
output_dir <- opt$output_dir
chr <- opt$chr

############

all_vcf_files <- list.files(paired_vcf_dir, pattern = "\\.vcf.gz$")

dir.create(file.path(output_dir), showWarnings = FALSE)

num_files <- length(all_vcf_files)

for(id in 1:num_files){
  ptm <- proc.time()
  
  cat("Paired file #", id, "\n")
  out.filename <- gsub(".vcf.gz", paste0("_chr", chr, ".vcf.gz"), all_vcf_files[id])
  system(paste0("tabix -h ", paired_vcf_dir, all_vcf_files[id], " chr", chr," | bgzip > ", output_dir, out.filename))
  system(paste0("cd ", output_dir, "; tabix -p vcf ", out.filename))
  cat(paste0("tabix -h ", paired_vcf_dir, all_vcf_files[id], " chr", chr," | bgzip > ", output_dir, out.filename), "\n")
  cat(paste0("cd ", output_dir, "; tabix -p vcf ", out.filename), "\n")
  print(proc.time()-ptm)
}
cat("Extraction done! \n")
all_vcf_files <- list.files(output_dir, pattern = "\\.vcf.gz$")
file_list <- paste(all_vcf_files, collapse = " ")

ptm <- proc.time()
system(paste0("cd ", output_dir, "; vcf-merge ", file_list, " | bgzip -c > all_chr", chr,".vcf.gz"))
system(paste0("cd ", output_dir, "; tabix -p vcf all_chr", chr,".vcf.gz"))
cat(paste0("cd ", output_dir, "; vcf-merge ", file_list, " | bgzip -c > all_chr", chr,".vcf.gz"), "\n")
cat(paste0("cd ", output_dir, "; tabix -p vcf all_chr", chr,".vcf.gz"), "\n")
print(proc.time() - ptm)
