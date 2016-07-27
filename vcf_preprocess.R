require(vcfR)

# cloud
VCF_file_dir <- "/mnt/scratch/hhuang/hli_vcf_renamed/"
output_dir <- "/mnt/scratch/hhuang/hli_vcf_preprocessed/"

# # local
# VCF_file_dir <- "../HLI_VCF_files/
# output_dir <- "../HLI_VCF_preprocessed/"

all_files <- list.files(VCF_file_dir, pattern = "\\.vcf.gz$")

num_files <- length(all_files)
for(id in 1:num_files){
  # donor's VCF - Chromosome 
  vcf_file <- paste0(VCF_file_dir, all_files[id])
  vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  total_num <- dim(vcf_info@fix)[1]
  rm_id_noPass <- which(vcf_info@fix[, 7] != "PASS")
  rm_id_noAlle <- which(!grepl("/", vcf_info@gt[, 2]))
  
  rm_id <- unique(c(rm_id_noPass, rm_id_noAlle))
  
  new_vcf_info <- vcf_info
  new_vcf_info@fix <- new_vcf_info@fix[-rm_id, ]
  new_vcf_info@gt <- new_vcf_info@gt[-rm_id,]
  
  write.vcf(new_vcf_info, file = paste0(output_dir, gsub(".vcf.gz", "_preprocess.vcf.gz", all_files[id])))
  
  cat(all_files[id], ": Removed ", length(rm_id), "(", round(length(rm_id)/total_num, digits = 4)*100, "%) (total); ", length(rm_id_noPass), " (LowQ); ", length(rm_id_noAlle), " (missing allele) \n")
}
