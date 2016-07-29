vcf_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/hli_vcf_annotated_preprocessed/"

all_vcf_files <- list.files(vcf_file_dir, pattern = "\\.vcf.gz$")

num_files <- length(all_vcf_files)

for(id in 1:num_files){
  
  vcf_name <- gsub(".gz", "", all_vcf_files[id])
  system(paste0("cd ", vcf_file_dir,"; gunzip ", all_vcf_files[id]))
  system(paste0("cd ", vcf_file_dir,"; bgzip ", vcf_file_dir, vcf_name))
  system(paste0("cd ", vcf_file_dir,"; tabix -p vcf ", vcf_file_dir, all_vcf_files[id]))
  
  cat(all_vcf_files[id], " done! \n")
  
}