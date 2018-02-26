vcf_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/hli_vcf_annotated_preprocessed/"

all_vcf_files <- list.files(vcf_file_dir, pattern = "\\.vcf.gz$")
indexed_fies <- list.files(vcf_file_dir, pattern = "\\.vcf.gz.tbi$")
indexed_fies <- gsub(".tbi", "", indexed_fies)

num_files <- length(all_vcf_files)

for(id in 1:num_files){
  if(!is.element(all_vcf_files[id], indexed_fies)){
    vcf_name <- gsub(".gz", "", all_vcf_files[id])
    system(paste0("cd ", vcf_file_dir,"; gunzip ", all_vcf_files[id]))
    system(paste0("cd ", vcf_file_dir,"; bgzip ", vcf_file_dir, vcf_name))
    system(paste0("cd ", vcf_file_dir,"; tabix -p vcf ", vcf_file_dir, all_vcf_files[id]))
    
    cat(all_vcf_files[id], " done! \n")
  }
  
}