vcf_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/hli_vcf_preprocessed/"

all_vcf_files <- list.files(vcf_file_dir, pattern = "\\.vcf.gz$")

num_files <- length(all_vcf_files)

for(id in 1:num_files){
  
  vcf_name <- gsub(".gz", "", all_vcf_files[id])
  eval(parse(text = paste0("gunzip ", vcf_file_dir, all_vcf_files[id])))
  eval(parse(text = paste0("bgzip ", vcf_file_dir, vcf_name)))
  eval(parse(text = paste0("tabix -p vcf ", vcf_file_dir, all_vcf_files[id])))
  
  cat(all_vcf_files[id], " done! \n")
  
}