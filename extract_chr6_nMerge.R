
paired_vcf_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"

all_vcf_files <- list.files(paired_vcf_dir, pattern = "\\.vcf.gz$")

output_dir <- "/mnt/scratch/hhuang/vcf_chr6/"

num_files <- length(all_vcf_files)

for(id in 1:num_files){
  ptm <- proc.time()
  
  cat("Paired file #", id, "\n")
  out.filename <- gsub(".vcf.gz", "_chr6.vcf.gz", all_vcf_files[id])
  system(paste0("tabix -h ", paired_vcf_dir, all_vcf_files[id], " chr6 | bgzip > ", output_dir, out.filename))
  cat(paste0("tabix -h ", paired_vcf_dir, all_vcf_files[id], " chr6 | bgzip > ", output_dir, out.filename), "\n")
  
  print(proc.time()-ptm)
}
cat("Extraction done! \n")
all_vcf_files <- list.files(output_dir, pattern = "\\.vcf.gz$")
file_list <- paste(all_vcf_files, collapse = " ")

ptm <- proc.time()
system(paste0("cd ", output_dir, "; vcf-merge -R 0/0 ", file_list, " | bgzip > all_chr6.vcf.gz"))
cat(paste0("cd ", output_dir, "; vcf-merge -R 0/0 ", file_list, " | bgzip > all_chr6.vcf.gz"), "\n")
print(proc.time() - ptm)
