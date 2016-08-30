
paired_vcf_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"

all_vcf_files <- list.files(paired_vcf_dir, pattern = "\\.vcf.gz$")

chr <- 4
output_dir <- paste0("/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/vcf_chr", chr, "/")

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
system(paste0("cd ", output_dir, "; vcf-merge -R 0/0 ", file_list, " | bgzip -c > all_chr", chr,".vcf.gz"))
system(paste0("cd ", output_dir, "; tabix -p vcf all_chr", chr,".vcf.gz"))
cat(paste0("cd ", output_dir, "; vcf-merge -R 0/0 ", file_list, " | bgzip -c > all_chr", chr,".vcf.gz"), "\n")
cat(paste0("cd ", output_dir, "; tabix -p vcf all_chr", chr,".vcf.gz"), "\n")
print(proc.time() - ptm)
