# list .vcf.gz files
original_file_dir = "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/hli_vcf_annotated/"
destination_dir = "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/hli_vcf_StatsSummary/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
out.filenames <- gsub(".vcf.gz", "", filenames)

num_files <- length(filenames)

for(id in 1:num_files){
  system(paste0("vcf-stats ", original_file_dir, filenames[id], " > ", destination_dir, out.filenames[id], "_stats.txt"))
}