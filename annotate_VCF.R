# list .vcf.gz files
original_file_dir = "/mnt/cloudbiodata_nfs_1/hli_data/"
destination_dir = "/mnt/scratch/hli_vcf_annotated/"
snpEff_dir = "/home/hhuang/tools/snpEff/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
filenames <- gsub(".vcf.gz", "", filenames)

num_files <- length(filenames)

for(id in 1:num_files){
  t1 <- system(paste0("java -Xmx4g -jar ", snpEff_dir, "snpEff.jar -v GRCh38.82 ",
                      original_file_dir, filenames[id], ".vcf.gz > ", 
                      destination_dir, filenames[id],"_annotated.vcf"), 
               intern = TRUE)
  system(paste0("bgzip ", destination_dir, filenames[id],"_annotated.vcf"))
  system(paste0("tabix -p vcf ", destination_dir, filenames[id],"_annotated.vcf.gz"))
  system(paste0("mv snpEff_genes.txt ",  destination_dir, filenames[id], "_snpEff_genes.txt"))
  system(paste0("mv snpEff_summary.html ",  destination_dir, filenames[id], "_snpEff_summary.html"))
  save(t1, file = paste0(destination_dir, filenames[id],"_annotated.out.RData"))
}