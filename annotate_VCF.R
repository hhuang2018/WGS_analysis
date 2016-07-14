# list .vcf.gz files
original_file_dir = "/mnt/scratch/hhuang/hli_vcf_renamed/"
destination_dir = "/mnt/scratch/hhuang/hli_vcf_annotated/"
snpEff_dir = "/home/hhuang/tools/snpEff/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
filenames <- gsub(".vcf.gz", "", filenames)

annotated_filenames <- list.files(destination_dir, pattern = "\\.vcf.gz$")
annotated_filenames <- gsub("_annotated.vcf.gz", "", annotated_filenames)

Not_annotated_index <- which(!(filenames %in% annotated_filenames))

num_files <- length(Not_annotated_index)

for(id in 1:num_files){
  t1 <- system(paste0("java -Xmx4g -jar ", snpEff_dir, "snpEff.jar -v GRCh38.82 ",
                      original_file_dir, filenames[Not_annotated_index[id]], ".vcf.gz > ", 
                      destination_dir, filenames[Not_annotated_index[id]],"_annotated.vcf"), 
               intern = TRUE)
  system(paste0("bgzip ", destination_dir, filenames[Not_annotated_index[id]],"_annotated.vcf"))
  system(paste0("tabix -p vcf ", destination_dir, filenames[Not_annotated_index[id]],"_annotated.vcf.gz"))
  system(paste0("mv snpEff_genes.txt ",  destination_dir, filenames[Not_annotated_index[id]], "_snpEff_genes.txt"))
  system(paste0("mv snpEff_summary.html ",  destination_dir, filenames[Not_annotated_index[id]], "_snpEff_summary.html"))
  save(t1, file = paste0(destination_dir, filenames[Not_annotated_index[id]],"_annotated.out.RData"))
}