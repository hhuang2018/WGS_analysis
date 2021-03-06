# list .vcf.gz files
original_file_dir = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_renamed/"
destination_dir = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_annotated_RefSeq/"
snpEff_dir = "/home/hhuang/tools/snpEff/"
snpEff_database = "/mnt/common/data/reference/grch38/SnpEff_database/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
filenames <- gsub(".vcf.gz", "", filenames)

annotated_filenames <- list.files(destination_dir, pattern = "\\.vcf.gz$")
annotated_filenames <- gsub("_annotated.vcf.gz", "", annotated_filenames)

Not_annotated_index <- which(!(filenames %in% annotated_filenames))

num_files <- length(Not_annotated_index)

for(id in 1:num_files){
  ptm <- proc.time()
  t1 <- system(paste0("java -Xmx32g -jar ", snpEff_dir, "snpEff.jar -t -noStats -dataDir ", snpEff_database,
                      " -v GRCh38.p7.RefSeq ", 
                      # destination_dir, filenames[Not_annotated_index[id]], "_snpEff_RefSeq_summary.csv -v ",
                       original_file_dir, filenames[Not_annotated_index[id]], ".vcf.gz > ",
                      destination_dir, filenames[Not_annotated_index[id]],"_RefSeq_annotated.vcf"),
               intern = TRUE)
  proc.time() - ptm
  
  system(paste0("bgzip ", destination_dir, filenames[Not_annotated_index[id]],"_RefSeq_annotated.vcf"))
  system(paste0("tabix -p vcf ", destination_dir, filenames[Not_annotated_index[id]],"_RefSeq_annotated.vcf.gz"))
  # system(paste0("mv snpEff_genes.txt ",  destination_dir, filenames[Not_annotated_index[id]], "_snpEff_RefSeq_genes.txt"))
  # system(paste0("mv snpEff_summary.html ",  destination_dir, filenames[Not_annotated_index[id]], "_snpEff_RefSeq_summary.html"))
  save(t1, file = paste0(destination_dir, filenames[Not_annotated_index[id]],"_RefSeq_annotated.out.RData"))
}

##### need csv format of summary
# for(id in 1:num_files){
#   system(paste0("java -Xmx16g -jar ", snpEff_dir, "snpEff.jar ann -csvStats ",
#                       destination_dir, filenames[Not_annotated_index[id]], "_snpEff_summary.csv -v ",
#                       "GRCh38.82 ", original_file_dir, filenames[Not_annotated_index[id]], ".vcf.gz > ", 
#                       destination_dir, filenames[Not_annotated_index[id]],"_annotated.vcf"), 
#                intern = TRUE)
#   system(paste0("rm ", destination_dir, filenames[Not_annotated_index[id]],"_annotated.vcf"))
#   system(paste0("rm ", destination_dir, filenames[Not_annotated_index[id]], "_snpEff_summary.genes.txt"))
#   # system(paste0("mv snpEff_summary.csv ",  destination_dir, filenames[Not_annotated_index[id]], "_snpEff_summary.csv"))
#   #save(t1, file = paste0(destination_dir, filenames[Not_annotated_index[id]],"_annotated.out.RData"))
# }