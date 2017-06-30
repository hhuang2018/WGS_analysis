#####
####
#### "java -jar /home/hhuang/tools/snpEff/snpEff.jar 
####      -canon -dataDir /mnt/common/data/reference/grch38/SnpEff_database/ 
####      -csvStats /mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_annotated_RefSeq_geneStats/test.csv 
####      -v GRCh38.p7.RefSeq /mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_renamed/a_100_R_6394954.vcf.gz > test_a_100_R.vcf"


# list .vcf.gz files
original_file_dir = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_renamed/"
destination_dir = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_annotated_RefSeq_geneStats/"
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
  system(paste0("java -jar ", snpEff_dir, "snpEff.jar -canon -dataDir ", snpEff_database,
                      " -csvStats ", destination_dir, filenames[Not_annotated_index[id]], "_summary.csv",
                      " -v GRCh38.p7.RefSeq ", 
                      # destination_dir, filenames[Not_annotated_index[id]], "_snpEff_RefSeq_summary.csv -v ",
                      original_file_dir, filenames[Not_annotated_index[id]], ".vcf.gz "),
               intern = TRUE)
  proc.time() - ptm
  
  # system(paste0("mv snpEff_genes.txt ",  destination_dir, filenames[Not_annotated_index[id]], "_snpEff_RefSeq_genes.txt"))
  # system(paste0("mv snpEff_summary.html ",  destination_dir, filenames[Not_annotated_index[id]], "_snpEff_RefSeq_summary.html"))
  # save(t1, file = paste0(destination_dir, filenames[Not_annotated_index[id]],"_RefSeq_annotated.out_canon.RData"))
}
