
vcf_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/hli_vcf_annotated_preprocessed/"

paired_vcf_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"

alternative_paired_vcf_dir <- "/mnt/scratch/hhuang/paired_vcf/"

all_vcf_files <- list.files(vcf_file_dir, pattern = "\\.vcf.gz$")
# indexed_fies <- list.files(vcf_file_dir, pattern = "\\.vcf.gz.tbi$")
# indexed_fies <- gsub(".tbi", "", indexed_fies)
all_groupIDs <- sapply(1:length(all_vcf_files), function(x) paste0(unlist(strsplit(all_vcf_files[x], "_"))[1:2], collapse = "_"))

existed_files <- list.files(paired_vcf_dir, pattern = "\\.vcf.gz$")
groupIDs <- gsub(".vcf.gz", "", existed_files)

remain_ids <- which(!(all_groupIDs %in% groupIDs))
all_groupIDs <- all_groupIDs[remain_ids]
all_vcf_files <- all_vcf_files[remain_ids]

paired_IDs <- which(duplicated(all_groupIDs))

num_files <- length(paired_IDs)

for(id in 1:num_files){
  
  index <- which(all_groupIDs %in% all_groupIDs[paired_IDs[id]])
  
  if(grepl("R", all_groupIDs[index[1]])){ 
    Recipient_file <- all_vcf_files[index[1]]
    Donor_file <- all_vcf_files[index[2]]
  }else{
    Recipient_file <- all_vcf_files[index[2]]
    Donor_file <- all_vcf_files[index[1]]
  }
  
  cat("Group # ", id, "\nCommand: \n")
  
  ptm <- proc.time() 
  system(paste0("cd ", alternative_paired_vcf_dir,"; vcf-merge -R 0/0 ", vcf_file_dir, Recipient_file," ",
                vcf_file_dir, Donor_file, " | bgzip > ", all_groupIDs[paired_IDs[id]], ".vcf.gz"))
  system(paste0("cd ", alternative_paired_vcf_dir, "; tabix -p vcf ", all_groupIDs[paired_IDs[id]], ".vcf.gz"))
  
  cat(paste0("cd ", alternative_paired_vcf_dir,"; vcf-merge -R 0/0 ", vcf_file_dir, Recipient_file," ",
             vcf_file_dir, Donor_file, " | bgzip > ", all_groupIDs[paired_IDs[id]], ".vcf.gz"), "\n")
  cat(paste0("cd ", alternative_paired_vcf_dir, "; tabix -p vcf ", all_groupIDs[paired_IDs[id]], ".vcf.gz"), "\n")
  
  print(proc.time() - ptm)
  
  cat(all_groupIDs[paired_IDs[id]], " done! \n\n")
  
}