
vcf_file_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_vcf_renamed/"

paired_vcf_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_renamed_paired_vcf/"

alternative_paired_vcf_dir <- paired_vcf_dir # "/mnt/scratch/hhuang/paired_vcf/"

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

cat(">>>>>>>>>>>>>>>>>>>>>>>>>>\nParing is Finished! Starting Chr6 Extraction:\n")
########### Extract chromosome 6
# paired_vcf_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"

all_vcf_files <- list.files(paired_vcf_dir, pattern = "\\.vcf.gz$")

chr <- 6
output_dir <- paste0("/mnt/cloudbiodata_nfs_2/users/hhuang/hli_renamed_all_chr", chr, "/")

# dir.create(file.path(output_dir), showWarnings = FALSE)

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

cat("Merging all extracted Chr6 is done!! Start IBD calculation:\n")

########### Calculate IBD for chr6
IBD_out_fp <- "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_renamed_all_chr6/IBD/"

system(paste0("java -Xmx64g -jar ~/tools/ibdseq.r1206.jar gt=", output_dir, "all_chr", chr, ".vcf.gz out=", IBD_out_fp ,"ibdseq_output_all_chr6.gt nthread=16"))
cat(paste0("java -Xmx8g -jar ~/tools/ibdseq.r1206.jar gt=", output_dir, "all_chr", chr, ".vcf.gz out=", IBD_out_fp ,"ibdseq_output_all_chr6.gt nthread=16", "\n"))
