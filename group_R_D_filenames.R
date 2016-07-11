
data_source_dir = "/mnt/cloudbiodata_nfs_1/hli_data/"
data_dest_dir = "/mnt/scratch/hhuang/hli_vcf_renamed/"
  
# HLI_R_D_list <- read.csv("../Data/HLI_ID_pairs.csv")
# colnames(HLI_R_D_list) <- c("Index", "RID", "DID", "BMT")
# 
# HLI_R_D_list$RID <- as.integer(gsub("-", "", HLI_R_D_list$RID))
# HLI_R_D_list$DID <- as.integer(gsub("-", "", HLI_R_D_list$DID))
# 
# save(HLI_R_D_list, file = "../Data/HLI_ID_pairs_table.RData")

load("../Data/HLI_ID_pairs_table.RData")

# HLI_outcome_list1 <- read.csv("../Data/hli_discovery_cohort_list_9feb2015.csv")
# HLI_outcome_list2 <- read.csv("../Data/hli_discvry_and_validation_cohort_list_11feb2015.csv")
# HLI_outcome_list3 <- read.csv("../Data/hli_discvry_and_validation_cohort_list_17feb2015.csv")
# HLI_outcome_list4 <- read.csv("../Data/hli_discvry_and_validation_cohort_list_19feb2015.csv")
# 
# HLI_outcome_list <- rbind(HLI_outcome_list1[, 1:7], 
#                           HLI_outcome_list2[, 1:7],
#                           HLI_outcome_list3[, 1:7],
#                           HLI_outcome_list4[, 1:7])
# aa <- which(duplicated(HLI_outcome_list[, -1]))
# HLI_outcome_list <- HLI_outcome_list[-aa, ]
# HLI_outcome_list <- HLI_outcome_list[-dim(HLI_outcome_list)[1], ]
# save(HLI_outcome_list, file = "../Data/HLI_outcome_table.RData")
# write.csv(HLI_outcome_list, file = "../Data/HLI_outcome.csv", row.names = F)
load("../Data/HLI_outcome_table.RData")

source_filenames <- list.files(data_source_dir, pattern = "\\.vcf.gz$")
source_filenames <- gsub(".vcf.gz", "", source_filenames)
source_filenames2 <- gsub("-", "", source_filenames)

dest_filenames <- list.files(data_dest_dir, pattern = "\\.vcf.gz$")
dest_filenames <- gsub(".vcf.gz", "", dest_filenames)

dest_ids <- sapply(1:length(dest_filenames), function(x) unlist(strsplit(dest_filenames[x], "_"))[4])
unassigned_files_index <- which(!(source_filenames2 %in% dest_ids))
#filenames <- list.files("../HLI_VCF_files/", pattern = "\\.vcf.gz$")
#filenames <- gsub(".vcf.gz", "", filenames)

#unassigned_files_index <- which(!grepl("[[:alpha:]]", filenames))


for(id in 1:length(unassigned_files_index)){
  
  Rindex <- which(HLI_R_D_list$RID %in% as.integer(source_filenames2[unassigned_files_index[id]]))
  if(length(Rindex) == 0){
    Dindex <- which(HLI_R_D_list$DID %in% as.integer(source_filenames2[unassigned_files_index[id]]))
    R_D <- "D"
    Index <- Dindex
    
    outcomeDindex <- which(HLI_outcome_list$did %in% as.integer(source_filenames2[unassigned_files_index[id]]))
    if(length(outcomeDindex) == 1){
      if(HLI_outcome_list$newagvhdgrp[outcomeDindex] == 1) GVHD <- "n" else GVHD <- "a"
    } else GVHD <- "X"
    
    
  }else {
    R_D <- "R"
    Index <- Rindex
    
    outcomeDindex <- which(HLI_outcome_list$rid %in% as.integer(source_filenames2[unassigned_files_index[id]]))
    if(length(outcomeDindex) == 1){
      if(HLI_outcome_list$newagvhdgrp[outcomeDindex] == 1) GVHD <- "n" else GVHD <- "a"
    } else GVHD <- "X"
  }
  
  if(Index !=0){
#     file.rename(from = paste0("/mnt/cloudbiodata_nfs_1/hli_data/",filenames[unassigned_files_index[id]], ".vcf.gz.tbi"), 
#                 to = paste0("/mnt/cloudbiodata_nfs_1/hli_data/", GVHD, "_", Index,"_", R_D, "_", filenames[unassigned_files_index[id]], ".vcf.gz.tbi"))
    
#     file.rename(from = paste0("../HLI_VCF_files/",filenames[unassigned_files_index[id]], ".vcf.gz"), 
#                 to = paste0("../HLI_VCF_files/", GVHD, "_", Index,"_", R_D, "_", filenames[unassigned_files_index[id]], ".vcf.gz"))
    file.copy(from = paste0(data_source_dir, source_filenames[unassigned_files_index[id]], "vcf.gz"),
              to = paste0(data_dest_dir, GVHD, "_", Index, "_", R_D, "_", source_filenames2[unassigned_files_index[id]], ".vcf.gz"))
    
    file.copy(from = paste0(data_source_dir, source_filenames[unassigned_files_index[id]], "vcf.gz.tbi"),
              to = paste0(data_dest_dir, GVHD, "_", Index, "_", R_D, "_", source_filenames2[unassigned_files_index[id]], ".vcf.gz.tbi"))
    
  }
  Index <- 0
}

#### Undo
# filenames <- list.files("../HLI_VCF_files/")
# filenames <- gsub(".vcf.gz", "", filenames)
# 
# assigned_files_index <- which(grepl("[[:alpha:]]", filenames))
# 
# for(id in 1:length(assigned_files_index)){
#   new_filename <- strsplit(filenames[assigned_files_index[id]], "_")[[1]][4]
#   file.rename(from = paste0("../HLI_VCF_files/",filenames[assigned_files_index[id]], ".vcf.gz"), 
#               to = paste0("../HLI_VCF_files/", new_filename, ".vcf.gz"))
#   
# }

#####
# D_R_HLA_typing <- read.csv("../HLI_hla_mg_v1.csv")
# 
# aa <- intersect(D_R_HLA_typing$bmt_case_num, HLI_R_D_list$BMT)
