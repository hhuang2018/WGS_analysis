
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

filenames <- list.files("../HLI_VCF_files/")
filenames <- gsub(".vcf.gz", "", filenames)

unassigned_files_index <- which(!grepl("[[:alpha:]]", filenames))
for(id in 1:length(unassigned_files_index)){
  
  Rindex <- which(HLI_R_D_list$RID %in% as.integer(filenames[unassigned_files_index[id]]))
  if(length(Rindex) == 0){
    Dindex <- which(HLI_R_D_list$DID %in% as.integer(filenames[unassigned_files_index[id]]))
    R_D <- "D"
    Index <- Dindex
    
    outcomeDindex <- which(HLI_outcome_list$did %in% as.integer(filenames[unassigned_files_index[id]]))
    if(length(outcomeDindex) == 1){
      if(HLI_outcome_list$newagvhdgrp[outcomeDindex] == 1) GVHD <- "n" else GVHD <- "a"
    } else GVHD <- "X"
    
    
  }else {
    R_D <- "R"
    Index <- Rindex
    
    outcomeDindex <- which(HLI_outcome_list$rid %in% as.integer(filenames[unassigned_files_index[id]]))
    if(length(outcomeDindex) == 1){
      if(HLI_outcome_list$newagvhdgrp[outcomeDindex] == 1) GVHD <- "n" else GVHD <- "a"
    } else GVHD <- "X"
  }
  
  if(Index !=0){
    file.rename(from = paste0("../HLI_VCF_files/",filenames[unassigned_files_index[id]], ".vcf.gz"), 
                to = paste0("../HLI_VCF_files/", GVHD, "_", Index,"_", R_D, "_", filenames[unassigned_files_index[id]], ".vcf.gz"))
  }
  Index <- 0
}

#### Undo
filenames <- list.files("../HLI_VCF_files/")
filenames <- gsub(".vcf.gz", "", filenames)

assigned_files_index <- which(grepl("[[:alpha:]]", filenames))

for(id in 1:length(assigned_files_index)){
  new_filename <- strsplit(filenames[assigned_files_index[id]], "_")[[1]][4]
  file.rename(from = paste0("../HLI_VCF_files/",filenames[assigned_files_index[id]], ".vcf.gz"), 
              to = paste0("../HLI_VCF_files/", new_filename, ".vcf.gz"))
  
}
