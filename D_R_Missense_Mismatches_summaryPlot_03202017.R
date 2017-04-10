source("util.R")
# library(BSgenome.Hsapiens.UCSC.hg38)
library(ggplot2)

##### Build a comprehensive group table
# load("../Data/ID_table.RData")
# available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]
# ###### disease Types -- 203 groups
# disease_types <- read.csv("../ClinVar/HLI_DiseaseTypes_2.csv")
# HLI_ID <- read.csv("../HLI_ID_pairs.csv")
# did <- sapply(1:dim(HLI_ID)[1], function(x) as.numeric(gsub("-", "", as.character(HLI_ID$DID[x]))))
# rid <- sapply(1:dim(HLI_ID)[1], function(x) as.numeric(gsub("-", "", as.character(HLI_ID$RID[x]))))
# HLI_ID$DID_num <- did
# HLI_ID$RID_num <- rid 
# 
# num_groups <- 205
# Available_paired_table <- data.frame(GroupID = numeric(num_groups),
#                                      GroupType = character(num_groups),
#                                      BMT = numeric(num_groups),
#                                      DID = numeric(num_groups),
#                                      RID = numeric(num_groups),
#                                      GVHDGRP = numeric(num_groups),
#                                      Disease = character(num_groups),
#                                      stringsAsFactors = F)
# for(id in 1:num_groups){
#   
#   Available_paired_table$GroupID[id] <- available_IDs$GroupID[2*(id-1)+1]
#   Available_paired_table$GroupType[id] <- available_IDs$Group[2*(id-1)+1]
#   Available_paired_table$DID[id] <- available_IDs$R_D_ID[2*(id-1)+1]
#   Available_paired_table$RID[id] <- available_IDs$R_D_ID[2*id]
#   
#   Available_paired_table$BMT[id] <- HLI_ID$BMT[which(HLI_ID$DID_num %in% Available_paired_table$DID[id])]
#   
#   dis_ind <- which(disease_types$did %in% Available_paired_table$DID[id])
#   if(length(dis_ind)!=0){
#     
#     Available_paired_table$Disease[id] <- as.character(disease_types$disease[dis_ind])
#     Available_paired_table$GVHDGRP[id] <- disease_types$gvhdgrp[dis_ind]
#     
#   }
#   
# }
# save(Available_paired_table, file = "../Data/HLI_available_pairs_dis_table.RData")

############ check D-R missense mismatches
load("../Data/HLI_available_pairs_dis_table.RData")
Missense_mismatche_fp <- "../ClinVar/D_R_Missense_Mismatches/"
Missense_mismatche_files <- list.files(Missense_mismatche_fp)
num_files <- length(Missense_mismatche_files)
sink(file = "../ClinVar/D_R_Missense_Mismatches/Missense_numbers.txt")
for(id in 1:num_files){
  
  missense_table <- read.table(file = paste0(Missense_mismatche_fp,Missense_mismatche_files[id]))
  missense_table$V1 <- factor(missense_table$V1, levels = sapply(c(1:22, "X", "Y"),function(x) paste0("chr", x)))
  print(Missense_mismatche_files[id])
  print(table(missense_table$V1))
  
}
sink()

load("../Output/Missense_variant_stats/RefSeq_canon_0228/Matched_donor_missesense_stats_RefSeq_Canon_0228.RData")
