KnownMiHA_table <- read.csv("../ClinVar/Data/KnownMiHA_Table_Ref.csv", stringsAsFactors = F)
Known_MiHA_SNPs <- read.delim("../WW_MiHA/Known_MiHA_coordinates.txt")

############################ GWAS HLI overlapping SNPs #############
# load("../Data/GWAS_HLI_overlapping_SNPs.RData") # overlapping_SNPs
# overlapping_SNPs$SNP_currentID <- as.character(overlapping_SNPs$SNP_currentID)
# num_overlapping_SNPs <- dim(overlapping_SNPs)[1]
# overlapping_SNPs$HLA <- character(num_overlapping_SNPs)
# for(id in 1:num_overlapping_SNPs){
#   
#   ind <- which(KnownMiHA_table$rs.number %in% overlapping_SNPs$SNP_currentID[id])
#   
#   if(length(ind) == 1){
#     
#     overlapping_SNPs$HLA[id] <- KnownMiHA_table$HLA[ind]
#     
#   }else{
#     
#     tempHLA <- unique(KnownMiHA_table$HLA[ind])
#     if(length(tempHLA) == 1){
#       
#       overlapping_SNPs$HLA[id] <- tempHLA
#       
#     }else{
#       
#       cat(id)
#       
#     }
#     
#   }
#   
# }
# 
# tempHLA <- unique(KnownMiHA_table$HLA[KnownMiHA_table$rs.number %in% overlapping_SNPs$SNP_currentID[17]])
# overlapping_SNPs$HLA[17] <- tempHLA[1]
# overlapping_SNPs <- rbind(overlapping_SNPs, overlapping_SNPs[17,])
# overlapping_SNPs$HLA[dim(overlapping_SNPs)[1]] <- tempHLA[2]
# overlapping_SNPs$HLA_SNP <- sapply(1:dim(overlapping_SNPs)[1], function(x) paste0(overlapping_SNPs$HLA[x], "-", overlapping_SNPs$SNP_currentID[x]))
# 
# save(overlapping_SNPs, file = "../Data/GWAS_HLI_overlapping_SNPs_wHLA.RData")

################# All HLI cohort All MiHA
Individual_MiHA_table_fp <- "../Output/wwangMiHAIP_DtoR/RestrictedMiHA/"
Individual_MiHA_table_files <- list.files(Individual_MiHA_table_fp, pattern = "\\.RData")
num_files <- length(Individual_MiHA_table_files)

#### Reformat MiHA database
# KnownMiHA_table <- read.csv("../ClinVar/Data/KnownMiHA_Table_Ref.csv", stringsAsFactors = F) 
# KnownMiHA_table$HLA_SNP <- sapply(1:dim(KnownMiHA_table)[1], function(x) paste0(KnownMiHA_table$HLA[x], "-", KnownMiHA_table$rs.number[x]))
#   
# Restircted_MiHA_db <- KnownMiHA_table[, c(1,3,5,9,10)]
# rm_id <- c(which(Restircted_MiHA_db$HLA == ""), which(Restircted_MiHA_db$rs.number == ""))
# Restircted_MiHA_db <- Restircted_MiHA_db[-rm_id, ]
# 
# Restircted_MiHA_db$CHROM <- sapply(1:dim(Restircted_MiHA_db)[1], function(x) 
#   if(length(which(Known_MiHA_SNPs$SNPs %in% Restircted_MiHA_db$rs.number[x]))>0) {as.character(Known_MiHA_SNPs$Chr[which(Known_MiHA_SNPs$SNPs %in% Restircted_MiHA_db$rs.number[x])]) } else 0)
# rm_id <- which(Restircted_MiHA_db$CHROM == "0")
# Restircted_MiHA_db <- Restircted_MiHA_db[-rm_id, ]
# Restircted_MiHA_db$POS <- sapply(1:dim(Restircted_MiHA_db)[1], function(x)
#   Known_MiHA_SNPs$Pos[which(Known_MiHA_SNPs$SNPs %in% Restircted_MiHA_db$rs.number[x])])
# 
# Restircted_MiHA_db$CHROM_POS <- sapply(1:dim(Restircted_MiHA_db)[1], function(x)
#   paste0(Restircted_MiHA_db$CHROM[x], "-", Restircted_MiHA_db$POS[x]))
# 
# Restircted_MiHA_db <- Restircted_MiHA_db[-which(duplicated(Restircted_MiHA_db$HLA_SNP)), ]
# 
# save(Restircted_MiHA_db, file = "../Data/Restricted_MiHA_database_RMduplicated.RData")

load("../Data/Restricted_MiHA_database_RMduplicated.RData") ## Restircted_MiHA_db
load("../Data/ID_table_wCaseID.RData")          ## ID_table
# load("../Data/HLI_reformatted_HLA_table.RData") ## reformated_HLA_typing_list
# reformated_HLA_typing_list$HLA_A1 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("A*", reformated_HLA_typing_list$HLA_A1[x]))
# reformated_HLA_typing_list$HLA_A2 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("A*", reformated_HLA_typing_list$HLA_A2[x]))
# 
# reformated_HLA_typing_list$HLA_B1 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("B*", reformated_HLA_typing_list$HLA_B1[x]))
# reformated_HLA_typing_list$HLA_B2 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("B*", reformated_HLA_typing_list$HLA_B2[x]))
# 
# reformated_HLA_typing_list$HLA_C1 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("C*", reformated_HLA_typing_list$HLA_C1[x]))
# reformated_HLA_typing_list$HLA_C2 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("C*", reformated_HLA_typing_list$HLA_C2[x]))
# 
# reformated_HLA_typing_list$HLA_DRB11 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("DRB1*", reformated_HLA_typing_list$HLA_DRB11[x]))
# reformated_HLA_typing_list$HLA_DRB12 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("DRB1*", reformated_HLA_typing_list$HLA_DRB12[x]))
# 
# reformated_HLA_typing_list$HLA_DQB11 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("DQB1*", reformated_HLA_typing_list$HLA_DQB11[x]))
# reformated_HLA_typing_list$HLA_DQB12 <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x)
#   paste0("DQB1*", reformated_HLA_typing_list$HLA_DQB12[x]))
# save(reformated_HLA_typing_list, file = "../Data/HLI_reformatted_HLA_table_corrected.RData")
load("../Data/HLI_reformatted_HLA_table_corrected.RData") ## reformated_HLA_typing_list

num_Restricted_MiHA <- dim(Restircted_MiHA_db)[1]

aGVHD_Restricted_MiHA_table <- data.frame(MiHA_SNP = Restircted_MiHA_db$HLA_SNP, stringsAsFactors = F)
aGVHD_Restricted_MiHA_table <- cbind(aGVHD_Restricted_MiHA_table, matrix(0, nrow = num_Restricted_MiHA, ncol = 32))
aGVHD_Restricted_MiHA_table_count <- data.frame(HLA_SNP = Restircted_MiHA_db$HLA_SNP, 
                                                SNP =  Restircted_MiHA_db$rs.number, 
                                                Count = numeric(length(Restircted_MiHA_db$HLA_SNP)),
                                                stringsAsFactors = F)
nGVHD_Restricted_MiHA_table <- aGVHD_Restricted_MiHA_table
nGVHD_Restricted_MiHA_table_count <- aGVHD_Restricted_MiHA_table_count

for(id in 1:num_files){
  
  eval(parse(text=paste0("load (\"", Individual_MiHA_table_fp, Individual_MiHA_table_files[id], "\")"))) # MiHA_table
  
  groupID <- as.numeric(unlist(strsplit(x = Individual_MiHA_table_files[id], split = "_"))[1])
  if(groupID == 0) stop(id)
  HLA_table_id <- which(reformated_HLA_typing_list$GroupID == groupID)
  groupType <- reformated_HLA_typing_list$GroupType[HLA_table_id]
  num_pos <- dim(MiHA_table)[1]
  
  for(jd in 1:num_pos){
    
    MiHA_db_ind <- which(Restircted_MiHA_db$CHROM_POS %in% MiHA_table$CHROM_POS[jd])
    # if(length(MiHA_db_ind) > 1) stop(jd)
    Restricted_HLA <- unique(Restircted_MiHA_db$HLA[MiHA_db_ind])
    # length(which(reformated_HLA_typing_list[HLA_table_id, 8:17] %in% Restricted_HLA))
    matched_HLA_id <- which(Restricted_HLA %in% as.character(reformated_HLA_typing_list[HLA_table_id, 8:17]))
    
    if(length(matched_HLA_id) == 1){ # restricted MiHA Presence
    # if(is.element(Restricted_HLA, as.character(reformated_HLA_typing_list[HLA_table_id, 8:17]))){ # restricted MiHA Presence
      
      for(kd in 1:length(matched_HLA_id)){
        
        if(groupType =="a"){
          MiHA_db_ind_id <- which(Restircted_MiHA_db$HLA[MiHA_db_ind] %in% Restricted_HLA[matched_HLA_id[kd]])
          
          aGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] <- aGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] + 1
          firstZeroIndex1 <- which(aGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], -1] == 0)[1] + 1
          aGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], firstZeroIndex1] <- groupID
          
        }else{
          MiHA_db_ind_id <- which(Restircted_MiHA_db$HLA[MiHA_db_ind] %in% Restricted_HLA[matched_HLA_id[kd]])
          
          nGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] <- nGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] + 1
          firstZeroIndex2 <- which(nGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], -1] == 0)[1] + 1
          nGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], firstZeroIndex2] <- groupID
          
        }
        
      }
      
    }
    
  }
  
  
}

total_counts <- cbind(aGVHD_Restricted_MiHA_table_count, nGVHD_Restricted_MiHA_table_count$Count)
colnames(total_counts)[3:4] <- c("aGVHD", "nGVHD")
write.csv(total_counts, file = "../FirstPaper/Table/HLI_Restricted_MiHA_All_counts_corrected.csv", row.names = F)

total_counts_rm_0s <-total_counts[-which(rowSums(total_counts[,3:4]) == 0), ]
write.csv(total_counts_rm_0s, file = "../FirstPaper/Table/HLI_Restricted_MiHA_All_counts_rm0_CORRECTED.csv", row.names = F)

#### add MiHA-HLA-SNP-gene
num_MiHas <- dim(total_counts_rm_0s)[1]
comprehensive_MiHA_count_table <- data.frame(MiHA = character(num_MiHas),
                                             Restricted_HLA = character(num_MiHas),
                                             SNP = character(num_MiHas),
                                             Gene = character(num_MiHas),
                                             aGVHD = numeric(num_MiHas),
                                             nGVHD = numeric(num_MiHas),
                                             LLR = numeric(num_MiHas),
                                             stringsAsFactors = F)
for(id in 1:num_MiHas){
  
  index <- which(Restircted_MiHA_db$HLA_SNP %in% total_counts_rm_0s$HLA_SNP[id])
  comprehensive_MiHA_count_table$MiHA[id] <- as.character(Restircted_MiHA_db$KnownMiHA[index])
  comprehensive_MiHA_count_table$Restricted_HLA[id] <- as.character(Restircted_MiHA_db$HLA[index])
  comprehensive_MiHA_count_table$SNP[id] <- as.character(Restircted_MiHA_db$rs.number[index])
  comprehensive_MiHA_count_table$Gene[id] <- as.character(Restircted_MiHA_db$Gene[index])
  comprehensive_MiHA_count_table$aGVHD[id] <- as.character(total_counts_rm_0s$aGVHD[id])
  comprehensive_MiHA_count_table$nGVHD[id] <- total_counts_rm_0s$nGVHD[id]
  comprehensive_MiHA_count_table$LLR[id] <- log10(total_counts_rm_0s$aGVHD[id]/total_counts_rm_0s$nGVHD[id])
  
}

write.csv(comprehensive_MiHA_count_table, file = "../FirstPaper/Table/HLI_Restricted_MiHA_final_table.csv", row.names = F)

#######
load("../Data/GWAS_HLI_overlapping_SNPs_wHLA.RData") ## overlapping_SNPs
comprehensive_MiHA_count_table <- read.csv("../FirstPaper/Table/HLI_Restricted_MiHA_final_table.csv")
### overlapping SNPs table


###### beeswarm
n_cases <- num_files
miha_table <- data.frame(groupID = numeric(n_cases),
                         GVHD = character(n_cases),
                         SEX = character(n_cases),
                         # NumVar = numeric(n_cases),
                         HLA.Restricted = numeric(n_cases),
                         stringsAsFactors = F)

for(id in 1:n_cases){
  
  eval(parse(text=paste0("load (\"", Individual_MiHA_table_fp, Individual_MiHA_table_files[id], "\")"))) # MiHA_table
  
  miha_table$groupID[id] <- as.numeric(unlist(strsplit(x = Individual_MiHA_table_files[id], split = "_"))[1])

  temp_ID_table_ind <- which(reformated_HLA_typing_list$GroupID %in% miha_table$groupID[id])
  # temp_HLI_ind <- which(HLI_metadata$bmt_case_num %in% unique(ID_table$caseNumber[temp_ID_table_ind]))
  
  miha_table$SEX[id] <- paste0(reformated_HLA_typing_list$DSex[temp_ID_table_ind], ">", reformated_HLA_typing_list$RSex[temp_ID_table_ind])
  
  miha_table$GVHD[id] <- reformated_HLA_typing_list$GroupType[temp_ID_table_ind]
  
  # miha_table$NumVar[id] <- WGS_woXY[paste0(miha_table$groupID[id], ".txt"), ]

  # miha_table$HLA.Restricted[id] <- dim(MiHA_table)[1]
  HLA_table_id <- which(reformated_HLA_typing_list$GroupID == miha_table$groupID[id])
  groupType <- reformated_HLA_typing_list$GroupType[HLA_table_id]
  num_pos <- dim(MiHA_table)[1]
  
  for(jd in 1:num_pos){
    
    MiHA_db_ind <- which(Restircted_MiHA_db$CHROM_POS %in% MiHA_table$CHROM_POS[jd])
    # if(length(MiHA_db_ind) > 1) stop(jd)
    Restricted_HLA <- unique(Restircted_MiHA_db$HLA[MiHA_db_ind])
    # length(which(reformated_HLA_typing_list[HLA_table_id, 8:17] %in% Restricted_HLA))
    matched_HLA_id <- which(Restricted_HLA %in% as.character(reformated_HLA_typing_list[HLA_table_id, 8:17]))
    
    if(length(matched_HLA_id) == 1){ # restricted MiHA Presence
      # if(is.element(Restricted_HLA, as.character(reformated_HLA_typing_list[HLA_table_id, 8:17]))){ # restricted MiHA Presence
      
      for(kd in 1:length(matched_HLA_id)){
        
        miha_table$HLA.Restricted[id] <-  miha_table$HLA.Restricted[id] + 1
        
      }
      
    }
    
  } 
}
##### Restricted MiHA
library(beeswarm)
library(ggplot2)
library(plyr)
beeswarm_ResMiHA_all <- beeswarm(HLA.Restricted ~ GVHD, data = miha_table, 
                                 method = 'swarm',
                                 spacing = 0.5)[, c(1,2,6)]
colnames(beeswarm_ResMiHA_all) <- c("x", "y", "GVHD") 
# beeswarm_ResMiHA_all$x <- abs(beeswarm_ResMiHA_all$x)

ggplot(beeswarm_ResMiHA_all, aes(x, y)) +
  xlab("") +
  scale_y_continuous(expression("#Mismatched HLA Restricted Known MiHA SNPs")) + 
  # geom_point(aes(colour = GVHD), size = 3, alpha = 0.5) +
  geom_jitter(width = 0.02, aes(color = GVHD), size = 3, alpha = 0.8) +
  # geom_jitter() +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_continuous(breaks = c(1:2), 
                     labels = c("acute GVHD", "non aGVHD"), expand = c(0, 0.05)) +
  geom_boxplot(aes(x, y, group = round_any(beeswarm_ResMiHA_all$x, 1, round)), outlier.shape = NA, alpha = 0) +
  theme(legend.position = "none")

t.test(beeswarm_ResMiHA_all$y[which(beeswarm_ResMiHA_all$GVHD == "a")], beeswarm_ResMiHA_all$y[which(beeswarm_ResMiHA_all$GVHD == "n")])

wilcox.test(beeswarm_ResMiHA_all$y[which(beeswarm_ResMiHA_all$GVHD == "a")], beeswarm_ResMiHA_all$y[which(beeswarm_ResMiHA_all$GVHD == "n")])

wilcox.test(y ~ GVHD, beeswarm_ResMiHA_all)


###################### Overlapping SNPs ###################
#################
load("../Data/GWAS_HLI_overlapping_SNPs_wHLA.RData") ## overlapping_SNPs
overlapping_SNPs$CHROM_POS <- overlapping_SNPs$chr_POS
load("../Data/ID_table_wCaseID.RData")          ## ID_table

load("../Data/HLI_reformatted_HLA_table_corrected.RData") ## reformated_HLA_typing_list

num_Restricted_MiHA <- dim(overlapping_SNPs)[1]

aGVHD_Restricted_MiHA_table <- data.frame(MiHA_SNP = overlapping_SNPs$HLA_SNP, stringsAsFactors = F)
aGVHD_Restricted_MiHA_table <- cbind(aGVHD_Restricted_MiHA_table, matrix(0, nrow = num_Restricted_MiHA, ncol = 32))
aGVHD_Restricted_MiHA_table_count <- data.frame(HLA_SNP = overlapping_SNPs$HLA_SNP, 
                                                SNP =  overlapping_SNPs$SNP_currentID, 
                                                Count = numeric(length(overlapping_SNPs$HLA_SNP)),
                                                stringsAsFactors = F)
nGVHD_Restricted_MiHA_table <- aGVHD_Restricted_MiHA_table
nGVHD_Restricted_MiHA_table_count <- aGVHD_Restricted_MiHA_table_count

for(id in 1:num_files){
  
  eval(parse(text=paste0("load (\"", Individual_MiHA_table_fp, Individual_MiHA_table_files[id], "\")"))) # MiHA_table
  
  groupID <- as.numeric(unlist(strsplit(x = Individual_MiHA_table_files[id], split = "_"))[1])
  if(groupID == 0) stop(id)
  HLA_table_id <- which(reformated_HLA_typing_list$GroupID == groupID)
  groupType <- reformated_HLA_typing_list$GroupType[HLA_table_id]
  num_pos <- dim(MiHA_table)[1]
  
  for(jd in 1:num_pos){
    
    MiHA_db_ind <- which(overlapping_SNPs$CHROM_POS %in% MiHA_table$CHROM_POS[jd])
    # if(length(MiHA_db_ind) > 1) stop(jd)
    Restricted_HLA <- unique(overlapping_SNPs$HLA[MiHA_db_ind])
    # length(which(reformated_HLA_typing_list[HLA_table_id, 8:17] %in% Restricted_HLA))
    matched_HLA_id <- which(Restricted_HLA %in% as.character(reformated_HLA_typing_list[HLA_table_id, 8:17]))
    
    if(length(matched_HLA_id) == 1){ # restricted MiHA Presence
      # if(is.element(Restricted_HLA, as.character(reformated_HLA_typing_list[HLA_table_id, 8:17]))){ # restricted MiHA Presence
      
      for(kd in 1:length(matched_HLA_id)){
        
        if(groupType =="a"){
          MiHA_db_ind_id <- which(overlapping_SNPs$HLA[MiHA_db_ind] %in% Restricted_HLA[matched_HLA_id[kd]])
          
          aGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] <- aGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] + 1
          firstZeroIndex1 <- which(aGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], -1] == 0)[1] + 1
          aGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], firstZeroIndex1] <- groupID
          
        }else{
          MiHA_db_ind_id <- which(overlapping_SNPs$HLA[MiHA_db_ind[MiHA_db_ind_id]] %in% Restricted_HLA[matched_HLA_id[kd]])
          
          nGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] <- nGVHD_Restricted_MiHA_table_count$Count[MiHA_db_ind[MiHA_db_ind_id]] + 1
          firstZeroIndex2 <- which(nGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], -1] == 0)[1] + 1
          nGVHD_Restricted_MiHA_table[MiHA_db_ind[MiHA_db_ind_id], firstZeroIndex2] <- groupID
          
        }
        
      }
      
    }
    
  }
  
  
}

total_counts <- cbind(aGVHD_Restricted_MiHA_table_count, nGVHD_Restricted_MiHA_table_count$Count)
colnames(total_counts)[3:4] <- c("aGVHD", "nGVHD")
write.csv(total_counts, file = "../FirstPaper/Table/HLI_overlapping_Restricted_MiHA_counts.csv", row.names = F)


#####################
# known_MiHA_coordinates <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T, stringsAsFactors = F)
load("../Data/Restricted_MiHA_database_RMduplicated.RData") # Restircted_MiHA_db
num_rows <- dim(total_counts)[1] 
total_counts$MiHA_HLA_SNP_Gene <- sapply(1:num_rows, function(x) {ind <- which(Restircted_MiHA_db$HLA_SNP == total_counts$HLA_SNP[x])
paste0(Restircted_MiHA_db$KnownMiHA[ind], "<>", Restircted_MiHA_db$HLA[x], "<>", Restircted_MiHA_db$Gene[ind])}) 

num_unique_MiHAs <- dim(total_counts)[1]

MiHA_LLR <- data.frame(MiHA_HLA_SNP_Gene = character(num_unique_MiHAs),
                       aGVHD_count = numeric(num_unique_MiHAs),
                       nGVHD_count = numeric(num_unique_MiHAs),
                       LLR = numeric(num_unique_MiHAs),
                       stringsAsFactors = F)

for(id in 1:num_unique_MiHAs){
  
  MiHA_LLR$MiHA_HLA_SNP_Gene[id] <- total_counts$MiHA_HLA_SNP_Gene[id]
  MiHA_LLR$aGVHD_count[id] <- total_counts$aGVHD[id]
  MiHA_LLR$nGVHD_count[id] <- total_counts$nGVHD[id]
  MiHA_LLR$LLR[id] <- log10(MiHA_LLR$aGVHD_count[id] / MiHA_LLR$nGVHD_count[id])
  
}

write.csv(MiHA_LLR, file = "../FirstPaper/Table/HLI_MiHA_LLR_table.csv", row.names = F)

############################ Reformat HLI HLA table #############################################
# HLA_typings <- read.csv(file = "../HLI_hla_mg_v3.csv")
# 
# available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]
# 
# available_IDs_HLA_typing <- HLA_typings[(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID),]
# num_groups <- dim(available_IDs_HLA_typing)[1]
# group_type <- sapply(1:num_groups, function(x) available_IDs$Group[available_IDs$R_D_ID %in% available_IDs_HLA_typing$nmdp_rid[x]])
# group_type[group_type == "a"] <- "aGVHD"
# group_type[group_type == "n"] <- "non-GVHD"
# 
# reformated_HLA_typing_list <- available_IDs_HLA_typing[, 1:3]
# reformated_HLA_typing_list$RRace <- available_IDs_HLA_typing$rid.broad.race
# reformated_HLA_typing_list$DRace <- available_IDs_HLA_typing$did.broad.race
# reformated_HLA_typing_list$RSex <- available_IDs_HLA_typing$rid_sex
# reformated_HLA_typing_list$DSex <- available_IDs_HLA_typing$donor_sex
# 
# reformated_HLA_typing_list$HLA_A1 <- sapply(as.character(available_IDs_HLA_typing$r_a_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_A2 <- sapply(as.character(available_IDs_HLA_typing$r_a_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_B1 <- sapply(as.character(available_IDs_HLA_typing$r_b_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_B2 <- sapply(as.character(available_IDs_HLA_typing$r_b_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_C1 <- sapply(as.character(available_IDs_HLA_typing$r_c_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_C2 <- sapply(as.character(available_IDs_HLA_typing$r_c_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_DRB11 <- sapply(as.character(available_IDs_HLA_typing$r_drb1_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_DRB12 <- sapply(as.character(available_IDs_HLA_typing$r_drb1_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_DQB11 <- sapply(as.character(available_IDs_HLA_typing$r_dqb1_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
# reformated_HLA_typing_list$HLA_DQB12 <- sapply(as.character(available_IDs_HLA_typing$r_dqb1_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
# 
# get_two_field <- function(HLA_typing){
#   
#   reformat_typing <- unlist(strsplit(HLA_typing, ":"))
#   if(length(reformat_typing) > 2){
#     
#     two_fields <- paste0(reformat_typing[c(1,2)], collapse = ":")
#     
#   }else two_fields <- HLA_typing
#   
#   return(two_fields)
# }
# 
# reformated_HLA_typing_list$HLA_A1 <- sapply(reformated_HLA_typing_list$HLA_A1, get_two_field)
# reformated_HLA_typing_list$HLA_A2 <- sapply(reformated_HLA_typing_list$HLA_A2, get_two_field)
# reformated_HLA_typing_list$HLA_B1 <- sapply(reformated_HLA_typing_list$HLA_B1, get_two_field)
# reformated_HLA_typing_list$HLA_B2 <- sapply(reformated_HLA_typing_list$HLA_B2, get_two_field)
# reformated_HLA_typing_list$HLA_C1 <- sapply(reformated_HLA_typing_list$HLA_C1, get_two_field)
# reformated_HLA_typing_list$HLA_C2 <- sapply(reformated_HLA_typing_list$HLA_C2, get_two_field)
# reformated_HLA_typing_list$HLA_DRB11 <- sapply(reformated_HLA_typing_list$HLA_DRB11, get_two_field)
# reformated_HLA_typing_list$HLA_DRB12 <- sapply(reformated_HLA_typing_list$HLA_DRB12, get_two_field)
# reformated_HLA_typing_list$HLA_DQB11 <- sapply(reformated_HLA_typing_list$HLA_DQB11, get_two_field)
# reformated_HLA_typing_list$HLA_DQB12 <- sapply(reformated_HLA_typing_list$HLA_DQB12, get_two_field)
# 
# reformated_HLA_typing_list$GroupID <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x) 
#   unique(ID_table$GroupID[which(ID_table$caseID %in% reformated_HLA_typing_list$bmt_case_num[x])]))
# reformated_HLA_typing_list$GroupType <- sapply(1:dim(reformated_HLA_typing_list)[1], function(x) 
#   unique(ID_table$Group[which(ID_table$caseID %in% reformated_HLA_typing_list$bmt_case_num[x])]))
# 
# save(reformated_HLA_typing_list, file = "../Data/HLI_reformatted_HLA_table.RData")

####### Use All Table from WW #####################################################################################################
known_miha_freq <- read.delim("../FirstPaper/Data/Updated_Restricted_KnownMiHA.txt", header = F, stringsAsFactors = F)
colnames(known_miha_freq) <- c("SNP", "HLA", "GroupID", "CHROM", "POS", "V6", "REF", "ALT")
# table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq2 <- unique(known_miha_freq[, c(1:5, 8)])

known_miha_freq2$HLA_SNP <- sapply(1:dim(known_miha_freq2)[1],
                                  function(x) paste0(known_miha_freq2$HLA[x], "-", known_miha_freq2$SNP[x]))
unique_groups_groupID <- unique(known_miha_freq2[, 1:3])

table(known_miha_freq2[, c(1,3)])
groupID <- unique(known_miha_freq2$PID[known_miha_freq2$HLA == "B*07:02"])

load("../Data/ID_table_wCaseID.RData")
HLI_metadata <- read.csv("../HLI_hla_mg_v3.csv", stringsAsFactors = F)
caseID <- sapply(1:length(gorupID), function(x) unique(ID_table$caseID[which(ID_table$GroupID == groupID[x])]))
aa <- HLI_metadata[which(HLI_metadata$bmt_case_num %in% caseID), c(1, 85:92, 101, 102)]

#####
known_miha_freq <- read.delim("../WW_MiHA/All_restricted_MiHAs.txt", header = T, stringsAsFactors = F)
# colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")

# table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq$HLA_SNP <- sapply(1:dim(known_miha_freq)[1], 
                                  function(x) paste0(known_miha_freq$HLA[x], "-", known_miha_freq$SNP[x]))

#################
single_MiHA_stats <- as.data.frame(table(known_miha_freq[c("GVHD", "HLA_SNP")]))

GroupID_MiHA_stats <- as.data.frame(table(known_miha_freq[c("GVHD", "HLA_SNP")]))
GroupID_MiHA_stats <- GroupID_MiHA_stats[GroupID_MiHA_stats$Freq>0, ]

length(unique(single_MiHA_stats$HLA_SNP))
length(unique(GroupID_MiHA_stats$HLA_SNP))
length(unique(known_miha_freq$SNP)) # 37 MiHA SNPs

unique_IDs <- unique(known_miha_freq$PID) # 131 groups

SNP_mat <- matrix(data = 0, nrow = length(unique_IDs), ncol = length(unique(GroupID_MiHA_stats$HLA_SNP)) + 2)
colnames(SNP_mat) <- c("GroupID", "GroupType", as.character(unique(GroupID_MiHA_stats$HLA_SNP)))
SNP_mat <- as.data.frame(SNP_mat)

for(id in 1:length(unique_IDs)){
  
  selected_group <- known_miha_freq[known_miha_freq$PID == unique_IDs[id],]
  SNP_mat$GroupID[id] <- unique_IDs[id]
  SNP_mat$GroupType[id] <- ID_table$Group[which(ID_table$GroupID == unique_IDs[id])[1]]
  
  if(length(unique(selected_group$HLA_SNP)) == length(selected_group$HLA_SNP)){ # no duplicated SNPs for one group
    
    SNP_mat[id, selected_group$HLA_SNP] <- 1
    
  }else{ # one homozygous SNP (duplicated)
    
    SNP_counts <- as.data.frame(table(selected_group$HLA_SNP), stringsAsFactors = F)
    SNP_mat[id, SNP_counts[, 1]] <- SNP_counts[, 2]
    
  }
  
}
# write.csv(SNP_mat, file = "../WW_MiHA/SNP_mat.csv", row.names = FALSE)

sort(colSums(SNP_mat[, -c(1,2)]), decreasing = T)
sort(colSums(SNP_mat[which(SNP_mat$GroupType == "n"), -c(1,2)]), decreasing = T)
sort(colSums(SNP_mat[which(SNP_mat$GroupType == "a"), -c(1,2)]), decreasing = T)

# presence-absence 
presence_SNP_mat <- SNP_mat[, -c(1,2)]
presence_SNP_mat[presence_SNP_mat>0] <- 1
presence_SNP_mat <- cbind(SNP_mat[, c(1,2)], presence_SNP_mat)

nonZeros_rows <- which(rowSums(presence_SNP_mat[, -c(1,2)]) !=0)

##### Log Ratio
known_MiHA_coordinates <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T, stringsAsFactors = F)
num_rows <- dim(known_miha_freq)[1] 
known_miha_freq$MiHA_HLA_SNP_Gene <- sapply(1:num_rows, function(x) {ind <- which(known_MiHA_coordinates$SNPs == known_miha_freq$SNP[x])
paste0(known_MiHA_coordinates$MiHAs[ind], "<>", known_miha_freq$HLA_SNP[x], "<>", known_MiHA_coordinates$Gene[ind])}) 

unique_MiHAs <- unique(known_miha_freq$MiHA_HLA_SNP_Gene)
num_unique_MiHAs <- length(unique_MiHAs)

MiHA_LLR <- data.frame(MiHA_HLA_SNP_Gene = character(num_unique_MiHAs),
                       aGVHD_count = numeric(num_unique_MiHAs),
                       nGVHD_count = numeric(num_unique_MiHAs),
                       LLR = numeric(num_unique_MiHAs),
                       stringsAsFactors = F)

for(id in 1:num_unique_MiHAs){
  
  index <- which(known_miha_freq$MiHA_HLA_SNP_Gene %in% unique_MiHAs[id])
  
  count_summary <- as.data.frame(table(known_miha_freq$GVHD[index]))
  
  MiHA_LLR$MiHA_HLA_SNP_Gene[id] <- unique_MiHAs[id]
  MiHA_LLR$aGVHD_count[id] <- count_summary$Freq[which(count_summary$Var1 == "a")]
  MiHA_LLR$nGVHD_count[id] <- count_summary$Freq[which(count_summary$Var1 == "n")]
  MiHA_LLR$LLR[id] <- log10(MiHA_LLR$aGVHD_count[id] / MiHA_LLR$nGVHD_count[id])
  
}

write.csv(MiHA_LLR, file = "../FirstPaper/Table/HLI_All_MiHA_LLR_WWVersion.csv", row.names = F)


###############################################
####### Use All Table from WW #####################################################################################################
known_miha_freq <- read.delim("../FirstPaper/Data/Updated_Restricted_KnownMiHA.txt", header = F, stringsAsFactors = F)
colnames(known_miha_freq) <- c("SNP", "HLA", "GroupID", "CHROM", "POS", "V6", "REF", "ALT")
# table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq2 <- unique(known_miha_freq[, c(1:5, 8)])

known_miha_freq2$HLA_SNP <- sapply(1:dim(known_miha_freq2)[1],
                                   function(x) paste0(known_miha_freq2$HLA[x], "-", known_miha_freq2$SNP[x]))
unique_groups_groupID <- unique(known_miha_freq2[, 1:3])

table(known_miha_freq2[, c(1,3)])
groupID <- unique(known_miha_freq2$PID[known_miha_freq2$HLA == "B*07:02"])

load("../Data/ID_table_wCaseID.RData")
HLI_metadata <- read.csv("../HLI_hla_mg_v3.csv", stringsAsFactors = F)
caseID <- sapply(1:length(gorupID), function(x) unique(ID_table$caseID[which(ID_table$GroupID == groupID[x])]))
aa <- HLI_metadata[which(HLI_metadata$bmt_case_num %in% caseID), c(1, 85:92, 101, 102)]

#####
known_miha_freq <- read.delim("../WW_MiHA/All_restricted_MiHAs.txt", header = T, stringsAsFactors = F)
# colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")

# table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq$HLA_SNP <- sapply(1:dim(known_miha_freq)[1], 
                                  function(x) paste0(known_miha_freq$HLA[x], "-", known_miha_freq$SNP[x]))

#################
single_MiHA_stats <- as.data.frame(table(known_miha_freq[c("GVHD", "HLA_SNP")]))

GroupID_MiHA_stats <- as.data.frame(table(known_miha_freq[c("GVHD", "HLA_SNP")]))
GroupID_MiHA_stats <- GroupID_MiHA_stats[GroupID_MiHA_stats$Freq>0, ]

length(unique(single_MiHA_stats$HLA_SNP))
length(unique(GroupID_MiHA_stats$HLA_SNP))
length(unique(known_miha_freq$SNP)) # 37 MiHA SNPs

unique_IDs <- unique(known_miha_freq$PID) # 131 groups

SNP_mat <- matrix(data = 0, nrow = length(unique_IDs), ncol = length(unique(GroupID_MiHA_stats$HLA_SNP)) + 2)
colnames(SNP_mat) <- c("GroupID", "GroupType", as.character(unique(GroupID_MiHA_stats$HLA_SNP)))
SNP_mat <- as.data.frame(SNP_mat)

for(id in 1:length(unique_IDs)){
  
  selected_group <- known_miha_freq[known_miha_freq$PID == unique_IDs[id],]
  SNP_mat$GroupID[id] <- unique_IDs[id]
  SNP_mat$GroupType[id] <- ID_table$Group[which(ID_table$GroupID == unique_IDs[id])[1]]
  
  if(length(unique(selected_group$HLA_SNP)) == length(selected_group$HLA_SNP)){ # no duplicated SNPs for one group
    
    SNP_mat[id, selected_group$HLA_SNP] <- 1
    
  }else{ # one homozygous SNP (duplicated)
    
    SNP_counts <- as.data.frame(table(selected_group$HLA_SNP), stringsAsFactors = F)
    SNP_mat[id, SNP_counts[, 1]] <- SNP_counts[, 2]
    
  }
  
}
# write.csv(SNP_mat, file = "../WW_MiHA/SNP_mat.csv", row.names = FALSE)

sort(colSums(SNP_mat[, -c(1,2)]), decreasing = T)
sort(colSums(SNP_mat[which(SNP_mat$GroupType == "n"), -c(1,2)]), decreasing = T)
sort(colSums(SNP_mat[which(SNP_mat$GroupType == "a"), -c(1,2)]), decreasing = T)

# presence-absence 
presence_SNP_mat <- SNP_mat[, -c(1,2)]
presence_SNP_mat[presence_SNP_mat>0] <- 1
presence_SNP_mat <- cbind(SNP_mat[, c(1,2)], presence_SNP_mat)

nonZeros_rows <- which(rowSums(presence_SNP_mat[, -c(1,2)]) !=0)

##### Log Ratio
known_MiHA_coordinates <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T, stringsAsFactors = F)
num_rows <- dim(known_miha_freq)[1] 
known_miha_freq$MiHA_HLA_SNP_Gene <- sapply(1:num_rows, function(x) {ind <- which(known_MiHA_coordinates$SNPs == known_miha_freq$SNP[x])
paste0(known_MiHA_coordinates$MiHAs[ind], "<>", known_miha_freq$HLA_SNP[x], "<>", known_MiHA_coordinates$Gene[ind])}) 

unique_MiHAs <- unique(known_miha_freq$MiHA_HLA_SNP_Gene)
num_unique_MiHAs <- length(unique_MiHAs)

MiHA_LLR <- data.frame(MiHA_HLA_SNP_Gene = character(num_unique_MiHAs),
                       aGVHD_count = numeric(num_unique_MiHAs),
                       nGVHD_count = numeric(num_unique_MiHAs),
                       LLR = numeric(num_unique_MiHAs),
                       stringsAsFactors = F)

for(id in 1:num_unique_MiHAs){
  
  index <- which(known_miha_freq$MiHA_HLA_SNP_Gene %in% unique_MiHAs[id])
  
  count_summary <- as.data.frame(table(known_miha_freq$GVHD[index]))
  
  MiHA_LLR$MiHA_HLA_SNP_Gene[id] <- unique_MiHAs[id]
  MiHA_LLR$aGVHD_count[id] <- count_summary$Freq[which(count_summary$Var1 == "a")]
  MiHA_LLR$nGVHD_count[id] <- count_summary$Freq[which(count_summary$Var1 == "n")]
  MiHA_LLR$LLR[id] <- log10(MiHA_LLR$aGVHD_count[id] / MiHA_LLR$nGVHD_count[id])
  
}


############################ Wei's Updated Table 06/29/2017
# restricted_miha_table <- read.delim(file = "../FirstPaper/Data/Updated_Restricted_KnownMiHA.txt", header = F, stringsAsFactors = F)
restricted_miha_table <- read.delim(file = "../FirstPaper/Data/Restricted_MismatchMiHAs_with_correctedTyping(Wei_updated).txt", header = F, stringsAsFactors = F)
restricted_miha_table <- restricted_miha_table[, c(1:8, 10)]
colnames(restricted_miha_table) <- c("SNP", "Allele", "GroupID", "CHROM", "POS", "V6", "Donor", "Recipient", "Group")
restricted_miha_table$HLA <- sapply(1:dim(restricted_miha_table)[1], function(x) {
  temp_allele <- unlist(strsplit(restricted_miha_table$Allele[x], ""))
  temp_gene <-  paste0(temp_allele[1:(length(temp_allele)-4)], collapse = "")
  paste0(temp_gene, "*", temp_allele[(length(temp_allele)-3)], temp_allele[(length(temp_allele)-2)], ":", temp_allele[(length(temp_allele)-1)], temp_allele[(length(temp_allele)-0)])
  
}
)
load("../Data/HLI_reformatted_HLA_table_corrected.RData") ## reformated_HLA_typing_list

restricted_miha_table$Restricted_HLA <- sapply(1:dim(restricted_miha_table)[1], function(x) is.element(restricted_miha_table$HLA[x], reformated_HLA_typing_list[which(reformated_HLA_typing_list$GroupID == restricted_miha_table$GroupID[x]),8:17]))
which(!restricted_miha_table$Restricted_HLA)
restricted_miha_table$HLA_SNP <- sapply(1:dim(restricted_miha_table)[1], function(x) paste0(restricted_miha_table$HLA[x], "-", restricted_miha_table$SNP[x]))
#################
single_MiHA_stats <- as.data.frame(table(restricted_miha_table[c("Group", "HLA_SNP")]))

length(unique(single_MiHA_stats$HLA_SNP))
# length(unique(GroupID_MiHA_stats$HLA_SNP))
# length(unique(known_miha_freq$SNP)) # 37 MiHA SNPs

unique_IDs <- unique(restricted_miha_table$GroupID) # 170 groups

SNP_mat <- matrix(data = 0, nrow = length(unique_IDs), ncol = length(unique(GroupID_MiHA_stats$HLA_SNP)) + 2)
colnames(SNP_mat) <- c("GroupID", "GroupType", as.character(unique(GroupID_MiHA_stats$HLA_SNP)))
SNP_mat <- as.data.frame(SNP_mat)

for(id in 1:length(unique_IDs)){
  
  selected_group <- known_miha_freq[known_miha_freq$PID == unique_IDs[id],]
  SNP_mat$GroupID[id] <- unique_IDs[id]
  SNP_mat$GroupType[id] <- ID_table$Group[which(ID_table$GroupID == unique_IDs[id])[1]]
  
  if(length(unique(selected_group$HLA_SNP)) == length(selected_group$HLA_SNP)){ # no duplicated SNPs for one group
    
    SNP_mat[id, selected_group$HLA_SNP] <- 1
    
  }else{ # one homozygous SNP (duplicated)
    
    SNP_counts <- as.data.frame(table(selected_group$HLA_SNP), stringsAsFactors = F)
    SNP_mat[id, SNP_counts[, 1]] <- SNP_counts[, 2]
    
  }
  
}
# write.csv(SNP_mat, file = "../WW_MiHA/SNP_mat.csv", row.names = FALSE)

sort(colSums(SNP_mat[, -c(1,2)]), decreasing = T)
sort(colSums(SNP_mat[which(SNP_mat$GroupType == "n"), -c(1,2)]), decreasing = T)
sort(colSums(SNP_mat[which(SNP_mat$GroupType == "a"), -c(1,2)]), decreasing = T)

# presence-absence 
presence_SNP_mat <- SNP_mat[, -c(1,2)]
presence_SNP_mat[presence_SNP_mat>0] <- 1
presence_SNP_mat <- cbind(SNP_mat[, c(1,2)], presence_SNP_mat)

nonZeros_rows <- which(rowSums(presence_SNP_mat[, -c(1,2)]) !=0)

##### Log Ratio
known_MiHA_coordinates <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T, stringsAsFactors = F)
num_rows <- dim(known_miha_freq)[1] 
known_miha_freq$MiHA_HLA_SNP_Gene <- sapply(1:num_rows, function(x) {ind <- which(known_MiHA_coordinates$SNPs == known_miha_freq$SNP[x])
paste0(known_MiHA_coordinates$MiHAs[ind], "<>", known_miha_freq$HLA_SNP[x], "<>", known_MiHA_coordinates$Gene[ind])}) 

unique_MiHAs <- unique(known_miha_freq$MiHA_HLA_SNP_Gene)
num_unique_MiHAs <- length(unique_MiHAs)

MiHA_LLR <- data.frame(MiHA_HLA_SNP_Gene = character(num_unique_MiHAs),
                       aGVHD_count = numeric(num_unique_MiHAs),
                       nGVHD_count = numeric(num_unique_MiHAs),
                       LLR = numeric(num_unique_MiHAs),
                       stringsAsFactors = F)

for(id in 1:num_unique_MiHAs){
  
  index <- which(known_miha_freq$MiHA_HLA_SNP_Gene %in% unique_MiHAs[id])
  
  count_summary <- as.data.frame(table(known_miha_freq$GVHD[index]))
  
  MiHA_LLR$MiHA_HLA_SNP_Gene[id] <- unique_MiHAs[id]
  MiHA_LLR$aGVHD_count[id] <- count_summary$Freq[which(count_summary$Var1 == "a")]
  MiHA_LLR$nGVHD_count[id] <- count_summary$Freq[which(count_summary$Var1 == "n")]
  MiHA_LLR$LLR[id] <- log10(MiHA_LLR$aGVHD_count[id] / MiHA_LLR$nGVHD_count[id])
  
}
