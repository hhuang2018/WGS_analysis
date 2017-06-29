KnownMiHA_table <- read.csv("../ClinVar/Data/KnownMiHA_Table_Ref.csv", stringsAsFactors = F)
Known_MiHA_SNPs <- read.delim("../WW_MiHA/Known_MiHA_coordinates.txt")
load("../Data/GWAS_HLI_overlapping_SNPs.RData") # overlapping_SNPs
overlapping_SNPs$SNP_currentID <- as.character(overlapping_SNPs$SNP_currentID) 
num_overlapping_SNPs <- dim(overlapping_SNPs)[1]
overlapping_SNPs$HLA <- character(num_overlapping_SNPs)
for(id in 1:num_overlapping_SNPs){
  
  ind <- which(KnownMiHA_table$rs.number %in% overlapping_SNPs$SNP_currentID[id])
  
  if(length(ind) == 1){
    
    overlapping_SNPs$HLA[id] <- KnownMiHA_table$HLA[ind]
    
  }else{
    
    tempHLA <- unique(KnownMiHA_table$HLA[ind])
    if(length(tempHLA) == 1){
      
      overlapping_SNPs$HLA[id] <- tempHLA
      
    }else{
      
      cat(id)
      
    }
    
  }
  
}


# known_miha_freq <- read.delim("../WW_MiHA/All_restricted_MiHAs.txt", header = T, stringsAsFactors = F)
# # colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")
# # table(known_miha_freq[c("HLA_type","SNP")])
# 
# known_miha_freq2 <- unique(known_miha_freq[, c(1:5, 8)])
# 
# known_miha_freq2$HLA_SNP <- sapply(1:dim(known_miha_freq2)[1], 
#                                   function(x) paste0(known_miha_freq2$HLA[x], "-", known_miha_freq2$SNP[x]))
# unique_groups_groupID <- unique(known_miha_freq2[, 1:3])
# 
# table(known_miha_freq2[, c(1,3)])
# groupID <- unique(known_miha_freq2$PID[known_miha_freq2$HLA == "B*07:02"])
# 
# load("../Data/ID_table_wCaseID.RData")
# HLI_metadata <- read.csv("../HLI_hla_mg_v3.csv", stringsAsFactors = F)
# caseID <- sapply(1:length(gorupID), function(x) unique(ID_table$caseID[which(ID_table$GroupID == groupID[x])]))
# aa <- HLI_metadata[which(HLI_metadata$bmt_case_num %in% caseID), c(1, 85:92, 101, 102)]

file_list_fp <- "../Output/MiHA_missense_Summary/"
files <- list.files(file_list_fp, pattern = "\\.RData")
# all_vcf_files <- list.files(renamed_vcf_dir, pattern = "\\.vcf.gz$")
num_files <- length(files)
for(id in 1:num_files){
  if(grepl("Summary", files[id])){
    eval(parse(text=paste0("load (\"", file_list_fp, files[id], "\")")))
  }
  
  
}