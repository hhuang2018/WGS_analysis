
# file_fp <- "/mnt/cloudbiodata_nfs_1/hli_scratch/wwang/MiHAIP/DtoR/"
# KnownMiHA_table <- read.csv("../ClinVar/Data/KnownMiHA_Table_Ref.csv", stringsAsFactors = F)
Known_MiHA_SNPs <- read.delim("../Known_MiHA_coordinates.txt")
Known_MiHA_SNPs_Pos <- unique(Known_MiHA_SNPs[, -1])
Known_MiHA_SNPs_Pos$CHROM_POS <- sapply(1:dim(Known_MiHA_SNPs_Pos)[1], function(x) paste0(Known_MiHA_SNPs_Pos$Chr[x], "-", Known_MiHA_SNPs_Pos$Pos[x]))

file_fp <- "../Output/wwangMiHAIP_DtoR/"

file_list <- list.files(file_fp, pattern = "\\.txt")

num_files <- length(file_list)

output_fp <- "/mnt/cloudbiodata_nfs_2/users/hhuang/wwMiHA_summary/"

for(id in 1:num_files){
  
  temp_tab <- read.delim(paste0(file_fp, file_list[id]), header = F, stringsAsFactors = F)
  colnames(temp_tab) <- c("CHROM", "POS", "V3", "Donor", "Recipient", "Qual", "QFilter", "V8", "V9", "V10")
  temp_tab$CHROM_POS <- sapply(1:dim(temp_tab)[1], function(x) paste0(temp_tab$CHROM[x], "-", temp_tab$POS[x]))

  ind <- which(temp_tab$CHROM_POS %in% Known_MiHA_SNPs_Pos$CHROM_POS)
  if(length(ind) > 0){
    
    MiHA_table <- temp_tab[ind, c(1,2,4,5,11)]
    save(MiHA_table, file = paste0(output_fp, gsub(".txt", "_Restricted_MiHA.RData",file_list[id])))
    
  }
    
}