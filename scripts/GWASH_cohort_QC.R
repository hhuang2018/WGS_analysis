
#source("utils/util.R")
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(ggplot2)
GWAS_samp_table_fp <- "../ClinVar/GWASH/Metadata/newfnldataib1203.csv"
PLILNK_file_fp <- "../../2018WGS/Data/GWASH/unimputed/"

GWAS_sample_table <-read.csv(GWAS_samp_table_fp, header = T, stringsAsFactors = F)
sexCheck_table <- read.table(paste0(PLILNK_file_fp,'check_sex.sexcheck'), header = T, stringsAsFactors = F)

num2rid <- function(ID){
  
  strID <- as.character(ID)
  strID_list <- strsplit(strID, '')[[1]]
  
  if(7-nchar(strID) < 3){
    addZeros <- paste0(rep('0', 7-nchar(strID)), collapse = '')
    
  }else if(7-nchar(strID) == 3){
    addZeros <- paste0(rep('0', 7-nchar(strID)), collapse = '')
    addZeros <- paste0(addZeros, '-', collapse = '-')
    
  }else{
    addZeros <- ''
    for(idxx2 in 1:(7-nchar(strID))){
      
      if(idxx2 == 3 | idxx2 == 6){
        addZeros <- paste0(addZeros, '0-', collapse = '')
      }else{
        addZeros <- paste0(addZeros, '0', collapse = '')
      }
      
    }
    
  }
  
  
  new_rid <- ''
  for(idxx in nchar(strID):1){
    
    if(idxx == nchar(strID) - 1 | idxx == nchar(strID) - 4){
      new_rid <- paste0(strID_list[idxx], '-', new_rid, collapse = '')
    }else{
      new_rid <- paste0(strID_list[idxx], new_rid, collapse = '')
    }
    
  }
  new_rid <- paste0(addZeros, new_rid, collapse = '')
  
  return(new_rid)
}

num2did <- function(ID){
  
  total_num_len = 9
  block_digt_len = 4
  
  strID <- as.character(ID)
  strID_list <- strsplit(strID, '')[[1]]
  
  if(total_num_len-nchar(strID) < block_digt_len){
    addZeros <- paste0(rep('0', total_num_len-nchar(strID)), collapse = '')
    
  }else if(total_num_len-nchar(strID) == block_digt_len){
    addZeros <- paste0(rep('0', total_num_len-nchar(strID)), collapse = '')
    addZeros <- paste0(addZeros, '-', collapse = '-')
    
  }else{
    addZeros <- ''
    for(idxx2 in 1:(total_num_len-nchar(strID))){
      
      if(idxx2 == block_digt_len | idxx2 == block_digt_len*2){
        addZeros <- paste0(addZeros, '0-', collapse = '')
      }else{
        addZeros <- paste0(addZeros, '0', collapse = '')
      }
      
    }
    
  }
  
  new_did <- ''
  for(idxx in nchar(strID):1){
    
    if(idxx == nchar(strID) - 1 | idxx == nchar(strID) - block_digt_len-1){
      new_did <- paste0(strID_list[idxx], '-', new_did, collapse = '')
    }else{
      new_did <- paste0(strID_list[idxx], new_did, collapse = '')
    }
    
  }
  new_did <- paste0(addZeros, new_did, collapse = '')

  return(new_did)
}

## sex check
# Donor/Recipient sex matching: 1=Male/Male, 2=Male/Female, 3=Female/Male, 4=Female/Female.

num_case <- dim(sexCheck_table)[1]
sexCheck_table$metaSex <- integer(num_case)
sexCheck_table$inferedSex <- integer(num_case)
for(id in 1:num_case){
  
  IID <- gsub("D","",sexCheck_table$IID[id])
  
  if(nchar(IID) == 11){ # donor
    #did <- as.integer(gsub("-", "",sexCheck_table$IID[id]))
    did <- as.numeric(gsub("[^\\d]+", "", IID, perl=TRUE))
    ind <- which(GWAS_sample_table$did %in% did)
    if(length(ind) >= 1){
      if (GWAS_sample_table$sexmatch[ind] < 3){ # donor male=1
        sexCheck_table$metaSex[id] <- 1
      }else{ # donor female=2
        sexCheck_table$metaSex[id] <- 2
      }
    }else{
      cat(paste0(sexCheck_table$IID[id], " No such individual found!"))
    }
    
  }else if(nchar(IID) == 9){ # recipient
    rid <- as.integer(gsub("-", "", IID))
    ind2 <- which(GWAS_sample_table$rid %in% rid)
    if(length(ind2) >= 1){
      if(GWAS_sample_table$sexmatch[ind2] %% 2 == 1){ # recipient male=1
        
        sexCheck_table$metaSex[id] <- 1
        
      }else{ # recipient female=2
        sexCheck_table$metaSex[id] <- 2
      }
    }else{
      cat(paste0(sexCheck_table$IID[id], " No such individual found!"))
    }
  }
  
  if(sexCheck_table$F[id] <= 0.2){ # female=2
    sexCheck_table$inferedSex[id] <- 2
  }else if (sexCheck_table$F[id] >=0.8){ # male=1
    sexCheck_table$inferedSex[id] <- 1
  }
  
}

### 
discordance <- which(sexCheck_table$metaSex != sexCheck_table$SNPSEX)
sexCheck_table[discordance,] # n = 14 

discordance_inferred <- which(sexCheck_table$inferedSex != sexCheck_table$SNPSEX)
sexCheck_table[discordance_inferred,] # n = 0

discordance_inferred_meta <- which(sexCheck_table$inferedSex != sexCheck_table$metaSex)
sexCheck_table[discordance_inferred_meta,] # n = 14

#####
colnames(sexCheck_table) <- c('FID', 'IID', 'PEDSEX', 'SNPSEX', 'STATUS', 'F', 'BIO_SEX', 'INFERRED_F_SEX')
write.csv(sexCheck_table[discordance,], file = paste0(PLILNK_file_fp, 'sex_discordant_IDs.csv'), row.names = F)
write.table(sexCheck_table[discordance,], file = paste0(PLILNK_file_fp,'sex_discordant_IDs.txt'), row.names = F, quote = F)

#### 
# available cases in SNP data
####
num_meta_cases <- dim(GWAS_sample_table)[1]            # 1378 cases with metadata
GWAS_sample_table$pres_abs <- integer(num_meta_cases)
IDs_w_D <- integer(50)
counter <- 0

num_case <- dim(sexCheck_table)[1]
sexCheck_table$pres_abs <- integer(num_case)

for(id in 1:num_meta_cases){
  
  DID <- num2did(GWAS_sample_table$did[id])
  RID <- num2rid(GWAS_sample_table$rid[id])
  
  if(length(which(sexCheck_table$IID %in% DID)) > 0){ # donor exists - no 'D'
    
    if(length(which(sexCheck_table$IID %in% RID)) > 0){ # recipient exists - no 'D'
      
      GWAS_sample_table$pres_abs[id] <- 1

      GWAS_sample_table$did[id] <- DID
      GWAS_sample_table$rid[id] <- RID
      
      sexCheck_table$pres_abs[which(sexCheck_table$IID %in% DID)] <- 1
      sexCheck_table$pres_abs[which(sexCheck_table$IID %in% RID)] <- 1

    }else if(length(which(sexCheck_table$IID %in% paste0(RID, 'D'))) > 0){ # recipient exists - with 'D'
      
      GWAS_sample_table$pres_abs[id] <- 1
      
      GWAS_sample_table$did[id] <- DID
      GWAS_sample_table$rid[id] <- paste0(RID, 'D')
      
      counter <- counter + 1
      IDs_w_D[counter] <- paste0(RID, 'D')
      
      sexCheck_table$pres_abs[which(sexCheck_table$IID %in% DID)] <- 1
      sexCheck_table$pres_abs[which(sexCheck_table$IID %in% paste0(RID, 'D'))] <- 1
    }
    
  }else if(length(which(sexCheck_table$IID %in% paste0(DID, 'D'))) > 0){# donor exists - with 'D'
    
    if(length(which(sexCheck_table$IID %in% RID)) > 0){ # recipient exists - no 'D'
      
      GWAS_sample_table$pres_abs[id] <- 1
      
      GWAS_sample_table$did[id] <- paste0(DID, 'D')
      GWAS_sample_table$rid[id] <- RID
      
      counter <- counter + 1
      IDs_w_D[counter] <- paste0(DID, 'D')
      
      sexCheck_table$pres_abs[which(sexCheck_table$IID %in% paste0(DID, 'D'))] <- 1
      sexCheck_table$pres_abs[which(sexCheck_table$IID %in% RID)] <- 1
      
    }else if(length(which(sexCheck_table$IID %in% paste0(RID, 'D'))) > 0){ # recipient exists - with 'D'
      
      GWAS_sample_table$pres_abs[id] <- 1
      
      GWAS_sample_table$did[id] <- paste0(DID, 'D')
      GWAS_sample_table$rid[id] <- paste0(RID, 'D')
      
      counter <- counter + 1
      IDs_w_D[counter] <- paste0(DID, 'D')
      counter <- counter + 1
      IDs_w_D[counter] <- paste0(RID, 'D')
      
      sexCheck_table$pres_abs[which(sexCheck_table$IID %in% paste0(DID, 'D'))] <- 1
      sexCheck_table$pres_abs[which(sexCheck_table$IID %in% paste0(RID, 'D'))] <- 1
    }
  }
  
}

cat(IDs_w_D[1:counter], file = paste0(PLILNK_file_fp,'irregular_IDs.txt'), sep = '\n')

exclude_IDs <- sexCheck_table[sexCheck_table$pres_abs == 0, c('FID', 'IID')] # 7

write.table(exclude_IDs, file = paste0(PLILNK_file_fp,'unmatched_IDs.txt'), row.names = F, quote = F)

Avail_cases_table <- GWAS_sample_table[which(GWAS_sample_table$pres_abs == 1), ] # 988
# Total: 988 cases (987 pairs with available data)
# AML  - 332 cases 
# ALL  - 123 cases 
# CML  - 343 cases
# MDS  - 191 cases -- could be considered as AML

write.csv(Avail_cases_table, file = paste0(PLILNK_file_fp,'GWAS_avaialbe_cases_info.csv'), row.names = F)

# AML_cases_table <- Avail_cases_table[which(Avail_cases_table$disease == 10), ] # 332

## donors only IDs
donor_num <- length(Avail_cases_table$did)
donor_IDs <- data.frame(FID = integer(donor_num),
                        IID = character(donor_num),
                        stringsAsFactors = F)
for(id in 1:donor_num){
  
  index1 <- which(sexCheck_table$IID %in% Avail_cases_table$did[id])
  if(length(index1) == 1){
    
    donor_IDs$FID[id] <- sexCheck_table$FID[index1]
    donor_IDs$IID[id] <- sexCheck_table$IID[index1]
    
  }
}

write.table(donor_IDs, file = paste0(PLILNK_file_fp,'donor_IDs.txt'), row.names = F, quote = F, col.names = F)

## all availabel BMTcase individual IDs
individual_num <- length(Avail_cases_table$did) * 2
individual_IDs <- data.frame(FID = integer(individual_num),
                             IID = character(individual_num),
                             stringsAsFactors = F)
counter <- 0
for(id in 1:individual_num){
  
  index1 <- which(sexCheck_table$IID %in% Avail_cases_table$did[id])
  if(length(index1) == 1){
    counter <- counter + 1
    
    individual_IDs$FID[counter] <- sexCheck_table$FID[index1]
    individual_IDs$IID[counter] <- sexCheck_table$IID[index1]
    
  }
  
  index2 <- which(sexCheck_table$IID %in% Avail_cases_table$rid[id])
  if(length(index2) == 1){
    counter <- counter + 1
    
    individual_IDs$FID[counter] <- sexCheck_table$FID[index2]
    individual_IDs$IID[counter] <- sexCheck_table$IID[index2]
    
  }
  
}

write.table(individual_IDs, file = paste0(PLILNK_file_fp,'Available_individual_IDs.txt'), row.names = F, quote = F, col.names = F)

## remove Sex discordant individuals
sexDiscordant <- sexCheck_table[discordance,]
num_discordant <- dim(sexDiscordant)[1]

sexDiscordant$avail <- integer(num_discordant)

for(id in 1:num_discordant){
  
  if(length(which(individual_IDs$IID %in% sexDiscordant$IID[id]))){
    
    sexDiscordant$avail <- 1
    
  }
  
}

num_individual <- dim(individual_IDs)[1]
avail_individuals <- individual_IDs
avail_individuals$rm <- integer(num_individual)
for(id in 1:num_individual){
  
  if(length(which(sexDiscordant$IID %in% avail_individuals$IID[id]))>0){
    
    avail_individuals$rm[id] <- 1 
    
    if(nchar(avail_individuals$IID[id])>=11){ # donor
      cat('---------------------\n')
      cat(paste0('Removing donor: ', avail_individuals$IID[id], '\n'))
      matching_rid <- Avail_cases_table$rid[which(Avail_cases_table$did %in% avail_individuals$IID[id])]
      
      avail_individuals$rm[which(avail_individuals$IID %in% matching_rid)] <- 1
      cat(paste0('Removing matching recipient: ', avail_individuals$IID[which(avail_individuals$IID %in% matching_rid)], '\n'))
    }else{ # recipient
      cat('---------------------\n')
      cat(paste0('Removing Recipient: ', avail_individuals$IID[id], '\n'))
      matching_did <- Avail_cases_table$did[which(Avail_cases_table$rid %in% avail_individuals$IID[id])]
      cat(paste0('Removing matching donor: ', avail_individuals$IID[id], '\n'))
      avail_individuals$rm[which(avail_individuals$IID %in% matching_did)] <- 1
      
    }
    
  }
  
}

remaining_individuals <- avail_individuals[which(avail_individuals$rm==0),] # 1948 are remaining

write.table(remaining_individuals, file = paste0(PLILNK_file_fp,'Available_individual_IDs_clean.txt'), row.names = F, quote = F, col.names = F)

remove_individuals <- avail_individuals[which(avail_individuals$rm==1),] # 28 individuals will be removed
write.table(remove_individuals, file = paste0(PLILNK_file_fp,'exclude_avail_pairs_ambiguous_sex.txt'), row.names = F, quote = F, col.names = F)

remaining_donors <- remaining_individuals[which(sapply(1:dim(remaining_individuals)[1], function(x) nchar(remaining_individuals$IID[x])>=11)), ]
write.table(remaining_donors, file = paste0(PLILNK_file_fp,'available_donors.txt'), row.names = F, quote = F, col.names = F)

### cleaned_case_metainfo
Avail_cases_table$pres_abs <- integer(dim(Avail_cases_table)[1])
num_individuals <- dim(remaining_individuals)[1]
for(id in 1:num_individuals){
  
  if(nchar(remaining_individuals$IID[id])>=11){
    d_ind <- which(Avail_cases_table$did %in% remaining_individuals$IID[id])
    if(length(d_ind)> 0){
      if(length(which(remaining_individuals$IID %in% Avail_cases_table$rid[d_ind]))>0){
        Avail_cases_table$pres_abs[d_ind] <- 1
      }else{
        cat(paste0('No Matching recipient present for ',remaining_individuals$IID[id] ,' !\n'))
      }
      
    }
  }else{
    r_ind <- which(Avail_cases_table$rid %in% remaining_individuals$IID[id])
    if(length(r_ind)>0){
      if(length(which(remaining_individuals$IID %in% Avail_cases_table$did[r_ind]))>0){
        Avail_cases_table$pres_abs[r_ind] <- 1
      }else{
        cat(paste0('No Matching donor present for ',remaining_individuals$IID[id] ,' !\n'))
      }
    }
  }
  
}

cleaned_BMT_case_table <- Avail_cases_table[which(Avail_cases_table$pres_abs==1),]
write.csv(cleaned_BMT_case_table, file = paste0(PLILNK_file_fp,'GWAS_cleaned_BMT_cases_info.csv'), row.names = F)
###########
# HWE stats
###########
HWE_stats <- read.table(paste0(PLILNK_file_fp,'hardy_stats.hwe'), header = T, stringsAsFactors = F)
HWE_ordered <- HWE_stats[order(HWE_stats$P),]
length(which(HWE_ordered$P<0.001))

donors_only_HWE <- read.table(paste0(PLILNK_file_fp,'available_donors_only_hardy.hwe'), header = T, stringsAsFactors = F)
donors_only_HWE_ordered <- donors_only_HWE[order(donors_only_HWE$P),]
length(which(donors_only_HWE$P<=0.001)) # 23,730 # available clean donors - 21,866

rm_SNPs_hwe_001 <- donors_only_HWE$SNP[which(donors_only_HWE$P<=0.001)]
write.table(rm_SNPs_hwe_001, file = paste0(PLILNK_file_fp,'exclude_cleaned_SNPs_hwe_001.txt'), row.names = F, quote = F, col.names = F)

length(which(donors_only_HWE$P<=1e-6))  # 17,864 # available clean donors - 17,344
rm_SNPs_hwe_1e6 <- donors_only_HWE$SNP[which(donors_only_HWE$P<=1e-6)]
write.table(rm_SNPs_hwe_1e6, file = paste0(PLILNK_file_fp,'exclude_cleaned_SNPs_hwe_1e_6.txt'), row.names = F, quote = F, col.names = F)


length(which(donors_only_HWE$P<1e-4))  # 20,240 # available clean donors - 19,029
rm_SNPs_hwe_1e4 <- donors_only_HWE$SNP[which(donors_only_HWE$P<=1e-4)]
write.table(rm_SNPs_hwe_1e4, file = paste0(PLILNK_file_fp,'exclude_cleaned_SNPs_hwe_1e_4.txt'), row.names = F, quote = F, col.names = F)

