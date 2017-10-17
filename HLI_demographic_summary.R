load("../Data/ID_table_wCaseID.RData")
HLI_metadata <- read.csv("../HLI_hla_mg_v3.csv", stringsAsFactors = F)
HLI_cohort <- read.csv("../Data/HLI-pull.csv")

library(eeptools)
recipient_dob <- as.Date(HLI_cohort$Rbdte)
transplant_dt <- as.Date(HLI_cohort$transplant_dte)
donor_dob <- as.Date(HLI_cohort$Dbdte)

recipient_ages <- floor(age_calc(recipient_dob, enddate = transplant_dt, units = "years"))
donor_ages <- floor(age_calc(donor_dob, enddate = transplant_dt, units = "years"))

HLI_cohort$Rage <- recipient_ages
HLI_cohort$Dage <- donor_ages
HLI_cohort$tx_year <- as.numeric(format(transplant_dt, "%Y"))

##### Add BMTcase number to the ID table
# Recipient_ids <- which(ID_table$subjectType == "R")
# Donor_ids <- which(ID_table$subjectType == "D")
# num_Rcases <- length(Recipient_ids)
# num_Dcases <- length(Donor_ids)
# ID_table$caseID <- 0
# for(id in 1:num_Rcases){
#    
#    caseID1 <- HLI_cohort$pair_id[which(HLI_cohort$rid %in% ID_table$R_D_ID[Recipient_ids[id]])]
#    caseID2 <- HLI_metadata$bmt_case_num[which(HLI_metadata$nmdp_rid %in% ID_table$R_D_ID[Recipient_ids[id]])]
#    
#    if(length(caseID1) == 1){
#      
#      if(caseID1 == caseID2){
#        
#        ID_table$caseID[Recipient_ids[id]] <- caseID1
#        
#      }else cat("Recipient ", id ,"- Inconsistent IDs between HLI metadata and HLI cohort.\n")
#    }else cat("Recipient ", id ,"- Multiple recipient cases.\n")
# } # Recipient BMTcaseID
#  
# for(id in 1:num_Dcases){
#   
#   caseID1 <- HLI_cohort$pair_id[which(HLI_cohort$did %in% ID_table$R_D_ID[Donor_ids[id]])]
#   caseID2 <- HLI_metadata$bmt_case_num[which(HLI_metadata$nmdp_id %in% ID_table$R_D_ID[Donor_ids[id]])]
#   
#   if(length(caseID1) == 1){
#     
#     if(caseID1 == caseID2){
#       
#       ID_table$caseID[Donor_ids[id]] <- caseID1
#       
#     }else cat("Donor ", id ,"- Inconsistent IDs between HLI metadata and HLI cohort.\n")
#   }else cat("Donor ", id ,"- Multiple Donor cases.\n") # id = 70 is the multiple donor case
# } # Donor BMTcaseID

# write.csv(ID_table, file = "../Data/ID_table_wCaseID.csv", row.names = F)
# save(ID_table, file = "../Data/ID_table_wCaseID.RData")
################
groups_stats <- as.data.frame(table(ID_table$GroupID))
paired_ID_table <- ID_table[which(ID_table$GroupID %in% groups_stats$Var1[which(groups_stats$Freq == 2)]), ]
paired_case_num <- unique(paired_ID_table$caseID)

HLI_paired_cohort <- HLI_cohort[which(HLI_cohort$pair_id %in% paired_case_num),]

HLI_paired_cohort$Group <- sapply(1:length(HLI_paired_cohort$pair_id), function(x) unique(paired_ID_table$Group[which(paired_ID_table$caseID %in% HLI_paired_cohort$pair_id[x])]))

HLI_paired_cohort$R_sex <- sapply(1:length(HLI_paired_cohort$pair_id), function(x) HLI_metadata$rid_sex[which(HLI_metadata$bmt_case_num %in% HLI_paired_cohort$pair_id[x])])
HLI_paired_cohort$D_sex <- sapply(1:length(HLI_paired_cohort$pair_id), function(x) HLI_metadata$donor_sex[which(HLI_metadata$bmt_case_num %in% HLI_paired_cohort$pair_id[x])])

# for(id in 1:length(HLI_paired_cohort$pair_id)){
#   
#   if(unique(paired_ID_table$Group[which(paired_ID_table$caseNumber %in% HLI_paired_cohort$pair_id[id])]) != HLI_paired_cohort$Group[id]) cat(id)
#   
# }

#### recipient age stats 
recipient_age_stats <- as.data.frame(table(HLI_paired_cohort[c("Group", "Rage")]))
recipient_age_stats$Rage <- as.character(recipient_age_stats$Rage) 
recipient_age_stats$Rage <- as.numeric(recipient_age_stats$Rage)
recipient_age_stats_table <- t(data.frame("15" = vector(mode = "numeric", length = 2),
                                        "30" = vector(mode = "numeric", length = 2),
                                        "45" = vector(mode = "numeric", length = 2),
                                        "60" = vector(mode = "numeric", length = 2),
                                        "61" = vector(mode = "numeric", length = 2)))
colnames(recipient_age_stats_table) <- c("aGVHD", "nonGVHD")
rownames(recipient_age_stats_table) <- c("0~17", "18~30", "31~45", "45~60", "61+")

aGVHD_index <- which(recipient_age_stats$Group == "a")
nGVHD_index <- which(recipient_age_stats$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  if(recipient_age_stats$Rage[aGVHD_index[id]] <= 17){
    
    recipient_age_stats_table[1, 1] <-  recipient_age_stats_table[1, 1] + recipient_age_stats$Freq[aGVHD_index[id]]
  
  }else if (recipient_age_stats$Rage[aGVHD_index[id]] <= 30){
    
    recipient_age_stats_table[2, 1] <-  recipient_age_stats_table[2, 1] + recipient_age_stats$Freq[aGVHD_index[id]]
    
  }else if (recipient_age_stats$Rage[aGVHD_index[id]] <= 45){
    
    recipient_age_stats_table[3, 1] <-  recipient_age_stats_table[3, 1] + recipient_age_stats$Freq[aGVHD_index[id]]
    
  }else if (recipient_age_stats$Rage[aGVHD_index[id]] <= 60){
    
    recipient_age_stats_table[4, 1] <-  recipient_age_stats_table[4, 1] + recipient_age_stats$Freq[aGVHD_index[id]]
    
  }else{
    
    recipient_age_stats_table[5, 1] <-  recipient_age_stats_table[5, 1] + recipient_age_stats$Freq[aGVHD_index[id]]
    
  }
  
  ##### nonGVHD groups
  if(recipient_age_stats$Rage[nGVHD_index[id]] <= 17){
    
    recipient_age_stats_table[1, 2] <-  recipient_age_stats_table[1, 2] + recipient_age_stats$Freq[nGVHD_index[id]]
    
  }else if (recipient_age_stats$Rage[nGVHD_index[id]] <= 30){
    
    recipient_age_stats_table[2, 2] <-  recipient_age_stats_table[2, 2] + recipient_age_stats$Freq[nGVHD_index[id]]
    
  }else if (recipient_age_stats$Rage[nGVHD_index[id]] <= 45){
    
    recipient_age_stats_table[3, 2] <-  recipient_age_stats_table[3, 2] + recipient_age_stats$Freq[nGVHD_index[id]]
    
  }else if (recipient_age_stats$Rage[nGVHD_index[id]] <= 60){
    
    recipient_age_stats_table[4, 2] <-  recipient_age_stats_table[4, 2] + recipient_age_stats$Freq[nGVHD_index[id]]
    
  }else{
    
    recipient_age_stats_table[5, 2] <-  recipient_age_stats_table[5, 2] + recipient_age_stats$Freq[nGVHD_index[id]]
    
  }
  
}

#### donor age stats
donor_age_stats <- as.data.frame(table(HLI_paired_cohort[c("Group", "Dage")]))
donor_age_stats$Dage <- as.character(donor_age_stats$Dage) 
donor_age_stats$Dage <- as.numeric(donor_age_stats$Dage)
donor_age_stats_table <- t(data.frame("15" = vector(mode = "numeric", length = 2),
                                          "30" = vector(mode = "numeric", length = 2),
                                          "45" = vector(mode = "numeric", length = 2),
                                          "60" = vector(mode = "numeric", length = 2),
                                          "61" = vector(mode = "numeric", length = 2)))
colnames(donor_age_stats_table) <- c("aGVHD", "nonGVHD")
rownames(donor_age_stats_table) <- c("0~17", "18~30", "31~45", "45~60", "61+")

aGVHD_index <- which(donor_age_stats$Group == "a")
nGVHD_index <- which(donor_age_stats$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  if(donor_age_stats$Dage[aGVHD_index[id]] <= 17){
    
    donor_age_stats_table[1, 1] <-  donor_age_stats_table[1, 1] + donor_age_stats$Freq[aGVHD_index[id]]
    
  }else if (donor_age_stats$Dage[aGVHD_index[id]] <= 30){
    
    donor_age_stats_table[2, 1] <-  donor_age_stats_table[2, 1] + donor_age_stats$Freq[aGVHD_index[id]]
    
  }else if (donor_age_stats$Dage[aGVHD_index[id]] <= 45){
    
    donor_age_stats_table[3, 1] <-  donor_age_stats_table[3, 1] + donor_age_stats$Freq[aGVHD_index[id]]
    
  }else if (donor_age_stats$Dage[aGVHD_index[id]] <= 60){
    
    donor_age_stats_table[4, 1] <-  donor_age_stats_table[4, 1] + donor_age_stats$Freq[aGVHD_index[id]]
    
  }else{
    
    donor_age_stats_table[5, 1] <-  donor_age_stats_table[5, 1] + donor_age_stats$Freq[aGVHD_index[id]]
    
  }
  
  ##### nonGVHD groups
  if(donor_age_stats$Dage[nGVHD_index[id]] <= 17){
    
    donor_age_stats_table[1, 2] <-  donor_age_stats_table[1, 2] + donor_age_stats$Freq[nGVHD_index[id]]
    
  }else if (donor_age_stats$Dage[nGVHD_index[id]] <= 30){
    
    donor_age_stats_table[2, 2] <-  donor_age_stats_table[2, 2] + donor_age_stats$Freq[nGVHD_index[id]]
    
  }else if (donor_age_stats$Dage[nGVHD_index[id]] <= 45){
    
    donor_age_stats_table[3, 2] <-  donor_age_stats_table[3, 2] + donor_age_stats$Freq[nGVHD_index[id]]
    
  }else if (donor_age_stats$Dage[nGVHD_index[id]] <= 60){
    
    donor_age_stats_table[4, 2] <-  donor_age_stats_table[4, 2] + donor_age_stats$Freq[nGVHD_index[id]]
    
  }else{
    
    donor_age_stats_table[5, 2] <-  donor_age_stats_table[5, 2] + donor_age_stats$Freq[nGVHD_index[id]]
    
  }
  
}


##### recipient disease type
recipient_disease_stats <- as.data.frame(table(HLI_paired_cohort[c("Group", "disease")]))
recipient_disease_stats$disease <- as.character(recipient_disease_stats$disease) 
recipient_disease_stats_table <- t(data.frame("AML" = vector(mode = "numeric", length = 2),
                                          "ALL" = vector(mode = "numeric", length = 2),
                                          "MDS" = vector(mode = "numeric", length = 2),
                                          "OL" = vector(mode = "numeric", length = 2)))
colnames(recipient_disease_stats_table) <- c("aGVHD", "nonGVHD")

aGVHD_index <- which(recipient_disease_stats$Group == "a")
nGVHD_index <- which(recipient_disease_stats$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  switch(recipient_disease_stats$disease[aGVHD_index[id]],
         
         AML = {recipient_disease_stats_table[1, 1] <- recipient_disease_stats$Freq[aGVHD_index[id]]},
         ALL = {recipient_disease_stats_table[2, 1] <- recipient_disease_stats$Freq[aGVHD_index[id]]},
         MDS = {recipient_disease_stats_table[3, 1] <- recipient_disease_stats$Freq[aGVHD_index[id]]},
         {recipient_disease_stats_table[4, 1] <- recipient_disease_stats$Freq[aGVHD_index[id]]}
         )
    
  ##### nonGVHD groups
  switch(recipient_disease_stats$disease[nGVHD_index[id]],
         
         AML = {recipient_disease_stats_table[1, 2] <- recipient_disease_stats$Freq[nGVHD_index[id]]},
         ALL = {recipient_disease_stats_table[2, 2] <- recipient_disease_stats$Freq[nGVHD_index[id]]},
         MDS = {recipient_disease_stats_table[3, 2] <- recipient_disease_stats$Freq[nGVHD_index[id]]},
         {recipient_disease_stats_table[4, 2] <- recipient_disease_stats$Freq[nGVHD_index[id]]}
  )
  
}

######## Transplant year suumary 
recipient_transplant_year <- as.data.frame(table(HLI_paired_cohort[c("Group", "tx_year")]))
recipient_transplant_year$tx_year <- as.character(recipient_transplant_year$tx_year)
recipient_tx_year_table <- t(data.frame("2000" = vector(mode = "numeric", length = 2),
                                              "2001" = vector(mode = "numeric", length = 2),
                                              "2002" = vector(mode = "numeric", length = 2),
                                              "2003" = vector(mode = "numeric", length = 2),
                                              "2004" = vector(mode = "numeric", length = 2),
                                              "2005" = vector(mode = "numeric", length = 2),
                                              "2006" = vector(mode = "numeric", length = 2),
                                              "2007" = vector(mode = "numeric", length = 2),
                                              "2008" = vector(mode = "numeric", length = 2),
                                              "2009" = vector(mode = "numeric", length = 2),
                                              "2010" = vector(mode = "numeric", length = 2),
                                              "2011" = vector(mode = "numeric", length = 2)))
colnames(recipient_tx_year_table) <- c("aGVHD", "nonGVHD")

aGVHD_index <- which(recipient_transplant_year$Group == "a")
nGVHD_index <- which(recipient_transplant_year$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  switch(recipient_transplant_year$tx_year[aGVHD_index[id]],
         "2000" = {recipient_tx_year_table[1, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2001" = {recipient_tx_year_table[2, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2002" = {recipient_tx_year_table[3, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2003" = {recipient_tx_year_table[4, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2004" = {recipient_tx_year_table[5, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2005" = {recipient_tx_year_table[6, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2006" = {recipient_tx_year_table[7, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2007" = {recipient_tx_year_table[8, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2008" = {recipient_tx_year_table[9, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2009" = {recipient_tx_year_table[10, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2010" = {recipient_tx_year_table[11, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]},
         "2011" = {recipient_tx_year_table[12, 1] <- recipient_transplant_year$Freq[aGVHD_index[id]]}
  )

  ##### nonGVHD groups
  switch(recipient_transplant_year$tx_year[nGVHD_index[id]],
         "2000" = {recipient_tx_year_table[1, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2001" = {recipient_tx_year_table[2, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2002" = {recipient_tx_year_table[3, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2003" = {recipient_tx_year_table[4, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2004" = {recipient_tx_year_table[5, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2005" = {recipient_tx_year_table[6, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2006" = {recipient_tx_year_table[7, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2007" = {recipient_tx_year_table[8, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2008" = {recipient_tx_year_table[9, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2009" = {recipient_tx_year_table[10, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2010" = {recipient_tx_year_table[11, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]},
         "2011" = {recipient_tx_year_table[12, 2] <- recipient_transplant_year$Freq[nGVHD_index[id]]}
  )
  
}

write.csv(recipient_tx_year_table, file="../FirstPaper/Table/tx_year.csv")

1####### Ethnicity ##################################################
# Broad Race - AFA; API; CAU; HIS; NAM
num_cases <- dim(HLI_paired_cohort)[1]
HLI_paired_cohort$R_BroadRace <- character(num_cases)
HLI_paired_cohort$D_BroadRace <- character(num_cases)

AFA <- c("AAFA", "AFB", "CARB", "SCAMB", "AFA")
API <- c("AINDI", "FILII", "HAWI", "JAPI", "KORI", "NCHI", "SCSEAI", "VIET", "API")
CAU <- c("EURCAU", "MENAFC", "CAU", "EURWRC", "MEDIT", "EEURO", "NEURO", "WEURO", "UNK", "DEC", "EURWRC", "MULTI", "MCAU", "WSCA","WCARIB", "NAMER")
HIS <- c("MSWHIS", "SCAHIS", "CARHIS", "HIS")
NAM <- c("CARIBI", "AMIND", "AISC", "ALANAM", "NAM")
for(id in 1:num_cases){
  
  ### Recipient 
  if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Rprim_race_cde[id])), CAU)){
    
    HLI_paired_cohort$R_BroadRace[id] <- "CAU"
    
  }else if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Rprim_race_cde[id])), AFA)){
    
    HLI_paired_cohort$R_BroadRace[id] <- "AFA"
    
  }else if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Rprim_race_cde[id])), API)){
    
    HLI_paired_cohort$R_BroadRace[id] <- "API"
    
  }else if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Rprim_race_cde[id])), HIS)){
    
    HLI_paired_cohort$R_BroadRace[id] <- "HIS"
    
  }else if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Rprim_race_cde[id])), NAM)){
    
    HLI_paired_cohort$R_BroadRace[id] <- "NAM"
    
  }
  
  ### Donor 
  if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Dprim_race_cde[id])), CAU)){
    
    HLI_paired_cohort$D_BroadRace[id] <- "CAU"
    
  }else if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Dprim_race_cde[id])), AFA)){
    
    HLI_paired_cohort$D_BroadRace[id] <- "AFA"
    
  }else if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Dprim_race_cde[id])), API)){
    
    HLI_paired_cohort$D_BroadRace[id] <- "API"
    
  }else if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Dprim_race_cde[id])), HIS)){
    
    HLI_paired_cohort$D_BroadRace[id] <- "HIS"
    
  }else if(is.element(gsub(" ", "", as.character(HLI_paired_cohort$Dprim_race_cde[id])), NAM)){
    
    HLI_paired_cohort$D_BroadRace[id] <- "NAM"
    
  }
}

HLI_paired_cohort$D2Rrace <- paste0(HLI_paired_cohort$D_BroadRace, ">", HLI_paired_cohort$R_BroadRace)
HLI_SIRE <- as.data.frame(table(HLI_paired_cohort[c("Group", "D2Rrace")]))

HLI_SIRE$D2Rrace <- as.character(HLI_SIRE$D2Rrace) 
HLI_SIRE_stats_table <- t(data.frame("CAU>CAU" = vector(mode = "numeric", length = 2),
                                     "CAU>NAM" = vector(mode = "numeric", length = 2),
                                     "NAM>CAU" = vector(mode = "numeric", length = 2),
                                     "HIS>CAU" = vector(mode = "numeric", length = 2),
                                     "NAM>NAM" = vector(mode = "numeric", length = 2),
                                     "API>API" = vector(mode = "numeric", length = 2),
                                     "CAU>AFA" = vector(mode = "numeric", length = 2),
                                     "AFA>CAU" = vector(mode = "numeric", length = 2),
                                     "AFA>AFA" = vector(mode = "numeric", length = 2),
                                     "API>CAU" = vector(mode = "numeric", length = 2),
                                     "HIS>AFA" = vector(mode = "numeric", length = 2)))
colnames(HLI_SIRE_stats_table) <- c("aGVHD", "nonGVHD")
rownames(HLI_SIRE_stats_table) <- c("CAU>CAU", "CAU>NAM", "NAM>CAU", "HIS>CAU", "NAM>NAM", "API>API",
                                          "CAU>AFA", "AFA>CAU", "AFA>AFA", "API>CAU", "HIS>AFA")

aGVHD_index <- which(HLI_SIRE$Group == "a")
nGVHD_index <- which(HLI_SIRE$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  switch(HLI_SIRE$D2Rrace[aGVHD_index[id]],
         
         `CAU>CAU` = {HLI_SIRE_stats_table[1, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `CAU>NAM` = {HLI_SIRE_stats_table[2, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `NAM>CAU` = {HLI_SIRE_stats_table[3, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `HIS>CAU` = {HLI_SIRE_stats_table[4, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `NAM>NAM` = {HLI_SIRE_stats_table[5, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `API>API` = {HLI_SIRE_stats_table[6, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `CAU>AFA` = {HLI_SIRE_stats_table[7, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `AFA>CAU` = {HLI_SIRE_stats_table[8, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `AFA>AFA` = {HLI_SIRE_stats_table[9, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `API>CAU` = {HLI_SIRE_stats_table[10, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]},
         `HIS>AFA` = {HLI_SIRE_stats_table[11, 1] <- HLI_SIRE$Freq[aGVHD_index[id]]}
  )
  ##### nonGVHD groups
  switch(HLI_SIRE$D2Rrace[nGVHD_index[id]],
         
         `CAU>CAU` = {HLI_SIRE_stats_table[1, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `CAU>NAM` = {HLI_SIRE_stats_table[2, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `NAM>CAU` = {HLI_SIRE_stats_table[3, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `HIS>CAU` = {HLI_SIRE_stats_table[4, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `NAM>NAM` = {HLI_SIRE_stats_table[6, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `API>API` = {HLI_SIRE_stats_table[7, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `CAU>AFA` = {HLI_SIRE_stats_table[8, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `AFA>CAU` = {HLI_SIRE_stats_table[5, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `AFA>AFA` = {HLI_SIRE_stats_table[9, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `API>CAU` = {HLI_SIRE_stats_table[10, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]},
         `HIS>AFA` = {HLI_SIRE_stats_table[11, 2] <- HLI_SIRE$Freq[nGVHD_index[id]]}
  )

}


# individual race
#### recipient SIRE stats  HLI_SIRE
recipient_SIRE_stats_ <- as.data.frame(table(HLI_paired_cohort[c("Group", "R_BroadRace")])) 
recipient_SIRE_stats_$R_BroadRace <- as.character(recipient_SIRE_stats_$R_BroadRace) 
recipient_SIRE_stats_table <- t(data.frame("CAU" = vector(mode = "numeric", length = 2),
                                           "NAM" = vector(mode = "numeric", length = 2),
                                           "AFA" = vector(mode = "numeric", length = 2),
                                           "API" = vector(mode = "numeric", length = 2),
                                           "HIS" = vector(mode = "numeric", length = 2)))

colnames(recipient_SIRE_stats_table) <- c("aGVHD", "nonGVHD")

aGVHD_index <- which(recipient_SIRE_stats$Group == "a")
nGVHD_index <- which(recipient_SIRE_stats$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  switch(recipient_SIRE_stats$R_BroadRace[aGVHD_index[id]],
         CAU = {recipient_SIRE_stats_table[1, 1] <- recipient_SIRE_stats$Freq[aGVHD_index[id]]},
         NAM = {recipient_SIRE_stats_table[2, 1] <- recipient_SIRE_stats$Freq[aGVHD_index[id]]},
         AFA = {recipient_SIRE_stats_table[3, 1] <- recipient_SIRE_stats$Freq[aGVHD_index[id]]},
         API = {recipient_SIRE_stats_table[4, 1] <- recipient_SIRE_stats$Freq[aGVHD_index[id]]},
         HIS = {recipient_SIRE_stats_table[5, 1] <- recipient_SIRE_stats$Freq[aGVHD_index[id]]})
  
  ##### nonGVHD groups
  switch(recipient_SIRE_stats$R_BroadRace[nGVHD_index[id]],
         CAU = {recipient_SIRE_stats_table[1, 2] <- recipient_SIRE_stats$Freq[nGVHD_index[id]]},
         NAM = {recipient_SIRE_stats_table[2, 2] <- recipient_SIRE_stats$Freq[nGVHD_index[id]]},
         AFA = {recipient_SIRE_stats_table[3, 2] <- recipient_SIRE_stats$Freq[nGVHD_index[id]]},
         API = {recipient_SIRE_stats_table[4, 2] <- recipient_SIRE_stats$Freq[nGVHD_index[id]]},
         HIS = {recipient_SIRE_stats_table[5, 2] <- recipient_SIRE_stats$Freq[nGVHD_index[id]]})
  
}

#### donor SIRE stats
donor_SIRE_stats <- as.data.frame(table(HLI_paired_cohort[c("Group", "D_BroadRace")])) 
donor_SIRE_stats$D_BroadRace <- as.character(donor_SIRE_stats$D_BroadRace) 
donor_SIRE_stats_table <- t(data.frame("CAU" = vector(mode = "numeric", length = 2),
                                       "NAM" = vector(mode = "numeric", length = 2),
                                       "AFA" = vector(mode = "numeric", length = 2),
                                       "API" = vector(mode = "numeric", length = 2),
                                       "HIS" = vector(mode = "numeric", length = 2)))

colnames(donor_SIRE_stats_table) <- c("aGVHD", "nonGVHD")

aGVHD_index <- which(donor_SIRE_stats$Group == "a")
nGVHD_index <- which(donor_SIRE_stats$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  switch(donor_SIRE_stats$D_BroadRace[aGVHD_index[id]],
         CAU = {donor_SIRE_stats_table[1, 1] <- donor_SIRE_stats$Freq[aGVHD_index[id]]},
         NAM = {donor_SIRE_stats_table[2, 1] <- donor_SIRE_stats$Freq[aGVHD_index[id]]},
         AFA = {donor_SIRE_stats_table[3, 1] <- donor_SIRE_stats$Freq[aGVHD_index[id]]},
         API = {donor_SIRE_stats_table[4, 1] <- donor_SIRE_stats$Freq[aGVHD_index[id]]},
         HIS = {donor_SIRE_stats_table[5, 1] <- donor_SIRE_stats$Freq[aGVHD_index[id]]})
  
  ##### nonGVHD groups
  switch(donor_SIRE_stats$D_BroadRace[nGVHD_index[id]],
         CAU = {donor_SIRE_stats_table[1, 2] <- donor_SIRE_stats$Freq[nGVHD_index[id]]},
         NAM = {donor_SIRE_stats_table[2, 2] <- donor_SIRE_stats$Freq[nGVHD_index[id]]},
         AFA = {donor_SIRE_stats_table[3, 2] <- donor_SIRE_stats$Freq[nGVHD_index[id]]},
         API = {donor_SIRE_stats_table[4, 2] <- donor_SIRE_stats$Freq[nGVHD_index[id]]}, 
         HIS = {donor_SIRE_stats_table[5, 2] <- donor_SIRE_stats$Freq[nGVHD_index[id]]})
  
}



##############################
## Sex mismatch
HLI_paired_cohort$D2RSex <- paste0(HLI_paired_cohort$D_sex, ">", HLI_paired_cohort$R_sex)
HLI_Sex_match_table <- as.data.frame(table(HLI_paired_cohort[c("Group", "D2RSex")]))
HLI_Sex_match_table$D2RSex <- as.character(HLI_Sex_match_table$D2RSex) 
HLI_Sex_match_table_stats <- t(data.frame("Female>Female" = vector(mode = "numeric", length = 2),
                                          "Female>Male" = vector(mode = "numeric", length = 2),
                                          "Male>Female" = vector(mode = "numeric", length = 2),
                                          "Male>Male" = vector(mode = "numeric", length = 2)))
colnames(HLI_Sex_match_table_stats) <- c("aGVHD", "nonGVHD")
rownames(HLI_Sex_match_table_stats) <- c("Female>Female", "Female>Male", "Male>Female", "Male>Male")

aGVHD_index <- which(HLI_Sex_match_table$Group == "a")
nGVHD_index <- which(HLI_Sex_match_table$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  switch(HLI_Sex_match_table$D2RSex[aGVHD_index[id]],
         
         "F>F" = {HLI_Sex_match_table_stats[1, 1] <- HLI_Sex_match_table$Freq[aGVHD_index[id]]},
         "F>M" = {HLI_Sex_match_table_stats[2, 1] <- HLI_Sex_match_table$Freq[aGVHD_index[id]]},
         "M>F" = {HLI_Sex_match_table_stats[3, 1] <- HLI_Sex_match_table$Freq[aGVHD_index[id]]},
         "M>M" = {HLI_Sex_match_table_stats[4, 1] <- HLI_Sex_match_table$Freq[aGVHD_index[id]]}
  )
  ##### nonGVHD groups
  switch(HLI_Sex_match_table$D2RSex[nGVHD_index[id]],
         
         "F>F" = {HLI_Sex_match_table_stats[1, 2] <- HLI_Sex_match_table$Freq[nGVHD_index[id]]},
         "F>M" = {HLI_Sex_match_table_stats[2, 2] <- HLI_Sex_match_table$Freq[nGVHD_index[id]]},
         "M>F" = {HLI_Sex_match_table_stats[3, 2] <- HLI_Sex_match_table$Freq[nGVHD_index[id]]},
         "M>M" = {HLI_Sex_match_table_stats[4, 2] <- HLI_Sex_match_table$Freq[nGVHD_index[id]]}
  )
  
}


# individual SEX
#### recipient Sex stats  
recipient_SEX_stats <- as.data.frame(table(HLI_paired_cohort[c("Group", "R_sex")])) 
recipient_SEX_stats$R_sex <- as.character(recipient_SEX_stats$R_sex) 
recipient_SEX_stats_table <- t(data.frame("Female" = vector(mode = "numeric", length = 2),
                                           "Male" = vector(mode = "numeric", length = 2)))
colnames(recipient_SEX_stats_table) <- c("aGVHD", "nonGVHD")

aGVHD_index <- which(recipient_SEX_stats$Group == "a")
nGVHD_index <- which(recipient_SEX_stats$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  switch(recipient_SEX_stats$R_sex[aGVHD_index[id]],
         "F" = {recipient_SEX_stats_table[1, 1] <- recipient_SEX_stats$Freq[aGVHD_index[id]]},
         "M" = {recipient_SEX_stats_table[2, 1] <- recipient_SEX_stats$Freq[aGVHD_index[id]]}
         )
  
  ##### nonGVHD groups
  switch(recipient_SEX_stats$R_sex[nGVHD_index[id]],
         "F" = {recipient_SEX_stats_table[1, 2] <- recipient_SEX_stats$Freq[nGVHD_index[id]]},
         "M" = {recipient_SEX_stats_table[2, 2] <- recipient_SEX_stats$Freq[nGVHD_index[id]]})
  
}

#### donor SIRE stats
donor_SEX_stats <- as.data.frame(table(HLI_paired_cohort[c("Group", "D_sex")])) 
donor_SEX_stats$D_sex <- as.character(donor_SEX_stats$D_sex) 
donor_SEX_stats_table <- t(data.frame("Female" = vector(mode = "numeric", length = 2),
                                      "Male" = vector(mode = "numeric", length = 2)))
colnames(donor_SEX_stats_table) <- c("aGVHD", "nonGVHD")

aGVHD_index <- which(donor_SEX_stats$Group == "a")
nGVHD_index <- which(donor_SEX_stats$Group == "n")
for(id in 1:length(aGVHD_index)){
  
  ##### aGVHD groups
  switch(donor_SEX_stats$D_sex[aGVHD_index[id]],
         "F" = {donor_SEX_stats_table[1, 1] <- donor_SEX_stats$Freq[aGVHD_index[id]]},
         "M" = {donor_SEX_stats_table[2, 1] <- donor_SEX_stats$Freq[aGVHD_index[id]]}
  )
  
  ##### nonGVHD groups
  switch(donor_SEX_stats$D_sex[nGVHD_index[id]],
         "F" = {donor_SEX_stats_table[1, 2] <- donor_SEX_stats$Freq[nGVHD_index[id]]},
         "M" = {donor_SEX_stats_table[2, 2] <- donor_SEX_stats$Freq[nGVHD_index[id]]})
  
}


######
all_table <- rbind(recipient_age_stats_table, donor_age_stats_table, 
                   recipient_disease_stats_table,  HLI_SIRE_stats_table, 
                   recipient_SIRE_stats_table, donor_SIRE_stats_table, 
                   HLI_Sex_match_table_stats, recipient_SEX_stats_table,
                   donor_SEX_stats_table)

write.csv(all_table, file = "../FirstPaper/HLI_demographics_summary.csv")


######### HLA-Table
# reformated_HLA_typing_list <- read.csv("../ClinVar/reformated_HLA_typing.csv", stringsAsFactors = F)
load("../Data/HLI_reformatted_HLA_table_corrected.RData") ## reformated_HLA_typing_list
num_cases <- dim(reformated_HLA_typing_list)[1]
# 
# reformated_HLA_typing_list$GL_HLA_A1 <- sapply(1:num_cases, function(x) paste0("A*", reformated_HLA_typing_list$HLA_A1[x]))
# reformated_HLA_typing_list$GL_HLA_A2 <- sapply(1:num_cases, function(x) paste0("A*", reformated_HLA_typing_list$HLA_A2[x]))
# reformated_HLA_typing_list$GL_HLA_B1 <- sapply(1:num_cases, function(x) paste0("B*", reformated_HLA_typing_list$HLA_B1[x]))
# reformated_HLA_typing_list$GL_HLA_B2 <- sapply(1:num_cases, function(x) paste0("B*", reformated_HLA_typing_list$HLA_B2[x]))
# reformated_HLA_typing_list$GL_HLA_C1 <- sapply(1:num_cases, function(x) paste0("C*", reformated_HLA_typing_list$HLA_C1[x]))
# reformated_HLA_typing_list$GL_HLA_C2 <- sapply(1:num_cases, function(x) paste0("C*", reformated_HLA_typing_list$HLA_C2[x]))
# reformated_HLA_typing_list$GL_HLA_DRB11 <- sapply(1:num_cases, function(x) paste0("DRB1*", reformated_HLA_typing_list$HLA_DRB11[x]))
# reformated_HLA_typing_list$GL_HLA_DRB12 <- sapply(1:num_cases, function(x) paste0("DRB1*", reformated_HLA_typing_list$HLA_DRB12[x]))
# reformated_HLA_typing_list$GL_HLA_DQB11 <- sapply(1:num_cases, function(x) paste0("DQB1*", reformated_HLA_typing_list$HLA_DQB11[x]))
# reformated_HLA_typing_list$GL_HLA_DQB12 <- sapply(1:num_cases, function(x) paste0("DQB1*", reformated_HLA_typing_list$HLA_DQB12[x]))
# 
# reformated_HLA_typing_list$Group <- sapply(1:num_cases, function(x) unique(paired_ID_table$Group[which(paired_ID_table$caseID %in% reformated_HLA_typing_list$bmt_case_num[x])]))


# restricted_MiHA_table <- read.delim("../WW_MiHA/Restricted_known_MiHAs.txt", header = F)
restricted_MiHA_table <- read.csv("../FirstPaper/Table/HLI_Restricted_MiHA_final_table.csv", stringsAsFactors = F)
# colnames(restricted_MiHA_table) <- c("Group", "GroupID", "HLA", "MiHA", "CHROM", "Donor", "Recipient")
restricted_MiHA_table$HLA <- restricted_MiHA_table$Restricted_HLA

Restricted_HLA <- unique(restricted_MiHA_table$HLA)
Restricted_HLA <- Restricted_HLA[-(14:17)]
num_cases <- dim(reformated_HLA_typing_list)[1]
num_HLA <- length(Restricted_HLA)

HLI_cohort_HLA_case_presAbs <- as.data.frame(matrix(0, nrow = num_cases, ncol = num_HLA))
colnames(HLI_cohort_HLA_case_presAbs) <- Restricted_HLA
HLI_cohort_HLA_case_presAbs$total <- 0
for(id in 1:num_cases){
  
  for(jd in 1:num_HLA){
    
    pres_ind <- which(as.character(reformated_HLA_typing_list[id, 8:17]) %in% Restricted_HLA[jd])
    if(length(pres_ind) > 0){
      HLI_cohort_HLA_case_presAbs[id, jd] <- HLI_cohort_HLA_case_presAbs[id, jd] + 1
    }
    
  }
  if(sum(HLI_cohort_HLA_case_presAbs[id, ]) > 0) HLI_cohort_HLA_case_presAbs$total[id] <- HLI_cohort_HLA_case_presAbs$total[id] + 1
  
  
}

# Restricted_HLA_summary <-colSums(HLI_cohort_HLA_case_presAbs)

aGVHD_id <- which(reformated_HLA_typing_list$GroupType == "a")
nGVHD_id <- which(reformated_HLA_typing_list$GroupType == "n")

Restricted_HLA_summary_aGVHD <- colSums(HLI_cohort_HLA_case_presAbs[aGVHD_id, ])
Restricted_HLA_summary_aGVHD <- as.matrix(Restricted_HLA_summary_aGVHD)
colnames(Restricted_HLA_summary_aGVHD) <- "aGVHD"
Restricted_HLA_summary_nGVHD <- colSums(HLI_cohort_HLA_case_presAbs[nGVHD_id, ])
Restricted_HLA_summary_nGVHD <- as.matrix(Restricted_HLA_summary_nGVHD)
colnames(Restricted_HLA_summary_nGVHD) <- "nGVHD"
Restricted_HLA_summary_all <- cbind(Restricted_HLA_summary_aGVHD, Restricted_HLA_summary_nGVHD)
write.csv(Restricted_HLA_summary_all, file = "../FirstPaper/Table/Restricted_HLA_summary_all_ClassIonly_0630.csv")

Restricted_HLA_CountTable <- as.matrix(Restricted_HLA_summary)
colnames(Restricted_HLA_CountTable) <- "Counts"
write.csv(Restricted_HLA_CountTable, file = "../FirstPaper/Table/Restricted_HLA_counts_table.csv")

length(unique(restricted_MiHA_table$GroupID))
