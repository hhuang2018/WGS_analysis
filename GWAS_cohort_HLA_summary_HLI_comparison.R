source("util.R") 

# convert HLA NMDP MAC into IMGT/HLA nomenclature format
MAC_table <- read.delim(file = "../../Database/alpha.v3.txt", header  = F, stringsAsFactors = F)
MAC_table <- MAC_table[-1, ]
colnames(MAC_table) <- c("FLAG", MAC_table[1, c(2,3)])
MAC_table <- MAC_table[-1, ]

GWAS_HLA_MAC <- read.csv(file="../ClinVar/GWASH/Metadata/GWASH_HLA.csv", stringsAsFactors = F)

GWAS_HLA_IMGT <- GWAS_HLA_MAC
HLA_type_colIND <- 5:48
num_types <- length(HLA_type_colIND)
num_cases <- dim(GWAS_HLA_IMGT)[1]
reformatted_HLA_Type <- matrix(data = NA, nrow = num_cases, ncol = num_types + 1)
reformatted_HLA_Type <- as.data.frame(reformatted_HLA_Type)
colnames(reformatted_HLA_Type) <- c("bmt_caseID", colnames(GWAS_HLA_IMGT)[HLA_type_colIND])

for(id in 1:num_cases){
  
  reformatted_HLA_Type$bmt_caseID[id] <- GWAS_HLA_IMGT$bmt_case_num[id]
  
  for(jd in 1:num_types){
    
    if(!grepl("^[[:digit:]]", unlist(strsplit(GWAS_HLA_IMGT[id, HLA_type_colIND[jd]], ":"))[2])){
      
      temp_code <- unlist(strsplit(GWAS_HLA_IMGT[id, HLA_type_colIND[jd]], ":"))
      
      MAC_ind <- which(MAC_table$CODE %in% temp_code[2])
      if(length(MAC_ind) > 0){
        
        if(MAC_table$FLAG[MAC_ind] == "*"){ ## flagged allele -- whole allele code
          
          decoded_allele <- paste0(unlist(strsplit(temp_code[1], "\\*"))[1], "*", MAC_table$SUBTYPE[MAC_ind])
          two_field_allele <- unlist(strsplit(decoded_allele, "/"))[1]
          
        }else{ # only second field coded
          
          decoded_allele <- paste0(temp_code[1], ":", MAC_table$SUBTYPE[MAC_ind])
          two_field_allele <- unlist(strsplit(decoded_allele, "/"))[1]

        }
        
        GWAS_HLA_IMGT[id, HLA_type_colIND[jd]] <- decoded_allele
        reformatted_HLA_Type[id, jd+1] <- two_field_allele
        
      }else stop("Code does not exsit!")
      # if the code is non-numeric code, then convert
    }else{ # if not the code, then take the first two fields
      
      two_field_allele <- paste0(unlist(strsplit(GWAS_HLA_IMGT[id, HLA_type_colIND[jd]], ":"))[c(1:2)], collapse = ":")
      
      GWAS_HLA_IMGT[id, HLA_type_colIND[jd]] <- two_field_allele
      reformatted_HLA_Type[id, jd+1] <- two_field_allele
      
    }
    
    
  }
  
  
}

GWAS_HLA_IMGT$r_c_typ1 <- gsub("Cw", "C", GWAS_HLA_IMGT$r_c_typ1)
GWAS_HLA_IMGT$r_c_typ2 <- gsub("Cw", "C", GWAS_HLA_IMGT$r_c_typ2)
GWAS_HLA_IMGT$d_c_typ1 <- gsub("Cw", "C", GWAS_HLA_IMGT$d_c_typ1)
GWAS_HLA_IMGT$d_c_typ2 <- gsub("Cw", "C", GWAS_HLA_IMGT$d_c_typ2)

reformatted_HLA_Type$r_c_typ1 <- gsub("Cw", "C", reformatted_HLA_Type$r_c_typ1)
reformatted_HLA_Type$r_c_typ2 <- gsub("Cw", "C", reformatted_HLA_Type$r_c_typ2)
reformatted_HLA_Type$d_c_typ1 <- gsub("Cw", "C", reformatted_HLA_Type$d_c_typ1)
reformatted_HLA_Type$d_c_typ2 <- gsub("Cw", "C", reformatted_HLA_Type$d_c_typ2)

save(GWAS_HLA_MAC, GWAS_HLA_IMGT, reformatted_HLA_Type, file = "../Data/GWAS_cohort_HLA_typing.RData")

############ Check number of HLA cases
# Avail_cases_table <- read.csv(file = "../ClinVar/GWASH/GWASH_available_IDs.csv", stringsAsFactors = F)
Aa <- Avail_cases_table
Avail_cases_table <- Avail_cases_table[-which(Avail_cases_table$bmt_case %in% overlapping_case_IDs), ] # remove overlapping cases
reformatted_HLA_Type$agvhi24 <- sapply(1:dim(reformatted_HLA_Type)[1], function(x) Avail_cases_table$agvhi24[which(Avail_cases_table$bmt_case %in% reformatted_HLA_Type$bmt_caseID[x])])

restricted_MiHA_table <- read.delim("../WW_MiHA/Restricted_known_MiHAs.txt", header = F)
colnames(restricted_MiHA_table) <- c("Group", "GroupID", "HLA", "MiHA", "CHROM", "Donor", "Recipient")
restricted_MiHA_table$HLA <- as.character(restricted_MiHA_table$HLA)

Restricted_HLA <- unique(restricted_MiHA_table$HLA)
num_cases <- dim(reformatted_HLA_Type)[1]
num_HLA <- length(Restricted_HLA)

GWAS_cohort_HLA_case_presAbs <- as.data.frame(matrix(0, nrow = num_cases, ncol = num_HLA))
colnames(GWAS_cohort_HLA_case_presAbs) <- Restricted_HLA
GWAS_cohort_HLA_case_presAbs$total <- 0
for(id in 1:num_cases){
  
  for(jd in 1:num_HLA){
    
    pres_ind <- which(as.character(reformatted_HLA_Type[id, 19:28]) %in% Restricted_HLA[jd])
    if(length(pres_ind) > 0){
      GWAS_cohort_HLA_case_presAbs[id, jd] <- GWAS_cohort_HLA_case_presAbs[id, jd] + 1
    }
    
  }
  if(sum(GWAS_cohort_HLA_case_presAbs[id, ]) > 0) GWAS_cohort_HLA_case_presAbs$total[id] <- GWAS_cohort_HLA_case_presAbs$total[id] + 1
  
  
}

aGVHD_id <- which(reformatted_HLA_Type$agvhi24 == 1)
nGVHD_id <- which(reformatted_HLA_Type$agvhi24 == 0)

Restricted_HLA_summary_aGVHD <- colSums(GWAS_cohort_HLA_case_presAbs[aGVHD_id, ])
Restricted_HLA_summary_aGVHD <- as.matrix(Restricted_HLA_summary_aGVHD)
colnames(Restricted_HLA_summary_aGVHD) <- "aGVHD"
Restricted_HLA_summary_nGVHD <- colSums(GWAS_cohort_HLA_case_presAbs[nGVHD_id, ])
Restricted_HLA_summary_nGVHD <- as.matrix(Restricted_HLA_summary_nGVHD)
colnames(Restricted_HLA_summary_nGVHD) <- "nGVHD"
Restricted_HLA_summary_all <- cbind(Restricted_HLA_summary_aGVHD, Restricted_HLA_summary_nGVHD)
# write.csv(Restricted_HLA_summary_all, file = "../FirstPaper/Table/GWAS_Restricted_HLA_summary_all.csv")
write.csv(Restricted_HLA_summary_all, file = "../FirstPaper/Table/GWAS_noOverlappingSamples_Restricted_HLA_summary_all.csv")

Restricted_HLA_CountTable <- as.matrix(colSums(GWAS_cohort_HLA_case_presAbs))
colnames(Restricted_HLA_CountTable) <- "Counts"
# write.csv(Restricted_HLA_CountTable, file = "../FirstPaper/Table/GWAS_Restricted_HLA_counts_table.csv")
write.csv(Restricted_HLA_CountTable, file = "../FirstPaper/Table/GWAS_noOverlappingSamples_Restricted_HLA_counts_table.csv")

############
# Avail_cases_table <- GWAS_sample_table[which(GWAS_sample_table$pres_abs == 1), ]
# Total: 989 cases (987 pairs with available data)
# AML  - 332 cases 
# ALL  - 123 cases 
# CML  - 343 cases
# MDS  - 191 cases -- could be considered as AML

# aGVHD Grades II-IV  
# Yes = 1 - 146 cases
# No = 0  - 186 cases
#####################
# aGVHD Grades III-IV 
# Yes = 1 - 64 cases
# No = 0  - 268 cases


######################################
# Omni Express Chip SNPs - GRCh37.p13
######################################
# GWASH_snps <- read.delim(file = "../ClinVar/GWASH/labcorp0814.map", header = FALSE)
# colnames(GWASH_snps) <- c("chr", "SNP", "GeneticPos(M)", "start")
# GWASH_snps$end <- GWASH_snps$start
# 
# GWASH_snps$chr <- paste0("chr", GWASH_snps$chr)
# # GWASH_snps$chr[which(GWASH_snps$chr == 0)] <- "chrM"
# 
# Known_MiHA_SNPs <- read.delim("../WW_MiHA/Known_MiHA_coordinates.txt")
# 
# ## liftover 
# # source("http://bioconductor.org/biocLite.R")
# # biocLite("rtracklayer")
# library(rtracklayer)
# # ??liftOver
# 
# chain <- import.chain("../Data/hg19ToHg38.over.chain") ## liftOver chain file
# 
# GWASH_snps_GRanges <- makeGRangesFromDataFrame(GWASH_snps, ignore.strand = TRUE, keep.extra.columns = TRUE)
# 
# GWASH_snps_hg38 <- liftOver(GWASH_snps_GRanges, chain)
# 
# GWASH_snps_hg38_dataFrame <- as.data.frame(GWASH_snps_hg38@unlistData)
# 
# # check overlapped SNPs 
# Known_MiHA_SNPs$chr_POS <- sapply(1:dim(Known_MiHA_SNPs)[1], function(x) paste0(Known_MiHA_SNPs$Chr[x], "-", Known_MiHA_SNPs$Pos[x]))
# GWASH_snps_hg38_dataFrame$chr_POS <- sapply(1:dim(GWASH_snps_hg38_dataFrame)[1], function(x) paste0(GWASH_snps_hg38_dataFrame$seqnames[x], "-", GWASH_snps_hg38_dataFrame$start[x]))
# 
# overlapping_SNPs <- GWASH_snps_hg38_dataFrame[which(GWASH_snps_hg38_dataFrame$chr_POS %in% Known_MiHA_SNPs$chr_POS), ]
# 
# overlapping_SNPs$SNP_currentID <- sapply(1:length(overlapping_SNPs$chr_POS), function(x) Known_MiHA_SNPs$SNPs[which(Known_MiHA_SNPs$chr_POS %in% overlapping_SNPs$chr_POS[x])])

save(overlapping_SNPs, file = "../Data/GWAS_HLI_overlapping_SNPs.RData")

load("../Data/GWAS_HLI_overlapping_SNPs.RData")

#################################################################
## GWAS single MiHA SNPs - HLA restricted
#################################################################

load("../Data/ID_table.RData")

KnownMiHA_table <- read.csv("../ClinVar/Data/KnownMiHA_Table_Ref.csv", stringsAsFactors = F)    

GWASH_HLA_typing_table <- read.csv("../ClinVar/GWASH/Metadata/GWASH_HLA.csv", stringsAsFactors = F)

num_avail_cases <- dim(mismatch_MiHA_table)[1]
SNP_ind <- seq(from=3, to = dim(mismatch_MiHA_table)[2], by = 2)
SNPs <- gsub("_PresAbs", "", colnames(mismatch_MiHA_table)[SNP_ind])

HLA_SNPs <- sapply(1:length(SNPs), function(x) KnownMiHA_table$HLA[which(KnownMiHA_table$rs.number %in% SNPs[x])])
Restricted_MiHA_summary <- data.frame(matrix(0, nrow = 2, ncol = length(SNPs)))
colnames(Restricted_MiHA_summary) <- SNPs
rownames(Restricted_MiHA_summary) <- c("aGVHD", "non-GVHD")
GWASH_MiHA_freq_tab <- data.frame(GroupType = character(num_avail_cases*5),
                                  GroupID = numeric(num_avail_cases*5),
                                  HLA_type = character(num_avail_cases*5),
                                  SNP = character(num_avail_cases*5),
                                  HLA_SNP = character(num_avail_cases*5),
                                  stringsAsFactors = F)
counter <- 0
for(id in 1:num_avail_cases){
  
  pres_SNPs <- which(mismatch_MiHA_table[id, SNP_ind] != 0)
  if(length(pres_SNPs) > 0){
    
    for(jd in 1:length(pres_SNPs)){
      
      restricted_HLA <- KnownMiHA_table$HLA[which(KnownMiHA_table$rs.number %in% SNPs[pres_SNPs[jd]])] 
      
      case_ind <- which(GWASH_HLA_typing_table$bmt_case_num %in% mismatch_MiHA_table$CaseNumber[id])
      
      if(is.element(restricted_HLA, GWASH_HLA_typing_table[case_ind, c(5:14, 21, 22)])){ # Restricted
        if(!is.na(mismatch_MiHA_table$aGVHD[id])){
          counter <- counter + 1
          GWASH_MiHA_freq_tab$GroupID[counter] <- mismatch_MiHA_table$CaseNumber[id]
          
          GWASH_MiHA_freq_tab$HLA_type[counter] <- restricted_HLA
          GWASH_MiHA_freq_tab$SNP[counter] <- SNPs[pres_SNPs[jd]]
          GWASH_MiHA_freq_tab$HLA_SNP[counter] <- paste0(restricted_HLA, "-", SNPs[pres_SNPs[jd]]) 
          
          if(mismatch_MiHA_table$aGVHD[id] == 1){ # aGVHD
            
            Restricted_MiHA_summary[1, SNPs[pres_SNPs[jd]]] <- Restricted_MiHA_summary[1, SNPs[pres_SNPs[jd]]] + 1
            GWASH_MiHA_freq_tab$GroupType[counter] <- "aGVHD"
            
          }else{ # non-aGVHD
            
            Restricted_MiHA_summary[2, SNPs[pres_SNPs[jd]]] <- Restricted_MiHA_summary[2, SNPs[pres_SNPs[jd]]] + 1
            GWASH_MiHA_freq_tab$GroupType[counter] <- "non-aGVHD"
            
          }
          
        }
      }
    }
    
  }
  
}

GWASH_Restricted_MiHA_table <- GWASH_MiHA_freq_tab[which(GWASH_MiHA_freq_tab$GroupID != 0), ]
GWASH_single_MiHA_stats <- as.data.frame(table(GWASH_Restricted_MiHA_table[c("GroupType","HLA_SNP")]))

###### HLI 
# load("../Data/ID_table.RData")
# known_miha_freq <- read.csv("../WW_MiHA/HLA_restricted_knownMiHAs.csv", header = F)
known_miha_freq <- read.delim("../WW_MiHA/Restricted_known_MiHAs.txt", header = F)
colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")
table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq$HLA_SNP <- sapply(1:dim(known_miha_freq)[1], 
                                  function(x) paste0(known_miha_freq$HLA_type[x], "-", known_miha_freq$SNP[x]))

single_MiHA_stats <- as.data.frame(table(known_miha_freq[c("GroupType", "HLA_SNP")]))

GroupID_MiHA_stats <- as.data.frame(table(known_miha_freq[c("GroupID", "HLA_SNP")]))
GroupID_MiHA_stats <- GroupID_MiHA_stats[GroupID_MiHA_stats$Freq>0, ]

unique_IDs <- unique(known_miha_freq$GroupID) # 131 groups - restricted; 205 groups - unrestricted

SNP_mat <- matrix(data = 0, nrow = length(unique_IDs), ncol = length(unique(GroupID_MiHA_stats$HLA_SNP)) + 2)
colnames(SNP_mat) <- c("GroupID", "GroupType", as.character(unique(GroupID_MiHA_stats$HLA_SNP)))
SNP_mat <- as.data.frame(SNP_mat)

for(id in 1:length(unique_IDs)){
  
  selected_group <- known_miha_freq[known_miha_freq$GroupID == unique_IDs[id],]
  SNP_mat$GroupID[id] <- unique_IDs[id]
  SNP_mat$GroupType[id] <- ID_table$Group[which(ID_table$GroupID == unique_IDs[id])[1]]
  
  if(length(unique(selected_group$HLA_SNP)) == length(selected_group$HLA_SNP)){ # no duplicated SNPs for one group
    
    SNP_mat[id, selected_group$HLA_SNP] <- 1
    
  }else{ # one homozygous SNP (duplicated)
    
    SNP_counts <- as.data.frame(table(selected_group$HLA_SNP), stringsAsFactors = F)
    SNP_mat[id, SNP_counts[, 1]] <- SNP_counts[, 2]
    
  }
  
}


#######


###### save HLA restricted MiHA table
# write.csv(GWASH_single_MiHA_table, file = "../ClinVar/GWASH/GWASH_cohort_Restricted_MiHA_summary.csv", row.names = F)
#####
# write.csv(HLI_single_MiHA_table, file = "../ClinVar/GWASH/HLI_cohort_Overlapping_Restricted_MiHA_summary.csv", row.names = F)

# ######### HLI GWASH overlapping cases
# HLA_typings <- read.csv(file = "../HLI_hla_mg_v3.csv")
# 
# load("../Data/ID_table.RData")
# 
# available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]
# 
# avaialble_IDs_HLA_typing <- HLA_typings[(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID),]
# num_groups <- dim(avaialble_IDs_HLA_typing)[1]
# 
# overlapping_case_IDs <- intersect(Avail_cases_table$bmt_case, avaialble_IDs_HLA_typing$bmt_case_num)
## 57 cases
save(overlapping_case_IDs, file = "../Data/GWAS_HLI_Overlapping_Cases_availablePairs.RData")
# overlapping_case_IDs <- intersect(GWAS_sample_table$bmt_case, avaialble_IDs_HLA_typing$bmt_case_num)
# ## 69 cases
# save(overlapping_case_IDs, file = "../Data/GWAS_HLI_Overlapping_Cases_All.RData")

###### GWASH table
KnownMiHA_table <- read.csv("../ClinVar/Data/KnownMiHA_Table_Ref.csv", stringsAsFactors = F)
KnownMiHA_table$HLA_SNP <- sapply(1:dim(KnownMiHA_table)[1], function(x) paste0(KnownMiHA_table$HLA[x], "-", KnownMiHA_table$rs.number[x]))

GWASH_HLA_typing_table <- read.csv("../ClinVar/GWASH/Metadata/GWASH_HLA.csv", stringsAsFactors = F)

# overlappling_mimatch_table <- mismatch_MiHA_table[which(mismatch_MiHA_table$CaseNumber %in% overlapping_case_IDs), ]
overlappling_mimatch_table <- mismatch_MiHA_table[-which(mismatch_MiHA_table$CaseNumber %in% overlapping_case_IDs), ] # non-overlapping
num_avail_cases <- dim(overlappling_mimatch_table)[1]
SNP_ind <- seq(from=3, to = dim(overlappling_mimatch_table)[2], by = 2)
SNPs <- gsub("_PresAbs", "", colnames(overlappling_mimatch_table)[SNP_ind])

HLA_SNPs <- sapply(1:length(SNPs), function(x) KnownMiHA_table$HLA[which(KnownMiHA_table$rs.number %in% SNPs[x])])
Restricted_MiHA_summary <- data.frame(matrix(0, nrow = 2, ncol = length(SNPs)))
colnames(Restricted_MiHA_summary) <- SNPs
rownames(Restricted_MiHA_summary) <- c("aGVHD", "non-GVHD")
overlapping_GWASH_MiHA_freq_tab <- data.frame(GroupType = character(num_avail_cases*5),
                                              GroupID = numeric(num_avail_cases*5),
                                              HLA_type = character(num_avail_cases*5),
                                              SNP = character(num_avail_cases*5),
                                              HLA_SNP = character(num_avail_cases*5),
                                              stringsAsFactors = F)
counter <- 0
for(id in 1:num_avail_cases){
  
  pres_SNPs <- which(overlappling_mimatch_table[id, SNP_ind] != 0)
  if(length(pres_SNPs) > 0){
    
    for(jd in 1:length(pres_SNPs)){
      
      restricted_HLA <- KnownMiHA_table$HLA[which(KnownMiHA_table$rs.number %in% SNPs[pres_SNPs[jd]])] 
      
      case_ind <- which(GWASH_HLA_typing_table$bmt_case_num %in% overlappling_mimatch_table$CaseNumber[id])
      
      if(is.element(restricted_HLA, GWASH_HLA_typing_table[case_ind, c(5:14, 21, 22)])){ # Restricted
        if(!is.na(overlappling_mimatch_table$aGVHD[id])){
          counter <- counter + 1
          overlapping_GWASH_MiHA_freq_tab$GroupID[counter] <- overlappling_mimatch_table$CaseNumber[id]
          
          overlapping_GWASH_MiHA_freq_tab$HLA_type[counter] <- restricted_HLA
          overlapping_GWASH_MiHA_freq_tab$SNP[counter] <- SNPs[pres_SNPs[jd]]
          overlapping_GWASH_MiHA_freq_tab$HLA_SNP[counter] <- paste0(restricted_HLA, "-", SNPs[pres_SNPs[jd]]) 
          
          if(overlappling_mimatch_table$aGVHD[id] == 1){ # aGVHD
            
            Restricted_MiHA_summary[1, SNPs[pres_SNPs[jd]]] <- Restricted_MiHA_summary[1, SNPs[pres_SNPs[jd]]] + 1
            overlapping_GWASH_MiHA_freq_tab$GroupType[counter] <- "aGVHD"
            
          }else{ # non-aGVHD
            
            Restricted_MiHA_summary[2, SNPs[pres_SNPs[jd]]] <- Restricted_MiHA_summary[2, SNPs[pres_SNPs[jd]]] + 1
            overlapping_GWASH_MiHA_freq_tab$GroupType[counter] <- "non-aGVHD"
            
          }
          
        }
      }
    }
    
  }
  
}

overlapping_GWASH_Restricted_MiHA_table <- overlapping_GWASH_MiHA_freq_tab[which(overlapping_GWASH_MiHA_freq_tab$GroupID != 0), ]
overlapping_GWASH_single_MiHA_stats <- as.data.frame(table(overlapping_GWASH_Restricted_MiHA_table[c("GroupType","HLA_SNP")]))


###### save HLA restricted MiHA table
num_MiHAs <- dim(GWASH_single_MiHA_stats)[1]/2
overlap_GWASH_single_MiHA_table <- data.frame(Gene = character(num_MiHAs), 
                                              MiHA = character(num_MiHAs), 
                                              rsNumber = character(num_MiHAs), 
                                              HLA = character(num_MiHAs),
                                              aGVHD = numeric(num_MiHAs), 
                                              nGVHD = numeric(num_MiHAs), 
                                              LLR10_avg = numeric(num_MiHAs),
                                              stringsAsFactors = F)

for(miha_id in 1:num_MiHAs){
  
  index <- which(overlapping_GWASH_single_MiHA_stats$HLA_SNP %in% unique_MiHAs[miha_id])
  if(length(index) >0){
    
    overlap_GWASH_single_MiHA_table$aGVHD[miha_id] <- overlapping_GWASH_single_MiHA_stats[index[which(overlapping_GWASH_single_MiHA_stats[index, "GroupType"] == "aGVHD")], "Freq"]
    overlap_GWASH_single_MiHA_table$nGVHD[miha_id] <- abs(overlapping_GWASH_single_MiHA_stats[index[which(overlapping_GWASH_single_MiHA_stats[index, "GroupType"] == "non-aGVHD")], "Freq"])
    
    iid <- which(KnownMiHA_table$HLA_SNP %in% unique(overlapping_GWASH_single_MiHA_stats$HLA_SNP[index]))
    # Gene_MiHA_name <- paste0(as.character(KnownMiHA_table[iid, c("Gene", "KnownMiHA", "rs.number")]), collapse = " <> ")
    overlap_GWASH_single_MiHA_table$Gene[miha_id] <- paste0(unique(KnownMiHA_table[iid, "Gene"]), collapse = "/")
    overlap_GWASH_single_MiHA_table$MiHA[miha_id] <- paste0(unique(KnownMiHA_table[iid, "KnownMiHA"]), collapse = "/")
    overlap_GWASH_single_MiHA_table$rsNumber[miha_id] <- paste0(unique(KnownMiHA_table[iid, "rs.number"]), collapse = "/")
    overlap_GWASH_single_MiHA_table$HLA[miha_id] <- paste0(unique(KnownMiHA_table[iid, "HLA"]), collapse = "/")
    
    if(aGVHD_count > 0 | nGVHD_count >0){
      
      overlap_GWASH_single_MiHA_table$LLR10_avg[miha_id] <- log10(overlap_GWASH_single_MiHA_table$aGVHD[miha_id]/overlap_GWASH_single_MiHA_table$nGVHD[miha_id])
      # GWASH_single_MiHA_table$name[miha_id] <- Gene_MiHA_name
    }
    
  }
}

overlap_GWASH_single_MiHA_table <- overlap_GWASH_single_MiHA_table[which(overlap_GWASH_single_MiHA_table$LLR10_avg != 0 ), ]

write.csv(overlap_GWASH_single_MiHA_table, file = "../ClinVar/GWASH/GWASH_overlapping_cohort_Overlapping_Restricted_MiHA_summary.csv", 
          row.names = F)

write.csv(overlap_GWASH_single_MiHA_table, file = "../ClinVar/GWASH/GWASH_Non-overlapping_cohort_Overlapping_Restricted_MiHA_summary.csv", 
          row.names = F)

overlappling_mimatch_table$GroupID <- unlist(sapply(1:length(overlappling_mimatch_table$CaseNumber), function(x) if (length(which(ID_table$CaseNum %in%  overlappling_mimatch_table$CaseNumber[x])) ==1 ) unique(ID_table$GroupID[which(ID_table$CaseNum %in%  overlappling_mimatch_table$CaseNumber[x])])  
                                                    else "0"))
write.csv(overlappling_mimatch_table, file = "../ClinVar/GWASH/GWASH_overlapping_cohort_count_table.csv", row.names = F)


write.csv(overlappling_mimatch_table, file = "../ClinVar/GWASH/GWASH_Non-overlapping_cohort_count_table.csv", row.names = F)

###### HLI 
load("../Data/GWAS_HLI_overlapping_SNPs.RData")
overlapping_MiHA_ref <- KnownMiHA_table[which(KnownMiHA_table$rs.number %in% overlapping_SNPs$SNP_currentID), ]
unique_MiHAs <- unique(overlapping_MiHA_ref$HLA_SNP)

known_miha_freq <- read.delim("../WW_MiHA/Restricted_known_MiHAs.txt", header = F)
colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")
table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq$HLA_SNP <- sapply(1:dim(known_miha_freq)[1], 
                                  function(x) paste0(known_miha_freq$HLA_type[x], "-", known_miha_freq$SNP[x]))
# overlap_known_miha_freq <- known_miha_freq[which(known_miha_freq$GroupID %in% over_HLI_groupID), ]

overlap_known_miha_freq <- known_miha_freq # all HLI cohort
overlap_single_MiHA_stats <- as.data.frame(table(overlap_known_miha_freq[c("GroupType", "HLA_SNP")]))

index <- which(overlap_single_MiHA_stats$HLA_SNP %in% unique_MiHAs)

overlap_HLI_single_MiHA_stats <- overlap_single_MiHA_stats[index, ]

num_MiHAs <- length(unique_MiHAs)
overlappling_HLI_single_MiHA_table <- data.frame(Gene = character(num_MiHAs), 
                                                 MiHA = character(num_MiHAs), 
                                                 rsNumber = character(num_MiHAs), 
                                                 HLA = character(num_MiHAs),
                                                 aGVHD = numeric(num_MiHAs), 
                                                 nGVHD = numeric(num_MiHAs), 
                                                 LLR10_avg = numeric(num_MiHAs),
                                                 stringsAsFactors = F)

for(miha_id in 1:num_MiHAs){
  
  index <- which(overlap_HLI_single_MiHA_stats$HLA_SNP %in% unique_MiHAs[miha_id])
  if(length(index) > 0){
    overlappling_HLI_single_MiHA_table$aGVHD[miha_id] <- overlap_HLI_single_MiHA_stats[index[which(overlap_HLI_single_MiHA_stats[index, "GroupType"] == "a")], "Freq"]
    overlappling_HLI_single_MiHA_table$nGVHD[miha_id] <- abs(overlap_HLI_single_MiHA_stats[index[which(overlap_HLI_single_MiHA_stats[index, "GroupType"] == "n")], "Freq"])
    
    iid <- which(KnownMiHA_table$HLA_SNP %in% overlap_HLI_single_MiHA_stats$HLA_SNP[index])
    # Gene_MiHA_name <- paste0(as.character(KnownMiHA_table[iid, c("Gene", "KnownMiHA", "rs.number")]), collapse = " <> ")
    overlappling_HLI_single_MiHA_table$Gene[miha_id] <- paste0(unique(KnownMiHA_table[iid, "Gene"]), collapse = "/")
    overlappling_HLI_single_MiHA_table$MiHA[miha_id] <- paste0(unique(KnownMiHA_table[iid, "KnownMiHA"]), collapse = "/")
    overlappling_HLI_single_MiHA_table$rsNumber[miha_id] <- paste0(unique(KnownMiHA_table[iid, "rs.number"]), collapse = "/")
    overlappling_HLI_single_MiHA_table$HLA[miha_id] <- paste0(unique(KnownMiHA_table[iid, "HLA"]), collapse = "/")
    
    if(overlappling_HLI_single_MiHA_table$aGVHD[miha_id] > 0 | overlappling_HLI_single_MiHA_table$nGVHD[miha_id] >0){
      
      overlappling_HLI_single_MiHA_table$LLR10_avg[miha_id] <- log10(overlappling_HLI_single_MiHA_table$aGVHD[miha_id]/overlappling_HLI_single_MiHA_table$nGVHD[miha_id])
      # GWASH_single_MiHA_table$name[miha_id] <- Gene_MiHA_name
    }
  }
  
}

overlappling_HLI_single_MiHA_table <- overlappling_HLI_single_MiHA_table[which(overlappling_HLI_single_MiHA_table$LLR10_avg != 0), ]

write.csv(overlappling_HLI_single_MiHA_table, file = "../ClinVar/GWASH/HLI_overlapping_cohort_Overlapping_Restricted_MiHA_summary.csv", 
          row.names = F)

write.csv(overlappling_HLI_single_MiHA_table, file = "../ClinVar/GWASH/HLI_Non-overlapping_cohort_Overlapping_Restricted_MiHA_summary.csv", 
          row.names = F)

overlap_SNP_mat <- SNP_mat[which(SNP_mat$GroupID %in% over_HLI_groupID), ]
write.csv(overlap_SNP_mat, file = "../ClinVar/GWASH/HLI_overlapping_cohort_count_table.csv", row.names = F)

