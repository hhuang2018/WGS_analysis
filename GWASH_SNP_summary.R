######################################
# Omni Express Chip SNPs - GRCh37.p13
######################################
GWASH_snps <- read.delim(file = "../ClinVar/GWASH/labcorp0814.map", header = FALSE)
colnames(GWASH_snps) <- c("chr", "SNP", "GeneticPos(M)", "start")
GWASH_snps$end <- GWASH_snps$start

GWASH_snps$chr <- paste0("chr", GWASH_snps$chr)
# GWASH_snps$chr[which(GWASH_snps$chr == 0)] <- "chrM"


Known_MiHA_SNPs <- read.delim("../WW_MiHA/Known_MiHA_coordinates.txt")

## liftover 
# source("http://bioconductor.org/biocLite.R")
# biocLite("rtracklayer")
library(rtracklayer)
# ??liftOver

chain <- import.chain("../Data/hg19ToHg38.over.chain") ## liftOver chain file

GWASH_snps_GRanges <- makeGRangesFromDataFrame(GWASH_snps, ignore.strand = TRUE, keep.extra.columns = TRUE)

GWASH_snps_hg38 <- liftOver(GWASH_snps_GRanges, chain)

GWASH_snps_hg38_dataFrame <- as.data.frame(GWASH_snps_hg38@unlistData)

# check overlapped SNPs 
Known_MiHA_SNPs$chr_POS <- sapply(1:dim(Known_MiHA_SNPs)[1], function(x) paste0(Known_MiHA_SNPs$Chr[x], "-", Known_MiHA_SNPs$Pos[x]))
GWASH_snps_hg38_dataFrame$chr_POS <- sapply(1:dim(GWASH_snps_hg38_dataFrame)[1], function(x) paste0(GWASH_snps_hg38_dataFrame$seqnames[x], "-", GWASH_snps_hg38_dataFrame$start[x]))

overlapping_SNPs <- GWASH_snps_hg38_dataFrame[which(GWASH_snps_hg38_dataFrame$chr_POS %in% Known_MiHA_SNPs$chr_POS), ]

overlapping_SNPs$SNP_currentID <- sapply(1:length(overlapping_SNPs$chr_POS), function(x) Known_MiHA_SNPs$SNPs[which(Known_MiHA_SNPs$chr_POS %in% overlapping_SNPs$chr_POS[x])])

##############
# Double-check the overlappin cases with overlapping MiHA SNPs
#############
source("util.R")

GWAS_samp_table_fp <- "../ClinVar/GWASH/Metadata/newfnldataib1203.csv"
GWAS_file_fp <- "../ClinVar/GWASH/MIHA/data/"

GWAS_sample_table <-read.csv(GWAS_samp_table_fp, header = T, stringsAsFactors = F)

MiHA_SNPs <- read.table(paste0(GWAS_file_fp, "MiHAgeno.map"), stringsAsFactors = F) # overlapping MiHA 
colnames(MiHA_SNPs) <- c("Map_ID", "SNP_ID", "Location", "POS")
num_SNPs <- length(MiHA_SNPs$SNP_ID)

MiHA_genotypes <- read.table(paste0(GWAS_file_fp, "MiHAgeno.ped"), stringsAsFactors = F)
colnames(MiHA_genotypes) <- c("MapID", "SubjectID", "Ped1", "Ped2", "Ped3", "Ped4", 
                              as.vector(sapply(1:num_SNPs, function(x) paste0(MiHA_SNPs$SNP_ID[x], c("_Maj", "_Min")))))
SNP_index <- 1:(2*num_SNPs) + 6


num_samples <- dim(MiHA_genotypes)[1]

MiHA_genotypes$R_D_ID <- sapply(1:length(MiHA_genotypes$SubjectID), function(x) as.numeric(gsub("-", "", MiHA_genotypes$SubjectID[x])))

num_cases <- dim(GWAS_sample_table)[1]            # 1378 cases with metadata
GWAS_sample_table$pres_abs <- numeric(num_cases)
MiHA_genotypes$R_D <- character(num_samples)
MiHA_genotypes$case_id <- numeric(num_samples)

for(id in 1:num_samples){
  
  if(MiHA_genotypes$R_D[id] == ""){
    
    num_char <- nchar(MiHA_genotypes$SubjectID[id]) # 11 - donor; 9 - recipient 
    if(num_char == 11){ # donor
      
      MiHA_genotypes$R_D[id] <- "donor"
      donor_ind <- which(GWAS_sample_table$did %in% MiHA_genotypes$R_D_ID[id])
      
      if(length(donor_ind) > 1){ # multiple donor for different cases
        
        r1_ind <- which(MiHA_genotypes$R_D_ID %in% GWAS_sample_table$rid[donor_ind])
        if(length(r1_ind) == 1){
          # GWAS_sample_table$pres_abs[donor_ind] <- 1
          MiHA_genotypes$case_id[id] <- GWAS_sample_table$bmt_case[which(GWAS_sample_table$rid == MiHA_genotypes$R_D_ID[r1_ind])]
          
          # r_index <- which(MiHA_genotypes$R_D_ID %in% GWAS_sample_table$rid[donor_ind])
          field3 <- GWAS_sample_table$rid[which(GWAS_sample_table$rid == MiHA_genotypes$R_D_ID[r1_ind])] %% 10
          field2 <- floor(GWAS_sample_table$rid[which(GWAS_sample_table$rid == MiHA_genotypes$R_D_ID[r1_ind])]/10) %% 1000
          if(nchar(field2) < 3) field2 <- paste0(paste0(rep("0", 3-nchar(field2)), collapse = ""), field2)
          
          field1 <- floor(GWAS_sample_table$rid[which(GWAS_sample_table$rid == MiHA_genotypes$R_D_ID[r1_ind])]/10000) 
          if(nchar(field1) < 3) field1 <- paste0(paste0(rep("0", 3-nchar(field1)), collapse = ""), field1)
        }else cat(id)
        
      }else { # only one Donor 
        MiHA_genotypes$case_id[id] <- GWAS_sample_table$bmt_case[donor_ind]
        # r_index <- which(MiHA_genotypes$R_D_ID %in% GWAS_sample_table$rid[donor_ind])
        field3 <- GWAS_sample_table$rid[donor_ind] %% 10
        field2 <- floor(GWAS_sample_table$rid[donor_ind]/10) %% 1000
        if(nchar(field2) < 3) field2 <- paste0(paste0(rep("0", 3-nchar(field2)), collapse = ""), field2)
        
        field1 <- floor(GWAS_sample_table$rid[donor_ind]/10000) 
        if(nchar(field1) < 3) field1 <- paste0(paste0(rep("0", 3-nchar(field1)), collapse = ""), field1)
      }
      rid <- paste0(field1, "-", field2, "-", field3)
      # if(nchar(rid) < 9) rid <- paste0(paste0(rep("0", 9-nchar(rid)), collapse = ""), rid)
      
      r_index <- which(MiHA_genotypes$SubjectID %in% rid)
      if(length(r_index) > 0){
        GWAS_sample_table$pres_abs[donor_ind] <- 1
        MiHA_genotypes$R_D[r_index] <- "recipient"
        if(MiHA_genotypes$case_id[r_index] == 0){
          
          MiHA_genotypes$case_id[r_index] <- MiHA_genotypes$case_id[id]
          
        }
        
      }
      
    }else if(num_char == 9){# recipient
      
      MiHA_genotypes$R_D[id] <- "recipient"
      recipient_ind <- which(GWAS_sample_table$rid %in% MiHA_genotypes$R_D_ID[id])
      
      MiHA_genotypes$case_id[id] <- GWAS_sample_table$bmt_case[recipient_ind]
      GWAS_sample_table$pres_abs[donor_ind] <- 1
      # d_index <- which(MiHA_genotypes$R_D_ID %in% GWAS_sample_table$rid[recipient_ind])
      field3 <- GWAS_sample_table$rid[recipient_ind] %% 10
      field2 <- floor(GWAS_sample_table$rid[recipient_ind]/10) %% 10000
      if(nchar(field2) < 4) field2 <- paste0(paste0(rep("0", 4-nchar(field2)), collapse = ""), field2)
      
      field1 <- floor(GWAS_sample_table$rid[recipient_ind]/100000) 
      if(nchar(field1) < 4) field2 <- paste0(paste0(rep("0", 4-nchar(field1)), collapse = ""), field1)
      
      did <- paste0(field1, "-", field2, "-", field3)
      # if(nchar(did) < 11) did <- paste0(rep("0", 11-nchar(did)), did)
      
      d_index <- which(MiHA_genotypes$SubjectID %in% did)
      if(length(d_index) > 0){
        GWAS_sample_table$pres_abs[recipient_ind] <- 1
        MiHA_genotypes$R_D[d_index] <- "donor"
        
        MiHA_genotypes$case_id[d_index] <- GWAS_sample_table$bmt_case[recipient_ind]
        
      }
      
      
    }
    # donor_ind <- which(GWAS_sample_table$did %in% MiHA_genotypes$R_D_ID[id])
    # recipient_ind <- which(GWAS_sample_table$rid %in% MiHA_genotypes$R_D_ID[id])
    
    # if(length(donor_ind) > 1) cat(id)
    # if(length(recipient_ind) > 1) cat(id)
    # 
    # if(length(donor_ind) > 0){
    #   MiHA_genotypes$R_D[id] <- "donor"
    #   
    #   if(length(donor_ind) > 1){
    #     
    #     r1_ind <- which(MiHA_genotypes$R_D_ID %in% GWAS_sample_table$rid[donor_ind])
    #     if(length(r1_ind) == 1){
    #       
    #       MiHA_genotypes$case_id[id] <- GWAS_sample_table$bmt_case[which(GWAS_sample_table$rid == MiHA_genotypes$R_D_ID[r1_ind])]
    #       
    #       # r_index <- which(MiHA_genotypes$R_D_ID %in% GWAS_sample_table$rid[donor_ind])
    #       field3 <- GWAS_sample_table$rid[which(GWAS_sample_table$rid == MiHA_genotypes$R_D_ID[r1_ind])] %% 10
    #       field2 <- floor(GWAS_sample_table$rid[which(GWAS_sample_table$rid == MiHA_genotypes$R_D_ID[r1_ind])]/10) %% 1000
    #       field1 <- floor(GWAS_sample_table$rid[which(GWAS_sample_table$rid == MiHA_genotypes$R_D_ID[r1_ind])]/10000) 
    #     }else cat(id)
    #     
    #   }else { 
    #     MiHA_genotypes$case_id[id] <- GWAS_sample_table$bmt_case[donor_ind]
    #     # r_index <- which(MiHA_genotypes$R_D_ID %in% GWAS_sample_table$rid[donor_ind])
    #     field3 <- GWAS_sample_table$rid[donor_ind] %% 10
    #     field2 <- floor(GWAS_sample_table$rid[donor_ind]/10) %% 1000
    #     field1 <- floor(GWAS_sample_table$rid[donor_ind]/10000) 
    #   }
    #   
    #   
    #   r_index <- which(MiHA_genotypes$SubjectID %in% paste0(field1, "-", field2, "-", field3))
    #   if(length(r_index) > 0){
    #     GWAS_sample_table$pres_abs[donor_ind] <- 1
    #     MiHA_genotypes$R_D[r_index] <- "recipient"
    #     if(MiHA_genotypes$case_id[r_index] == 0){
    #     
    #       MiHA_genotypes$case_id[r_index] <- MiHA_genotypes$case_id[id]
    #       
    #     }
    #     
    #   }
    # }
    # 
    #   if(length(recipient_ind) > 0){
    #     MiHA_genotypes$R_D[id] <- "recipient"
    #     
    #     MiHA_genotypes$case_id[id] <- GWAS_sample_table$bmt_case[recipient_ind]
    #     
    #     # d_index <- which(MiHA_genotypes$R_D_ID %in% GWAS_sample_table$rid[recipient_ind])
    #     field3 <- GWAS_sample_table$rid[recipient_ind] %% 10
    #     field2 <- floor(GWAS_sample_table$rid[recipient_ind]/10) %% 1000
    #     field1 <- floor(GWAS_sample_table$rid[recipient_ind]/10000) 
    #     
    #     d_index <- which(MiHA_genotypes$SubjectID %in% paste0(field1, "-", field2, "-", field3))
    #     if(length(d_index) > 0){
    #       GWAS_sample_table$pres_abs[recipient_ind] <- 1
    #       MiHA_genotypes$R_D[d_index] <- "donor"
    #       
    #       MiHA_genotypes$case_id[d_index] <- GWAS_sample_table$bmt_case[recipient_ind]
    #       
    #     }
    #   }
  }
  
} # add metadata to MIHA genotype table

##### Missing RID
#     BMT_ID Freq
# 264 10125    1
# 476 22879    1
# 583 39449    1
# 597 40223    1
# 874 47726    1
# 903 48546    1
# 961 52067    1

### ADD CASE index to the MiHA table
# MiHA_genotypes$case_id <- numeric(num_samples)
# for(id in 1:num_samples){
#   
#   if(MiHA_genotypes$R_D[id] == "donor"){
#     
#     index <- which(Avail_cases_table$did == MiHA_genotypes$R_D_ID[id])
#     MiHA_genotypes$case_id[id] <- Avail_cases_table$bmt_case[index]
#     
#   }else{
#     
#     index <- which(Avail_cases_table$rid == MiHA_genotypes$R_D_ID[id])
#     MiHA_genotypes$case_id[id] <- Avail_cases_table$bmt_case[index]
#     
#   }
#   
# }

## disease:
# 10 - AML;
# 20 - ALL;
# 40 - CML;
# 50 - MDS;

Avail_cases_table <- GWAS_sample_table[which(GWAS_sample_table$pres_abs == 1), ]
# Total: 988 cases (987 pairs with available data)
# AML  - 332 cases 
# ALL  - 123 cases 
# CML  - 343 cases
# MDS  - 191 cases -- could be considered as AML

write.csv(Avail_cases_table[, 1:3], file = "../ClinVar/GWASH/GWASH_available_IDs.csv", row.names = F)

AML_cases_table <- Avail_cases_table[which(Avail_cases_table$disease == 10), ] # 332
# aGVHD Grades II-IV  
# Yes = 1 - 146 cases
# No = 0  - 186 cases
#####################
# aGVHD Grades III-IV 
# Yes = 1 - 64 cases
# No = 0  - 268 cases

############################
## MiHA Combination
AML_cases_table <- Avail_cases_table

Num_AML_cases <- dim(AML_cases_table)[1]
mismatch_MiHA_table <- matrix(data = 0, nrow = Num_AML_cases, ncol = 2*num_SNPs+2)
colnames(mismatch_MiHA_table) <- c("CaseNumber", "aGVHD", 
                                   as.vector(sapply(1:num_SNPs, function(x) paste0(MiHA_SNPs$SNP_ID[x], c("_PresAbs", "_Number"))))) 
## presense-absence of the MihA mimatch + number of mismatched alleles


mismatch_MiHA_table <- as.data.frame(mismatch_MiHA_table)
counter <- 0
for(id in 1:Num_AML_cases){
  
  index <- which(MiHA_genotypes$case_id %in% AML_cases_table$bmt_case[id])
  
  if(length(index) == 2){ # matching cases
    counter <- counter + 1
    mismatch_MiHA_table$CaseNumber[counter] <- AML_cases_table$bmt_case[id]
    mismatch_MiHA_table$aGVHD[counter] <- AML_cases_table$agvhi24[id]
    
    for(jd in 1:num_SNPs){
      
      mismatch_MiHA_table[counter, 4+2*(jd-1)] <- identical(MiHA_genotypes[index[1], SNP_index[1+2*(jd-1)]], MiHA_genotypes[index[2], SNP_index[1+2*(jd-1)]]) + 
        identical(MiHA_genotypes[index[1], SNP_index[2*jd]], MiHA_genotypes[index[2], SNP_index[2*jd]])
      
      mismatch_MiHA_table[counter, 3+2*(jd-1)] <- as.integer(identical(MiHA_genotypes[index[1], SNP_index[1+2*(jd-1)]], MiHA_genotypes[index[2], SNP_index[1+2*(jd-1)]]) | 
                                                               identical(MiHA_genotypes[index[1], SNP_index[2*jd]], MiHA_genotypes[index[2], SNP_index[2*jd]]))
      
    }
    
    
  }
  
  
}

mismatch_MiHA_table <- mismatch_MiHA_table[1:counter, ]

all_possible_pairs_num <- num_SNPs * (num_SNPs -1) / 2
# all_possible_pairs_names <- character(all_possible_pairs_num)
a_pair_wise_combination <- as.data.frame(matrix(data = 0, nrow = all_possible_pairs_num, ncol = 5), stringsAsFactors = F) # col: abs - abs 0-0;  pres-abs 1-0;  abs-pres 0-1; pres - pres 1-1
colnames(a_pair_wise_combination) <- c("SNPComb", "abs-abs", "pres-abs", "abs-pres", "pres-pres")
n_pair_wise_combination <- as.data.frame(matrix(data = 0, nrow = all_possible_pairs_num, ncol = 5), stringsAsFactors = F)
colnames(n_pair_wise_combination) <- c("SNPComb", "abs-abs", "pres-abs", "abs-pres", "pres-pres")

aGVHD_index <- which(mismatch_MiHA_table$aGVHD == 1) # 146
nGVHD_index <- which(mismatch_MiHA_table$aGVHD == 0) # 186
counter <- 0 
for(id in 1:(num_SNPs-1)){
  
  for(jd in (id+1):num_SNPs){
    # if(id != jd){
    
    counter <- counter + 1
    
    sum_tab <- as.data.frame(table(mismatch_MiHA_table[aGVHD_index, c(3 + 2*(id-1), 3 + 2*(jd-1))]), stringsAsFactors = F)
    
    SNP_1 <- gsub("_PresAbs", "", colnames(sum_tab)[1])
    SNP_2 <- gsub("_PresAbs", "", colnames(sum_tab)[2])
    
    a_pair_wise_combination[counter, 1] <- paste0(SNP_1, "-", SNP_2)
    a_pair_wise_combination[counter, 2:5] <- as.numeric(sum_tab$Freq)
    
    if(!identical(as.numeric(sum_tab[2, c(1,2)]),  c(1, 0))) cat("aGVHD: id = ", id, "; jd = ", jd, "\n")
    
    nsum_tab <- as.data.frame(table(mismatch_MiHA_table[nGVHD_index, c(3 + 2*(id-1), 3 + 2*(jd-1))]), stringsAsFactors = F)
    
    nSNP_1 <- gsub("_PresAbs", "", colnames(nsum_tab)[1])
    nSNP_2 <- gsub("_PresAbs", "", colnames(nsum_tab)[2])
    
    n_pair_wise_combination[counter, 1] <- paste0(nSNP_1, "-", nSNP_2)
    n_pair_wise_combination[counter, 2:5] <- as.numeric(nsum_tab$Freq)
    
    if(!identical(as.numeric(nsum_tab[2, c(1,2)]),  c(1, 0))) cat("nGVHD: id = ", id, "; jd = ", jd, "\n")
    
    
  }
  
  # }
  
}


##### 45 MiHA SNPs --- Outcome vs rs-rs association test
a_joint_SNPs <- read.csv(file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_summary_AP.csv", stringsAsFactors = F)
n_joint_SNPs <- read.csv(file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_summary_AP.csv", stringsAsFactors = F)

a_joint_snps_ordered <- read.csv(file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_ordered_AP.csv", stringsAsFactors = F)
n_joint_snps_ordered <- read.csv(file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_ordered_AP.csv", stringsAsFactors = F)

a_overlapped_combination <- intersect(a_joint_SNPs$X, a_pair_wise_combination$SNPComb) # 102 
n_overlapped_combination <- intersect(n_joint_SNPs$X, n_pair_wise_combination$SNPComb) # 102

a_ordered_name <- sapply(1:length(a_joint_SNPs$X), function(x) paste0(sort(unlist(strsplit(a_joint_SNPs$X[x], "-"))), collapse = "-"))
n_ordered_name <- sapply(1:length(n_joint_SNPs$X), function(x) paste0(sort(unlist(strsplit(n_joint_SNPs$X[x], "-"))), collapse = "-"))

aaa_ordered_name <- sapply(1:length(a_pair_wise_combination$SNPComb), function(x) paste0(sort(unlist(strsplit(a_pair_wise_combination$SNPComb[x], "-"))), collapse = "-"))
nnn_ordered_name <- sapply(1:length(n_pair_wise_combination$SNPComb), function(x) paste0(sort(unlist(strsplit(n_pair_wise_combination$SNPComb[x], "-"))), collapse = "-"))

a_overlapped_combination_ordered <- intersect(a_ordered_name, aaa_ordered_name) # 190
n_overlapped_combination_ordered <- intersect(n_ordered_name, nnn_ordered_name) # 190

## verification of overlapping MiHA combination between the HLI cohort and GWAS cohort
aind <- which(!(aaa_ordered_name %in% a_overlapped_combination_ordered))
aind2 <- which(!(a_ordered_name %in% a_overlapped_combination_ordered))

a_num_pairs <- sapply(1:length(a_joint_SNPs$X), function(x) length(which((unlist(strsplit(a_joint_SNPs$X[x], "-"))) %in% MiHA_SNPs$SNP_ID))) # length(which(a_num_pairs == 2)): 190
n_num_pairs <- sapply(1:length(n_joint_SNPs$X), function(x) length(which((unlist(strsplit(n_joint_SNPs$X[x], "-"))) %in% MiHA_SNPs$SNP_ID))) # length(which(n_num_pairs == 2)): 190

########## verification ends here ####################

####### combine two tables

## aGVHD
a_GWAS_ind <- which((aaa_ordered_name %in% a_overlapped_combination_ordered))
a_HLI_ind <- which((a_ordered_name %in% a_overlapped_combination_ordered))

a_combined_table <- a_pair_wise_combination[a_GWAS_ind, ]
# abs-abs; pres-abs; abs-pres; pres-pres
num_overlapped <- length(a_GWAS_ind)
for(id in 1:num_overlapped){
  
  index <- which(a_ordered_name %in% aaa_ordered_name[a_GWAS_ind[id]])
  
  a_combined_table$`abs-abs`[id] <- a_combined_table$`abs-abs`[id] + a_joint_SNPs$absence.absence[index]
  a_combined_table$`pres-abs`[id] <- a_combined_table$`pres-abs`[id] + a_joint_snps_ordered$presence.HET..absence[index] + a_joint_snps_ordered$presence.HOM..absence[index]
  a_combined_table$`abs-pres`[id] <- a_combined_table$`abs-pres`[id] + a_joint_snps_ordered$absence.presence.HET.[index] + a_joint_snps_ordered$absence.presence.HOM.[index]
  a_combined_table$`pres-pres`[id] <- a_combined_table$`pres-pres`[id] + a_joint_snps_ordered$presence.HET..presence.HOM.[index] + a_joint_snps_ordered$presence.HOM..presence.HET.[index]
  
}
a_combined_table <- rbind(a_combined_table, a_pair_wise_combination[-a_GWAS_ind, ])

# nGVHD
n_GWAS_ind <- which((nnn_ordered_name %in% n_overlapped_combination_ordered))
n_HLI_ind <- which((n_ordered_name %in% n_overlapped_combination_ordered))

n_combined_table <- n_pair_wise_combination[n_GWAS_ind, ]
# abs-abs; pres-abs; abs-pres; pres-pres
num_overlapped <- length(n_GWAS_ind)
for(id in 1:num_overlapped){
  
  index <- which(n_ordered_name %in% nnn_ordered_name[n_GWAS_ind[id]])
  
  n_combined_table$`abs-abs`[id] <- n_combined_table$`abs-abs`[id] + n_joint_SNPs$absence.absence[index]
  n_combined_table$`pres-abs`[id] <- n_combined_table$`pres-abs`[id] + n_joint_snps_ordered$presence.HET..absence[index] + n_joint_snps_ordered$presence.HOM..absence[index]
  n_combined_table$`abs-pres`[id] <- n_combined_table$`abs-pres`[id] + n_joint_snps_ordered$absence.presence.HET.[index] + n_joint_snps_ordered$absence.presence.HOM.[index]
  n_combined_table$`pres-pres`[id] <- n_combined_table$`pres-pres`[id] + n_joint_snps_ordered$presence.HET..presence.HOM.[index] + n_joint_snps_ordered$presence.HOM..presence.HET.[index]
  
}
n_combined_table <- rbind(n_combined_table, n_pair_wise_combination[-n_GWAS_ind, ])


##################
# Statistic Test
##################
Num_joint <- dim(a_combined_table)[1]
p_values <- vector(mode = "numeric", length = Num_joint)
names(p_values) <- as.character(a_combined_table$SNPComb)
for(id in 1:Num_joint){
  
  rsNames <- gsub("-", ".", as.character(a_combined_table$SNPComb[id]))
  # rsNames <- unlist(strsplit(as.character(a_joint_SNPs$X[id]),"-"))
  # rs1_index <- which(grepl(rsNames[1], names(a_individual_SNP_freq)))
  # rs2_index <- which(grepl(rsNames[2], names(a_individual_SNP_freq)))
  # 
  # a_individual_SNP_freq[rs1_index]
  # a_individual_SNP_freq[rs2_index]
  # a_individual_SNP_freq[rs1_index] * a_individual_SNP_freq[rs2_index]
  
  #cat(a_joint_SNPs$presence.HET..presence.HET.[id]) #/dim(aSNP_mat)[1]
  if(a_combined_table$`pres-pres`[id] >0 || n_combined_table$`pres-pres`[id]>0){
    cat(id, "\n")
    eval(parse(text = paste0("test_matrix <- matrix(c(
                             a_combined_table$`pres-abs`[id],
                             n_combined_table$`pres-abs`[id],
                             a_combined_table$`abs-pres`[id],
                             n_combined_table$`abs-pres`[id],
                             a_combined_table$`pres-pres`[id],
                             n_combined_table$`pres-pres`[id]),
                             nrow = 2,
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"presence-absence\", \"absence-presence\", \"presence-presence\")
                             ))")))
  }else{
    eval(parse(text = paste0("test_matrix <- matrix(c( a_combined_table$`pres-abs`[id],
                             n_combined_table$`pres-abs`[id],
                             a_combined_table$`abs-pres`[id],
                             n_combined_table$`abs-pres`[id]),
                             nrow = 2,
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"presence-absence\", \"absence-presence\")
                             ))")))
  }
  
  
  # chisq.test(chi_mat)$statistic
  # chisq.test(chi_mat)$expected
  # 
  p_values[id] <- fisher.test(test_matrix)$"p.value"
  
  }

ordered_p_names <- names(p_values[order(p_values, decreasing = F)])
##########
# print to a file
Num_joint <- dim(a_combined_table)[1]
sink(file = "../ClinVar/All_disese_GWASH_HLI_combined_unrestricted_Fisher_test_pvalue.txt", append = F)
ordered_ID <- order(p_values, decreasing = F)
for(iid in 1:Num_joint){
  
  id <- ordered_ID[iid]
  rsNames <- gsub("-", ".", as.character(a_combined_table$SNPComb[id]))
  # rsNames <- unlist(strsplit(as.character(a_joint_SNPs$X[id]),"-"))
  # rs1_index <- which(grepl(rsNames[1], names(a_individual_SNP_freq)))
  # rs2_index <- which(grepl(rsNames[2], names(a_individual_SNP_freq)))
  # 
  # a_individual_SNP_freq[rs1_index]
  # a_individual_SNP_freq[rs2_index]
  # a_individual_SNP_freq[rs1_index] * a_individual_SNP_freq[rs2_index]
  
  # ordered_index <- which(a_joint_snps_ordered$X == a_joint_SNPs$X[id])
  
  #cat(a_joint_SNPs$presence.HET..presence.HET.[id]) #/dim(aSNP_mat)[1]
  if(a_combined_table$`pres-pres`[id] >0 || n_combined_table$`pres-pres`[id]>0){
    cat("ID = ", id, "\n")
    eval(parse(text = paste0("test_matrix <- matrix(c(
                             a_combined_table$`pres-abs`[id],
                             n_combined_table$`pres-abs`[id],
                             a_combined_table$`abs-pres`[id],
                             n_combined_table$`abs-pres`[id],
                             a_combined_table$`pres-pres`[id],
                             n_combined_table$`pres-pres`[id]), 
                             nrow = 2, 
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"presence-absence\", \"absence-presence\", \"presence-presence\")
                             ))")))
  }else{
    cat("ID = ", id, "\n")
    eval(parse(text = paste0("test_matrix <- matrix(c(
                             a_combined_table$`pres-abs`[id],
                             n_combined_table$`pres-abs`[id],
                             a_combined_table$`abs-pres`[id],
                             n_combined_table$`abs-pres`[id]), 
                             nrow = 2, 
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"presence-absence\", \"absence-presence\")
                             ))")))
  }
  
  print(test_matrix)
  cat("\n Fisher's exact test (two.sided) p-value: ")
  cat(p_values[id])
  cat("\n====================================\n")
  
  }
sink()

#### HLI only
####### combine two tables

## aGVHD
a_GWAS_ind <- which((aaa_ordered_name %in% a_overlapped_combination_ordered))
a_HLI_ind <- which((a_ordered_name %in% a_overlapped_combination_ordered))

a_combined_table <- a_pair_wise_combination[a_GWAS_ind, ]
# abs-abs; pres-abs; abs-pres; pres-pres
num_overlapped <- length(a_GWAS_ind)
for(id in 1:num_overlapped){
  
  index <- which(a_ordered_name %in% aaa_ordered_name[a_GWAS_ind[id]])
  
  a_combined_table$`abs-abs`[id] <-  a_joint_SNPs$absence.absence[index]
  a_combined_table$`pres-abs`[id] <- a_joint_snps_ordered$presence.HET..absence[index] + a_joint_snps_ordered$presence.HOM..absence[index]
  a_combined_table$`abs-pres`[id] <-  a_joint_snps_ordered$absence.presence.HET.[index] + a_joint_snps_ordered$absence.presence.HOM.[index]
  a_combined_table$`pres-pres`[id] <- a_joint_snps_ordered$presence.HET..presence.HOM.[index] + a_joint_snps_ordered$presence.HOM..presence.HET.[index]
  
}
# a_combined_table <- rbind(a_combined_table, a_pair_wise_combination[-a_GWAS_ind, ])

# nGVHD
n_GWAS_ind <- which((nnn_ordered_name %in% n_overlapped_combination_ordered))
n_HLI_ind <- which((n_ordered_name %in% n_overlapped_combination_ordered))

n_combined_table <- n_pair_wise_combination[n_GWAS_ind, ]
# abs-abs; pres-abs; abs-pres; pres-pres
num_overlapped <- length(n_GWAS_ind)
for(id in 1:num_overlapped){
  
  index <- which(n_ordered_name %in% nnn_ordered_name[n_GWAS_ind[id]])
  
  n_combined_table$`abs-abs`[id] <-  n_joint_SNPs$absence.absence[index]
  n_combined_table$`pres-abs`[id] <- n_joint_snps_ordered$presence.HET..absence[index] + n_joint_snps_ordered$presence.HOM..absence[index]
  n_combined_table$`abs-pres`[id] <- n_joint_snps_ordered$absence.presence.HET.[index] + n_joint_snps_ordered$absence.presence.HOM.[index]
  n_combined_table$`pres-pres`[id] <- n_joint_snps_ordered$presence.HET..presence.HOM.[index] + n_joint_snps_ordered$presence.HOM..presence.HET.[index]
  
}
# n_combined_table <- rbind(n_combined_table, n_pair_wise_combination[-n_GWAS_ind, ])

##################
# HLI only Statistic Test
##################
Num_joint <- dim(a_combined_table)[1]
p_values <- vector(mode = "numeric", length = Num_joint)
names(p_values) <- as.character(a_combined_table$SNPComb)
for(id in 1:Num_joint){
  
  rsNames <- gsub("-", ".", as.character(a_combined_table$SNPComb[id]))
  # rsNames <- unlist(strsplit(as.character(a_joint_SNPs$X[id]),"-"))
  # rs1_index <- which(grepl(rsNames[1], names(a_individual_SNP_freq)))
  # rs2_index <- which(grepl(rsNames[2], names(a_individual_SNP_freq)))
  # 
  # a_individual_SNP_freq[rs1_index]
  # a_individual_SNP_freq[rs2_index]
  # a_individual_SNP_freq[rs1_index] * a_individual_SNP_freq[rs2_index]
  
  #cat(a_joint_SNPs$presence.HET..presence.HET.[id]) #/dim(aSNP_mat)[1]
  if(a_combined_table$`pres-pres`[id] >0 || n_combined_table$`pres-pres`[id]>0){
    cat("ID = ", id, "\n")
    eval(parse(text = paste0("test_matrix <- matrix(c(
                             a_combined_table$`pres-abs`[id],
                             n_combined_table$`pres-abs`[id],
                             a_combined_table$`abs-pres`[id],
                             n_combined_table$`abs-pres`[id],
                             a_combined_table$`pres-pres`[id],
                             n_combined_table$`pres-pres`[id]), 
                             nrow = 2, 
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"presence-absence\", \"absence-presence\", \"presence-presence\")
                             ))")))
  }else{
    cat("ID = ", id, "\n")
    eval(parse(text = paste0("test_matrix <- matrix(c(
                             a_combined_table$`pres-abs`[id],
                             n_combined_table$`pres-abs`[id],
                             a_combined_table$`abs-pres`[id],
                             n_combined_table$`abs-pres`[id]), 
                             nrow = 2, 
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"presence-absence\", \"absence-presence\")
                             ))")))
  }
  
  # chisq.test(chi_mat)$statistic
  # chisq.test(chi_mat)$expected
  # 
  p_values[id] <- fisher.test(test_matrix)$"p.value"
  
  
  
  }

p_values[order(p_values, decreasing = F)]
##########
# print to a file
Num_joint <- dim(a_combined_table)[1]
sink(file = "../ClinVar/All_disease_HLI_only_combined_unrestricted_Fisher_test_pvalue.txt", append = F)
# ordered_ID <- order(p_values, decreasing = F)
# ordered_ID2 <- ordered_ID
ordered_ID <- unlist(sapply(1:length(ordered_p_names), function(x) which(a_combined_table$SNPComb %in% ordered_p_names[x])))
for(iid in 1:Num_joint){
  
  id <- ordered_ID[iid]
  if(id <= Num_joint){
    rsNames <- gsub("-", ".", as.character(a_combined_table$SNPComb[id]))
    # rsNames <- unlist(strsplit(as.character(a_joint_SNPs$X[id]),"-"))
    # rs1_index <- which(grepl(rsNames[1], names(a_individual_SNP_freq)))
    # rs2_index <- which(grepl(rsNames[2], names(a_individual_SNP_freq)))
    # 
    # a_individual_SNP_freq[rs1_index]
    # a_individual_SNP_freq[rs2_index]
    # a_individual_SNP_freq[rs1_index] * a_individual_SNP_freq[rs2_index]
    
    # ordered_index <- which(a_joint_snps_ordered$X == a_joint_SNPs$X[id])
    
    #cat(a_joint_SNPs$presence.HET..presence.HET.[id]) #/dim(aSNP_mat)[1]
    if(a_combined_table$`pres-pres`[id] >0 || n_combined_table$`pres-pres`[id]>0){
      cat("ID = ", id, "\n")
      eval(parse(text = paste0("test_matrix <- matrix(c(
                               a_combined_table$`pres-abs`[id],
                               n_combined_table$`pres-abs`[id],
                               a_combined_table$`abs-pres`[id],
                               n_combined_table$`abs-pres`[id],
                               a_combined_table$`pres-pres`[id],
                               n_combined_table$`pres-pres`[id]), 
                               nrow = 2, 
                               dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                               = c(\"presence-absence\", \"absence-presence\", \"presence-presence\")
                               ))")))
    }else{
      cat("ID = ", id, "\n")
      eval(parse(text = paste0("test_matrix <- matrix(c(
                               a_combined_table$`pres-abs`[id],
                               n_combined_table$`pres-abs`[id],
                               a_combined_table$`abs-pres`[id],
                               n_combined_table$`abs-pres`[id]), 
                               nrow = 2, 
                               dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                               = c(\"presence-absence\", \"absence-presence\")
                               ))")))
    }
    
    print(test_matrix)
    cat("\n Fisher's exact test (two.sided) p-value: ")
    cat(p_values[id])
    cat("\n====================================\n")
    
    }
  
  
}
sink()

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
KnownMiHA_table$HLA_SNP <- sapply(1:dim(KnownMiHA_table)[1], function(x) paste0(KnownMiHA_table$HLA[x], "-", KnownMiHA_table$rs.number[x]))
unique_MiHAs <- as.character(unique(GWASH_single_MiHA_stats$HLA_SNP))
num_MiHAs <- length(unique_MiHAs)
LLR_MiHAs <- data.frame(name = character(num_MiHAs),
                        hla = character(num_MiHAs),
                        LLR = numeric(num_MiHAs),
                        stringsAsFactors = F)
# all_MiHA_SNP$MiHAs <- as.character(all_MiHA_SNP$MiHAs)
# all_MiHA_SNP$Gene <- as.character(all_MiHA_SNP$Gene)
# all_MiHA_SNP$SNPs <- as.character(all_MiHA_SNP$SNPs)
for(miha_id in 1:num_MiHAs){
  
  index <- which(GWASH_single_MiHA_stats$HLA_SNP %in% unique_MiHAs[miha_id])
  
  aGVHD_count <- GWASH_single_MiHA_stats[index[which(GWASH_single_MiHA_stats[index, "GroupType"] == "aGVHD")], "Freq"]
  nGVHD_count <- abs(GWASH_single_MiHA_stats[index[which(GWASH_single_MiHA_stats[index, "GroupType"] == "non-aGVHD")], "Freq"])
  
  iid <- which(KnownMiHA_table$HLA_SNP %in% GWASH_single_MiHA_stats$HLA_SNP[index])
  Gene_MiHA_name <- paste0(as.character(KnownMiHA_table[iid, c("Gene", "KnownMiHA", "rs.number")]), collapse = " <> ")
  
  if(aGVHD_count > 0 | nGVHD_count >0){
    
    LLR_MiHAs$LLR[miha_id] <- log10(aGVHD_count/nGVHD_count)
    LLR_MiHAs$name[miha_id] <- Gene_MiHA_name
  }
  
}

mc_LLR <- LLR_MiHAs[order(LLR_MiHAs$LLR, decreasing = TRUE), ]
mc_LLR <- within(mc_LLR, name <- factor(name, levels=factor(mc_LLR$name)))

pp1 <- ggplot(GWASH_single_MiHA_stats, aes(x = HLA_SNP, ymax = Freq, ymin = 0, color = GroupType)) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  ggtitle("Known MiHAs") +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pp2 <- ggplot(mc_LLR, aes(x = name, ymax = LLR, ymin = 0, color = "#D55E00")) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  # ggtitle("Known MiHAs") +
  theme_bw() +
  xlab("") +
  coord_flip() +
  # scale_colour(values = "#D55E00") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") 

## 
index <- which(single_MiHA_stats$HLA_SNP %in% unique_MiHAs)

HLI_single_MiHA_stats <- single_MiHA_stats[index, ]

# num_MiHAs <- length(unique_MiHAs)
HLI_LLR_MiHAs <- data.frame(name = character(num_MiHAs),
                            hla = character(num_MiHAs),
                            LLR = numeric(num_MiHAs),
                            stringsAsFactors = F)

for(miha_id in 1:num_MiHAs){
  
  index <- which(HLI_single_MiHA_stats$HLA_SNP %in% unique_MiHAs[miha_id])
  if(length(index) > 0){
    aGVHD_count <- HLI_single_MiHA_stats[index[which(HLI_single_MiHA_stats[index, "GroupType"] == "a")], "Freq"]
    nGVHD_count <- abs(HLI_single_MiHA_stats[index[which(HLI_single_MiHA_stats[index, "GroupType"] == "n")], "Freq"])
    
    iid <- which(KnownMiHA_table$HLA_SNP %in% HLI_single_MiHA_stats$HLA_SNP[index])
    Gene_MiHA_name <- paste0(as.character(KnownMiHA_table[iid, c("Gene", "KnownMiHA", "rs.number")]), collapse = " <> ")
    
    if(aGVHD_count > 0 | nGVHD_count >0){
      
      HLI_LLR_MiHAs$LLR[miha_id] <- log10(aGVHD_count/nGVHD_count)
      HLI_LLR_MiHAs$name[miha_id] <- Gene_MiHA_name
    }
  }
}

HLI_mc_LLR <- HLI_LLR_MiHAs[order(HLI_LLR_MiHAs$LLR, decreasing = TRUE), ]
HLI_mc_LLR <- within(HLI_mc_LLR, name <- factor(name, levels=factor(HLI_mc_LLR$name)))

pp1 <- ggplot(HLI_single_MiHA_stats, aes(x = HLA_SNP, ymax = Freq, ymin = 0, color = GroupType)) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  ggtitle("Known MiHAs") +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

pp2 <- ggplot(HLI_mc_LLR, aes(x = name, ymax = LLR, ymin = 0, color = "#D55E00")) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  # ggtitle("Known MiHAs") +
  theme_bw() +
  xlab("") +
  coord_flip() +
  # scale_colour(values = "#D55E00") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") 


###### save HLA restricted MiHA table
num_MiHAs <- dim(GWASH_single_MiHA_stats)[1]/2
GWASH_single_MiHA_table <- data.frame(Gene = character(num_MiHAs), 
                                      MiHA = character(num_MiHAs), 
                                      rsNumber = character(num_MiHAs), 
                                      HLA = character(num_MiHAs),
                                      aGVHD = numeric(num_MiHAs), 
                                      nGVHD = numeric(num_MiHAs), 
                                      LLR10_avg = numeric(num_MiHAs),
                                      stringsAsFactors = F)

for(miha_id in 1:num_MiHAs){
  
  index <- which(GWASH_single_MiHA_stats$HLA_SNP %in% unique_MiHAs[miha_id])
  
  GWASH_single_MiHA_table$aGVHD[miha_id] <- GWASH_single_MiHA_stats[index[which(GWASH_single_MiHA_stats[index, "GroupType"] == "aGVHD")], "Freq"]
  GWASH_single_MiHA_table$nGVHD[miha_id] <- abs(GWASH_single_MiHA_stats[index[which(GWASH_single_MiHA_stats[index, "GroupType"] == "non-aGVHD")], "Freq"])
  
  iid <- which(KnownMiHA_table$HLA_SNP %in% GWASH_single_MiHA_stats$HLA_SNP[index])
  # Gene_MiHA_name <- paste0(as.character(KnownMiHA_table[iid, c("Gene", "KnownMiHA", "rs.number")]), collapse = " <> ")
  GWASH_single_MiHA_table$Gene[miha_id] <- paste0(unique(KnownMiHA_table[iid, "Gene"]), collapse = "/")
  GWASH_single_MiHA_table$MiHA[miha_id] <- paste0(unique(KnownMiHA_table[iid, "KnownMiHA"]), collapse = "/")
  GWASH_single_MiHA_table$rsNumber[miha_id] <- paste0(unique(KnownMiHA_table[iid, "rs.number"]), collapse = "/")
  GWASH_single_MiHA_table$HLA[miha_id] <- paste0(unique(KnownMiHA_table[iid, "HLA"]), collapse = "/")
  
  if(aGVHD_count > 0 | nGVHD_count >0){
    
    GWASH_single_MiHA_table$LLR10_avg[miha_id] <- log10(GWASH_single_MiHA_table$aGVHD[miha_id]/GWASH_single_MiHA_table$nGVHD[miha_id])
    # GWASH_single_MiHA_table$name[miha_id] <- Gene_MiHA_name
  }
  
}

write.csv(GWASH_single_MiHA_table, file = "../ClinVar/GWASH/GWASH_cohort_Restricted_MiHA_summary.csv", row.names = F)

num_MiHAs <- dim(GWASH_single_MiHA_stats)[1]/2
HLI_single_MiHA_table <- data.frame(Gene = character(num_MiHAs), 
                                    MiHA = character(num_MiHAs), 
                                    rsNumber = character(num_MiHAs), 
                                    HLA = character(num_MiHAs),
                                    aGVHD = numeric(num_MiHAs), 
                                    nGVHD = numeric(num_MiHAs), 
                                    LLR10_avg = numeric(num_MiHAs),
                                    stringsAsFactors = F)

for(miha_id in 1:num_MiHAs){
  
  index <- which(HLI_single_MiHA_stats$HLA_SNP %in% unique_MiHAs[miha_id])
  if(length(index) > 0){
    HLI_single_MiHA_table$aGVHD[miha_id] <- HLI_single_MiHA_stats[index[which(HLI_single_MiHA_stats[index, "GroupType"] == "a")], "Freq"]
    HLI_single_MiHA_table$nGVHD[miha_id] <- abs(HLI_single_MiHA_stats[index[which(HLI_single_MiHA_stats[index, "GroupType"] == "n")], "Freq"])
    
    iid <- which(KnownMiHA_table$HLA_SNP %in% HLI_single_MiHA_stats$HLA_SNP[index])
    # Gene_MiHA_name <- paste0(as.character(KnownMiHA_table[iid, c("Gene", "KnownMiHA", "rs.number")]), collapse = " <> ")
    HLI_single_MiHA_table$Gene[miha_id] <- paste0(unique(KnownMiHA_table[iid, "Gene"]), collapse = "/")
    HLI_single_MiHA_table$MiHA[miha_id] <- paste0(unique(KnownMiHA_table[iid, "KnownMiHA"]), collapse = "/")
    HLI_single_MiHA_table$rsNumber[miha_id] <- paste0(unique(KnownMiHA_table[iid, "rs.number"]), collapse = "/")
    HLI_single_MiHA_table$HLA[miha_id] <- paste0(unique(KnownMiHA_table[iid, "HLA"]), collapse = "/")
    
    if(aGVHD_count > 0 | nGVHD_count >0){
      
      HLI_single_MiHA_table$LLR10_avg[miha_id] <- log10(HLI_single_MiHA_table$aGVHD[miha_id]/HLI_single_MiHA_table$nGVHD[miha_id])
      # GWASH_single_MiHA_table$name[miha_id] <- Gene_MiHA_name
    }
  }
  
}

HLI_single_MiHA_table <- HLI_single_MiHA_table[which(HLI_single_MiHA_table$LLR10_avg != 0), ]

write.csv(HLI_single_MiHA_table, file = "../ClinVar/GWASH/HLI_cohort_Overlapping_Restricted_MiHA_summary.csv", row.names = F)

######### HLI GWASH overlapping cases
HLA_typings <- read.csv(file = "../HLI_hla_mg_v3.csv")

load("../Data/ID_table.RData")

available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]

avaialble_IDs_HLA_typing <- HLA_typings[(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID),]
num_groups <- dim(avaialble_IDs_HLA_typing)[1]

overlapping_case_IDs <- intersect(GWAS_sample_table$bmt_case, avaialble_IDs_HLA_typing$bmt_case_num)
## 69 cases

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
Recipients <- which(ID_table$subjectType == "R")
Donors <- which(ID_table$subjectType == "D")

ID_table$CaseNum <- numeric(dim(ID_table)[1]) 

ID_table$CaseNum[Recipients] <- unlist(sapply(1:length(Recipients), function(x) HLA_typings$bmt_case_num[which(HLA_typings$nmdp_rid == ID_table$R_D_ID[Recipients[x]])]))
ID_table$CaseNum[Donors] <- unlist(sapply(1:length(Donors), function(x) if(length(which(HLA_typings$nmdp_id == ID_table$R_D_ID[Donors[x]])) == 1) HLA_typings$bmt_case_num[which(HLA_typings$nmdp_id == ID_table$R_D_ID[Donors[x]])] 
                                          else 0))

over_HLI_groupID <- unique(ID_table$GroupID[which(ID_table$CaseNum %in% overlapping_case_IDs)])

non_overlapping_HLI_cases <- SNP_mat[-which(SNP_mat$GroupID %in% over_HLI_groupID), ]
write.csv(non_overlapping_HLI_cases, file = "../ClinVar/GWASH/HLI_Non-overlapping_cohort_count_table.csv", row.names = F)



known_miha_freq <- read.delim("../WW_MiHA/Restricted_known_MiHAs.txt", header = F)
colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")
table(known_miha_freq[c("HLA_type","SNP")])

known_miha_freq$HLA_SNP <- sapply(1:dim(known_miha_freq)[1], 
                                  function(x) paste0(known_miha_freq$HLA_type[x], "-", known_miha_freq$SNP[x]))
# overlap_known_miha_freq <- known_miha_freq[which(known_miha_freq$GroupID %in% over_HLI_groupID), ]

overlap_known_miha_freq <- known_miha_freq[-which(known_miha_freq$GroupID %in% over_HLI_groupID), ] # non-overlapping
overlap_single_MiHA_stats <- as.data.frame(table(overlap_known_miha_freq[c("GroupType", "HLA_SNP")]))

index <- which(overlap_single_MiHA_stats$HLA_SNP %in% unique_MiHAs)

overlap_HLI_single_MiHA_stats <- overlap_single_MiHA_stats[index, ]

num_MiHAs <- dim(GWASH_single_MiHA_stats)[1]/2
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
    
    if(aGVHD_count > 0 | nGVHD_count >0){
      
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

