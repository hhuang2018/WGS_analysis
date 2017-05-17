source("util.R")
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(ggplot2)
GWAS_samp_table_fp <- "../ClinVar/GWASH/Metadata/newfnldataib1203.csv"
GWAS_file_fp <- "../ClinVar/GWASH/MIHA/data/"

GWAS_sample_table <-read.csv(GWAS_samp_table_fp, header = T, stringsAsFactors = F)

MiHA_SNPs <- read.table(paste0(GWAS_file_fp, "MiHAgeno.map"), stringsAsFactors = F)
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
  
}

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

for(id in 1:Num_AML_cases){
  
  index <- which(MiHA_genotypes$case_id %in% AML_cases_table$bmt_case[id])
  
  if(length(index) == 2){ # matching cases
    
    mismatch_MiHA_table$CaseNumber[id] <- AML_cases_table$bmt_case[id]
    mismatch_MiHA_table$aGVHD[id] <- AML_cases_table$agvhi24[id]
    
    for(jd in 1:num_SNPs){
      
      mismatch_MiHA_table[id, 4+2*(jd-1)] <- identical(MiHA_genotypes[index[1], SNP_index[1+2*(jd-1)]], MiHA_genotypes[index[2], SNP_index[1+2*(jd-1)]]) + 
        identical(MiHA_genotypes[index[1], SNP_index[2*jd]], MiHA_genotypes[index[2], SNP_index[2*jd]])
      
      mismatch_MiHA_table[id, 3+2*(jd-1)] <- as.integer(identical(MiHA_genotypes[index[1], SNP_index[1+2*(jd-1)]], MiHA_genotypes[index[2], SNP_index[1+2*(jd-1)]]) | 
                                                          identical(MiHA_genotypes[index[1], SNP_index[2*jd]], MiHA_genotypes[index[2], SNP_index[2*jd]]))
      
    }
    
    
  }
  
  
}

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

KnownMiHA_table <- read.csv("../ClinVar/Data/KnowMiHA_db2tab.csv")    






