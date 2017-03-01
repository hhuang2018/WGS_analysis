################
# counts table
################
# library(corrgram)
load("../Data/ID_table.RData")

# known_miha_freq <- read.csv("../WW_MiHA/HLA_restricted_knownMiHAs.csv", header = F)
# known_miha_freq <- read.delim("../WW_MiHA/Restricted_known_MiHAs.txt", header = F)
# colnames(known_miha_freq) <- c("GroupType", "GroupID",  "HLA_type", "SNP", "CHROM", "REF", "ALT")

known_miha_freq <- read.csv("../ClinVar//unRestrictedMiHAs.csv", header = F)
# colnames(known_miha_freq) <- c("GroupType", "GroupID", "HLA_type", "SNP", "CHROM", "REF", "ALT", "GENE", "Peptide")  # restricted table
colnames(known_miha_freq) <- c("GroupType", "GroupID", "HLA_type", "SNP", "CHROM", "Donor", "Recipient") # unristricted table

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

##########
# correlation_mat <- cor(SNP_mat[, -c(1,2)], use = "pairwise.complete.obs", method = "pearson")
# 

SNP_mat[SNP_mat == 2] <- 1
aSNP_mat <- SNP_mat[which(SNP_mat$GroupType == "a"), ]
# aSNP_mat <- aSNP_mat[which(rowSums(aSNP_mat[, -c(1,2)]) > 0), ]
a_correlation_mat <- cor(aSNP_mat[, -c(1,2)], use = "pairwise.complete.obs", method = "pearson")
a_correlation_mat[upper.tri(a_correlation_mat, diag = T)] <- 0
# which(a_correlation_mat == max(a_correlation_mat, na.rm = T), arr.ind = T)

a_snp_num <- dim(aSNP_mat)[2] - 2
a_ind_num <- dim(aSNP_mat)[1]
a_count_table <- matrix(data = NA, nrow = 1, ncol = 1)
for(id in 1:(a_snp_num-1)){
  for(jd in (id+1):a_snp_num){
    
    ind_SNP_mat <- matrix(data = 0, nrow = a_ind_num, ncol = 3)
    Type_name <- paste0(sapply(colnames(aSNP_mat)[c(id, jd)+2], function(x) unlist(strsplit(x, "-"))[2]), collapse = "-")
    
    ind_SNP_mat <- cbind(aSNP_mat[, c(id, jd)+2], rowSums(aSNP_mat[, c(id, jd)+2]))
    colnames(ind_SNP_mat) <- c(colnames(aSNP_mat)[c(id, jd)+2], Type_name)
    
    if(dim(a_count_table)[1] == 1){
      a_count_table <- ind_SNP_mat
    } else {
      a_count_table <- cbind(a_count_table, ind_SNP_mat)
    }
    
  }
}


nSNP_mat <- SNP_mat[which(SNP_mat$GroupType == "n"), ]
n_correlation_mat <- cor(SNP_mat[which(SNP_mat$GroupType == "n"), -c(1,2)], use = "pairwise.complete.obs", method = "pearson")
n_correlation_mat[upper.tri(n_correlation_mat, diag = T)] <- 0
which(n_correlation_mat == max(n_correlation_mat, na.rm = T), arr.ind = T)

n_snp_num <- dim(nSNP_mat)[2] - 2
n_ind_num <- dim(nSNP_mat)[1]
n_count_table <- matrix(data = NA, nrow = 1, ncol = 1)
for(id in 1:(n_snp_num-1)){
  for(jd in (id+1):n_snp_num){
    
    ind_SNP_mat <- matrix(data = 0, nrow = n_ind_num, ncol = 3)
    Type_name <- paste0(sapply(colnames(nSNP_mat)[c(id, jd)+2], function(x) unlist(strsplit(x, "-"))[2]), collapse = "-")
    
    ind_SNP_mat <- cbind(nSNP_mat[, c(id, jd)+2], rowSums(nSNP_mat[, c(id, jd)+2]))
    colnames(ind_SNP_mat) <- c(colnames(nSNP_mat)[c(id, jd)+2], Type_name)
    
    if(dim(n_count_table)[1] == 1){
      n_count_table <- ind_SNP_mat
    } else {
      n_count_table <- cbind(n_count_table, ind_SNP_mat)
    }
    
  }
}

a_pairwise_SNPs <- a_count_table[, seq(from = 3, to = dim(a_count_table)[2], by = 3)] # 528 combinations
n_pairwise_SNPs <- n_count_table[, seq(from = 3, to = dim(n_count_table)[2], by = 3)] # 528 combinations

# identical(colnames(a_pairwise_SNPs), colnames(n_pairwise_SNPs) )
num_col <- dim(a_pairwise_SNPs)[2]
a_pairwise_tab <- matrix(data = 0, nrow = 6, ncol = num_col)
n_pairwise_tab <- matrix(data = 0, nrow = 6, ncol = num_col)
colnames(a_pairwise_tab) <- colnames(a_pairwise_SNPs)
colnames(n_pairwise_tab) <- colnames(n_pairwise_SNPs)

# sum - index
#   0 -     1: absence-absence
#   1 -     2: absence-presence or presence-absence (het)
#   2 -     3: absence-presence (hom)
#   2 -     4: presence-presence (het) == 2
#   3 -     5: presence-presence (hom+het) == 3
#   4 -     6: presence-presence (hom) == 4
category_names <- c("absence-absence", "absence-presence(HET)",
                    "absence-presence(HOM)", "presence(HET)-presence(HET)",
                    "presence(HOM)-presence(HET)", "presence(HOM)-presence(HOM)")

index <- seq(from = 1, to = dim(a_count_table)[2], by = 3)

for(id in 1:num_col){
  #### aGVHD groups
  a_tab <- as.data.frame(table(a_pairwise_SNPs[, id]), stringsAsFactors = F)
  a_ind_0 <- which(a_tab$Var1 == 0)
  a_ind_1 <- which(a_tab$Var1 == 1)
  a_ind_2 <- which(a_tab$Var1 == 2)
  a_ind_3 <- which(a_tab$Var1 == 3)
  a_ind_4 <- which(a_tab$Var1 == 4)
  
  if(length(a_ind_0) > 0 && length(a_ind_1) > 0){ # both 0 and 1 exists
    a_pairwise_tab[c(1,2), id] <- a_tab$Freq[c(a_ind_0, a_ind_1)]
  }else if(length(a_ind_0) == 0 && length(a_ind_1) > 0){ # 0 doesn't exist
    a_pairwise_tab[1, id] <- 0
    a_pairwise_tab[2, id] <- a_tab$Freq[a_ind_1]
  }else if(length(a_ind_0) > 0 && length(a_ind_1) == 0){ # 1 doesn't exist
    a_pairwise_tab[1, id] <- a_tab$Freq[a_ind_0]
    a_pairwise_tab[2, id] <- 0
  }else{ # both don't exist
    a_pairwise_tab[1, id] <- 0
    a_pairwise_tab[2, id] <- 0
  }
  
  if(length(a_ind_2)> 0){ # 2 presenses
     
    sum_two_index <- which(a_pairwise_SNPs[,id] == 2)
    hom_index <- which(a_count_table[sum_two_index, index[id]:(index[id]+1)] == 2, arr.ind = T)
    a_pairwise_tab[3, id] <- dim(hom_index)[1] # homozygotes; presense - absence
    a_pairwise_tab[4, id] <- length(sum_two_index) - dim(hom_index)[1] # heterozygotes: presence - presence

  }
  
  if(length(a_ind_3)> 0){ # 3 presense
    a_pairwise_tab[5, id] <- a_tab$Freq[a_ind_3]
  }
  
  if(length(a_ind_4)> 0){ # 4 presense
    a_pairwise_tab[6, id] <- a_tab$Freq[a_ind_4]
  }
  
  ##### non-aGVHD groups
  # n_tab <- as.data.frame(table(n_pairwise_SNPs[, id]), stringsAsFactors = F)
  # n_pairwise_tab[as.numeric(n_tab$Var1)+1, id] <- n_tab$Freq
  
  n_tab <- as.data.frame(table(n_pairwise_SNPs[, id]), stringsAsFactors = F)
  n_ind_0 <- which(n_tab$Var1 == 0)
  n_ind_1 <- which(n_tab$Var1 == 1)
  n_ind_2 <- which(n_tab$Var1 == 2)
  n_ind_3 <- which(n_tab$Var1 == 3)
  n_ind_4 <- which(n_tab$Var1 == 4)
  
  if(length(n_ind_0) > 0 && length(n_ind_1) > 0){ # both 0 and 1 exists
    n_pairwise_tab[c(1,2), id] <- n_tab$Freq[c(n_ind_0, n_ind_1)]
  }else if(length(n_ind_0) == 0 && length(n_ind_1) > 0){ # 0 doesn't exist
    n_pairwise_tab[1, id] <- 0
    n_pairwise_tab[2, id] <- n_tab$Freq[n_ind_1]
  }else if(length(n_ind_0) > 0 && length(n_ind_1) == 0){ # 1 doesn't exist
    n_pairwise_tab[1, id] <- n_tab$Freq[n_ind_0]
    n_pairwise_tab[2, id] <- 0
  }else{ # both don't exist
    n_pairwise_tab[1, id] <- 0
    n_pairwise_tab[2, id] <- 0
  }
  
  if(length(n_ind_2)> 0){ # 2 presenses
    
    sum_two_index <- which(n_pairwise_SNPs[,id] == 2)
    hom_index <- which(n_count_table[sum_two_index, index[id]:(index[id]+1)] == 2, arr.ind = T)
    n_pairwise_tab[3, id] <- dim(hom_index)[1] # homozygotes; presense - absence
    n_pairwise_tab[4, id] <- length(sum_two_index) - dim(hom_index)[1] # heterozygotes: presence - presence
  
  }
  
  if(length(n_ind_3)> 0){ # 3 presense
    n_pairwise_tab[5, id] <- n_tab$Freq[n_ind_3]
  }
  
  if(length(n_ind_4)> 0){ # 4 presense
    n_pairwise_tab[6, id] <- n_tab$Freq[n_ind_4]
  }
  
}

rownames(a_pairwise_tab) <- category_names
rownames(n_pairwise_tab) <- category_names

aa <- t(a_pairwise_tab)
nn <- t(n_pairwise_tab)

write.csv(aa, file = "../ClinVar/Data/aGVHD_pairwise_summary_AP.csv")
write.csv(nn, file = "../ClinVar/Data/nonGVHD_pairwise_summary_AP.csv")

write.csv(aa, file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_summary_AP.csv")
write.csv(nn, file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_summary_AP.csv")
# a_ordered_index <- order(a_pairwise_tab[3, ], decreasing = T)
# a_combination <- colnames(a_pairwise_tab)[a_ordered_index]
# a_pairwise_tab[, a_ordered_index]
# 
# n_ordered_index <- order(n_pairwise_tab[2, ], decreasing = T)
# n_combination <- colnames(n_pairwise_tab)[n_ordered_index]
# n_pairwise_tab[, n_ordered_index[1:10]]

###########
# Ordered 
#  absence - presence(het): 1
#  presence(het) - absence: 2
#  absence - presence(hom): 3
#  presence(hom) - absence: 4
#  presence(het) - presence(hom): 5
#  presence(hom) - presence(het): 6

# sum - index
#   1 -     1: absence - presence(het)
#   1 -     2: presence(het) - absence
#   2 -     3: absence - presence(hom)
#   2 -     4: presence(hom) - absence
#   3 -     5: presence(het) -presence(hom)
#   3 -     6: presence(hom) -presence(het)
category_names2 <- c("absence-presence(HET)", "presence(HET)-absence",
                     "absence-presence(HOM)", "presence(HOM)-absence",
                     "presence(HET)-presence(HOM)", "presence(HOM)-presence(HET)")

index <- seq(from = 1, to = dim(a_count_table)[2], by = 3)
num_col <- dim(a_pairwise_SNPs)[2]
a_pairwise_ordered <- matrix(data = 0, nrow = 6, ncol = num_col)
n_pairwise_ordered <- matrix(data = 0, nrow = 6, ncol = num_col)
colnames(a_pairwise_ordered) <- colnames(a_pairwise_SNPs)
colnames(n_pairwise_ordered) <- colnames(n_pairwise_SNPs)

for(id in 1:num_col){
  #### aGVHD groups
  a_tab <- as.data.frame(table(a_pairwise_SNPs[, id]), stringsAsFactors = F)
  # a_ind_0 <- which(a_tab$Var1 == 0)
  a_ind_1 <- which(a_tab$Var1 == 1)
  a_ind_2 <- which(a_tab$Var1 == 2)
  a_ind_3 <- which(a_tab$Var1 == 3)
  # a_ind_4 <- which(a_tab$Var1 == 4)
  
  if(length(a_ind_1) > 0){ 
    
    sum_one_index <- which(a_pairwise_SNPs[,id] == 1)
    a_pairwise_ordered[1, id] <- sum(a_count_table[sum_one_index, index[id]+1]) # absence-presense(HET)
    a_pairwise_ordered[2, id] <- sum(a_count_table[sum_one_index, index[id]]) # presense(HET)-absence

  }
  
  if(length(a_ind_2)> 0){ # 2 presenses
    
    sum_two_index <- which(a_pairwise_SNPs[,id] == 2)
    hom_index <- which(a_count_table[sum_two_index, index[id]:(index[id]+1)] == 2, arr.ind = T)
    if(dim(hom_index)[1] > 0){

      a_pairwise_ordered[3, id] <- length(unique(hom_index[hom_index[, 2] == 2, 1])) # absence-presense(HOM)
      a_pairwise_ordered[4, id] <- length(unique(hom_index[hom_index[, 2] == 1, 1])) # presense(HOM)-absence
      
    }
  }
  
  if(length(a_ind_3)> 0){ # 3 presense
    
    sum_three_index <- which(a_pairwise_SNPs[,id] == 3)
    hom_index <- which(a_count_table[sum_two_index, index[id]:(index[id]+1)] == 2, arr.ind = T)
    
    a_pairwise_ordered[5, id] <- length(unique(hom_index[hom_index[, 2] == 2, 1])) # presence(HET) -presence(HOM)
    a_pairwise_ordered[6, id] <- length(unique(hom_index[hom_index[, 2] == 1, 1])) # presence(HOM) -presence(HET)
    
  }
  
  
  ##### non-aGVHD groups
  # n_tab <- as.data.frame(table(n_pairwise_SNPs[, id]), stringsAsFactors = F)
  # n_pairwise_tab[as.numeric(n_tab$Var1)+1, id] <- n_tab$Freq
  
  n_tab <- as.data.frame(table(n_pairwise_SNPs[, id]), stringsAsFactors = F)
  n_ind_1 <- which(n_tab$Var1 == 1)
  n_ind_2 <- which(n_tab$Var1 == 2)
  n_ind_3 <- which(n_tab$Var1 == 3)
  
  if(length(n_ind_1) > 0){ 
    
    sum_one_index <- which(n_pairwise_SNPs[,id] == 1)
    n_pairwise_ordered[1, id] <- sum(n_count_table[sum_one_index, index[id]+1]) # absence-presense(HET)
    n_pairwise_ordered[2, id] <- sum(n_count_table[sum_one_index, index[id]]) # presense(HET)-absence
    
  }
  
  if(length(n_ind_2)> 0){ # 2 presenses
    
    sum_two_index <- which(n_pairwise_SNPs[,id] == 2)
    hom_index <- which(n_count_table[sum_two_index, index[id]:(index[id]+1)] == 2, arr.ind = T)
    if(dim(hom_index)[1] > 0){
      
      n_pairwise_ordered[3, id] <- length(unique(hom_index[hom_index[, 2] == 2, 1])) # absence-presense(HOM)
      n_pairwise_ordered[4, id] <- length(unique(hom_index[hom_index[, 2] == 1, 1])) # presense(HOM)-absence
      
    }
  }
  
  if(length(n_ind_3)> 0){ # 3 presense
    
    sum_three_index <- which(n_pairwise_SNPs[,id] == 3)
    hom_index <- which(n_count_table[sum_two_index, index[id]:(index[id]+1)] == 2, arr.ind = T)
    
    n_pairwise_ordered[5, id] <- length(unique(hom_index[hom_index[, 2] == 2, 1])) # presence(HET) -presence(HOM)
    n_pairwise_ordered[6, id] <- length(unique(hom_index[hom_index[, 2] == 1, 1])) # presence(HOM) -presence(HET)
    
  }
}

rownames(a_pairwise_ordered) <- category_names2
rownames(n_pairwise_ordered) <- category_names2

aa <- t(a_pairwise_ordered)
nn <- t(n_pairwise_ordered)

write.csv(aa, file = "../ClinVar/Data/aGVHD_pairwise_ordered_AP.csv")
write.csv(nn, file = "../ClinVar/Data/nonGVHD_pairwise_ordered_AP.csv")

write.csv(aa, file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_ordered_AP.csv")
write.csv(nn, file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_ordered_AP.csv")

#### independency test
a_individual_SNP_freq <- colSums(aSNP_mat[,-c(1,2)])#/dim(aSNP_mat)[1]
n_individual_SNP_freq <- colSums(nSNP_mat[, -c(1,2)])#/dim(nSNP_mat)[1]

a_joint_SNPs <- read.csv(file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_summary_AP.csv")
n_joint_SNPs <- read.csv(file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_summary_AP.csv")

a_joint_snps_ordered <- read.csv(file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_ordered_AP.csv")
n_joint_snps_ordered <- read.csv(file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_ordered_AP.csv")

# aGVHD groups
a_Num_joint <- dim(a_joint_SNPs)[1]
p_values <- vector(mode = "numeric", length = a_Num_joint)
for(id in 1:a_Num_joint){
  
  rsNames <- unlist(strsplit(as.character(a_joint_SNPs$X[id]),"-"))
  # rs1_index <- which(grepl(rsNames[1], names(a_individual_SNP_freq)))
  # rs2_index <- which(grepl(rsNames[2], names(a_individual_SNP_freq)))
  # 
  # a_individual_SNP_freq[rs1_index]
  # a_individual_SNP_freq[rs2_index]
  # a_individual_SNP_freq[rs1_index] * a_individual_SNP_freq[rs2_index]
  
  ordered_index <- which(a_joint_snps_ordered$X == a_joint_SNPs$X[id])
  
  a_joint_SNPs$presence.HET..presence.HET.[id] #/dim(aSNP_mat)[1]
  
  chi_mat <- matrix(c(a_joint_SNPs$presence.HET..presence.HET.[id], a_joint_snps_ordered$presence.HET..absence[ordered_index],
                      a_joint_snps_ordered$absence.presence.HET.[ordered_index], a_joint_SNPs$absence.absence[id]), nrow = 2)
  
  # chisq.test(chi_mat)$statistic
  # chisq.test(chi_mat)$expected
  # 
  p_values[id] <- fisher.test(chi_mat)$"p.value"
  
  
}

names(p_values) <- a_joint_SNPs$X

a_ind_rsPairs <- names(p_values)[which(p_values <= 0.05)]
a_cor_rsPairs <- names(p_values)[which(p_values > 0.05)]

# nonGVHD groups

n_Num_joint <- dim(n_joint_SNPs)[1]
n_p_values <- vector(mode = "numeric", length = n_Num_joint)
for(id in 1:n_Num_joint){
  
  rsNames <- unlist(strsplit(as.character(n_joint_SNPs$X[id]),"-"))
  # rs1_index <- which(grepl(rsNames[1], names(a_individual_SNP_freq)))
  # rs2_index <- which(grepl(rsNames[2], names(a_individual_SNP_freq)))
  # 
  # a_individual_SNP_freq[rs1_index]
  # a_individual_SNP_freq[rs2_index]
  # a_individual_SNP_freq[rs1_index] * a_individual_SNP_freq[rs2_index]
  
  ordered_index <- which(n_joint_snps_ordered$X == n_joint_SNPs$X[id])
  
  n_joint_SNPs$presence.HET..presence.HET.[id] #/dim(aSNP_mat)[1]
  
  chi_mat <- matrix(c(n_joint_SNPs$presence.HET..presence.HET.[id], n_joint_snps_ordered$presence.HET..absence[ordered_index],
                      n_joint_snps_ordered$absence.presence.HET.[ordered_index], n_joint_SNPs$absence.absence[id]), nrow = 2)
  
  # chisq.test(chi_mat)$statistic
  # chisq.test(chi_mat)$expected
  # 
  n_p_values[id] <- fisher.test(chi_mat)$"p.value"
  
}

names(n_p_values) <- n_joint_SNPs$X

n_ind_rsPairs <- names(n_p_values)[which(n_p_values <= 0.05)]
n_cor_rsPairs <- names(n_p_values)[which(n_p_values > 0.05)]

################
## independent variables
a_ind_rs <- do.call("rbind", strsplit(a_ind_rsPairs, "-"))
n_ind_rs <- do.call("rbind", strsplit(n_ind_rsPairs, "-"))

table(a_ind_rs)

table(n_ind_rs)

### dependent variables
a_cor_rs <- do.call("rbind", strsplit(a_cor_rsPairs, "-"))  # 488
n_cor_rs <- do.call("rbind", strsplit(n_cor_rsPairs, "-"))  # 483

table(a_cor_rs)

table(n_cor_rs)


#### Outcome vs rs-rs association test
a_joint_SNPs <- read.csv(file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_summary_AP.csv")
n_joint_SNPs <- read.csv(file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_summary_AP.csv")

a_joint_snps_ordered <- read.csv(file = "../ClinVar/Data/unrestricted_aGVHD_pairwise_ordered_AP.csv")
n_joint_snps_ordered <- read.csv(file = "../ClinVar/Data/unrestricted_nonGVHD_pairwise_ordered_AP.csv")

Num_joint <- dim(a_joint_SNPs)[1]
p_values <- vector(mode = "numeric", length = Num_joint)
names(p_values) <- as.character(a_joint_SNPs$X)
for(id in 1:Num_joint){
  
  rsNames <- gsub("-", ".", as.character(a_joint_SNPs$X[id]))
  # rsNames <- unlist(strsplit(as.character(a_joint_SNPs$X[id]),"-"))
  # rs1_index <- which(grepl(rsNames[1], names(a_individual_SNP_freq)))
  # rs2_index <- which(grepl(rsNames[2], names(a_individual_SNP_freq)))
  # 
  # a_individual_SNP_freq[rs1_index]
  # a_individual_SNP_freq[rs2_index]
  # a_individual_SNP_freq[rs1_index] * a_individual_SNP_freq[rs2_index]
  
  ordered_index <- which(a_joint_snps_ordered$X == a_joint_SNPs$X[id])
  
  #cat(a_joint_SNPs$presence.HET..presence.HET.[id]) #/dim(aSNP_mat)[1]
  if(a_joint_SNPs$presence.HET..presence.HET.[id] >0 || n_joint_SNPs$presence.HET..presence.HET.[id]>0){
    cat(id, "\n")
    eval(parse(text = paste0("test_matrix <- matrix(c(a_joint_snps_ordered$absence.presence.HET.[ordered_index], 
                      n_joint_snps_ordered$absence.presence.HET.[ordered_index], 
                             a_joint_snps_ordered$presence.HET..absence[ordered_index],
                             n_joint_snps_ordered$presence.HET..absence[ordered_index],
                             a_joint_SNPs$presence.HET..presence.HET.[id],
                             n_joint_SNPs$presence.HET..presence.HET.[id]), 
                             nrow = 2, 
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"absence-presence\", \"presence-absence\", \"presence-presence\")
                             ))")))
  }else{
    eval(parse(text = paste0("test_matrix <- matrix(c(a_joint_snps_ordered$absence.presence.HET.[ordered_index], 
                      n_joint_snps_ordered$absence.presence.HET.[ordered_index], 
                             a_joint_snps_ordered$presence.HET..absence[ordered_index],
                             n_joint_snps_ordered$presence.HET..absence[ordered_index]), 
                             nrow = 2, 
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"absence-presence\", \"presence-absence\")
                             ))")))
  }
  
  
  # chisq.test(chi_mat)$statistic
  # chisq.test(chi_mat)$expected
  # 
  p_values[id] <- fisher.test(test_matrix)$"p.value"
  
  
  
}

p_values[order(p_values, decreasing = F)]

write.csv(p_values[order(p_values, decreasing = F)], row.names = T, file = "../ClinVar/Data/unrestricted_fisher_test_pvalues.csv")

which(p_values < 0.01)

# print to a file
sink(file = "../ClinVar/unrestricted_Fisher_test_pvalue.txt", append = T)
ordered_ID <- order(p_values, decreasing = F)
for(iid in 1:Num_joint){
  
  id <- ordered_ID[iid]
  rsNames <- gsub("-", ".", as.character(a_joint_SNPs$X[id]))
  # rsNames <- unlist(strsplit(as.character(a_joint_SNPs$X[id]),"-"))
  # rs1_index <- which(grepl(rsNames[1], names(a_individual_SNP_freq)))
  # rs2_index <- which(grepl(rsNames[2], names(a_individual_SNP_freq)))
  # 
  # a_individual_SNP_freq[rs1_index]
  # a_individual_SNP_freq[rs2_index]
  # a_individual_SNP_freq[rs1_index] * a_individual_SNP_freq[rs2_index]
  
  ordered_index <- which(a_joint_snps_ordered$X == a_joint_SNPs$X[id])
  
  #cat(a_joint_SNPs$presence.HET..presence.HET.[id]) #/dim(aSNP_mat)[1]
  if(a_joint_SNPs$presence.HET..presence.HET.[id] >0 || n_joint_SNPs$presence.HET..presence.HET.[id]>0){
    cat("ID = ", id, "\n")
    eval(parse(text = paste0("test_matrix <- matrix(c(a_joint_snps_ordered$absence.presence.HET.[ordered_index], 
                             n_joint_snps_ordered$absence.presence.HET.[ordered_index], 
                             a_joint_snps_ordered$presence.HET..absence[ordered_index],
                             n_joint_snps_ordered$presence.HET..absence[ordered_index],
                             a_joint_SNPs$presence.HET..presence.HET.[id],
                             n_joint_SNPs$presence.HET..presence.HET.[id]), 
                             nrow = 2, 
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"absence-presence\", \"presence-absence\", \"presence-presence\")
                             ))")))
  }else{
    cat("ID = ", id, "\n")
    eval(parse(text = paste0("test_matrix <- matrix(c(a_joint_snps_ordered$absence.presence.HET.[ordered_index], 
                             n_joint_snps_ordered$absence.presence.HET.[ordered_index], 
                             a_joint_snps_ordered$presence.HET..absence[ordered_index],
                             n_joint_snps_ordered$presence.HET..absence[ordered_index]), 
                             nrow = 2, 
                             dimnames = list(outcome = c(\"aGVHD\", \"non-GVHD\"),", rsNames, "
                             = c(\"absence-presence\", \"presence-absence\")
                             ))")))
  }
  
  print(test_matrix)
  cat("\n Fisher's exact test (two.sided) p-value: ")
  cat(p_values[id])
  cat("\n====================================\n")
  
  }
sink()
#### RandomForest variable importance
library(randomForest)

Feat <- rbind(aSNP_mat[, -c(1,2)], nSNP_mat[, -c(1,2)])
Labels <- c(rep("a", dim(aSNP_mat)[1]), 
            rep("n", dim(nSNP_mat)[1]))
Labels <- as.factor(Labels)
feat_lab <- cbind(Feat, Labels)

RF_feat <- randomForest(x = Feat, y = Labels, importance=TRUE)
# round(importance(RF_feat), 2)
imp <- importance(RF_feat, type =1, scale=T) 
importances_order <- order(imp, decreasing = T)
ordered_feat_set <- colnames(Feat)[importances_order]
ordered_feat_imp <- imp[importances_order]

### 
# combinations of SNPs
aa <- which(grepl("rs7958311", as.character(a_joint_snps_ordered$X)))
a_joint_snps_ordered[aa, c(1, 2, 3)]

nn <- which(grepl("rs7958311", as.character(n_joint_snps_ordered$X)))
n_joint_snps_ordered[nn, c(1, 2, 3)]

###############
