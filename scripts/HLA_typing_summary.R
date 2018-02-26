source("util.R")

library(ggplot2)
# library(xlsx)
# HLA_typings <- read.xlsx2(file = "../HLI_hla_mg_v3.xlsx", sheetIndex = 1)
HLA_typings <- read.csv(file = "../HLI_hla_mg_v3.csv")

load("../Data/ID_table.RData")

available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]

avaialble_IDs_HLA_typing <- HLA_typings[(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID),]
num_groups <- dim(avaialble_IDs_HLA_typing)[1]
group_type <- sapply(1:num_groups, function(x) available_IDs$Group[available_IDs$R_D_ID %in% avaialble_IDs_HLA_typing$nmdp_rid[x]])
group_type[group_type == "a"] <- "aGVHD"
group_type[group_type == "n"] <- "non-GVHD"
# num_groups <- dim(avaialble_IDs_HLA_typing)[1]

HLA_A <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_a_typ1), as.character(avaialble_IDs_HLA_typing$r_a_typ2)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_B <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_b_typ1), as.character(avaialble_IDs_HLA_typing$r_b_typ2)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_C <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_c_typ1), as.character(avaialble_IDs_HLA_typing$r_c_typ2)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_DRB1 <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_drb1_typ1), as.character(avaialble_IDs_HLA_typing$r_drb1_typ2)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_DQB1 <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_dqb1_typ1), as.character(avaialble_IDs_HLA_typing$r_dqb1_typ2)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_DPB1 <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_dpb1_typ1), as.character(avaialble_IDs_HLA_typing$r_dpb1_typ2)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)

# HLA_B <- c(as.character(avaialble_IDs_HLA_typing$r_b_typ1), as.character(avaialble_IDs_HLA_typing$r_b_typ2))
# HLA_C <- c(as.character(avaialble_IDs_HLA_typing$r_c_typ1), as.character(avaialble_IDs_HLA_typing$r_c_typ2))
# HLA_DRB1 <- c(as.character(avaialble_IDs_HLA_typing$r_drb1_typ1), as.character(avaialble_IDs_HLA_typing$r_drb1_typ2))
# HLA_DQB1 <- c(as.character(avaialble_IDs_HLA_typing$r_dqb1_typ1), as.character(avaialble_IDs_HLA_typing$r_dqb1_typ2))
# HLA_DPB1 <- c(as.character(avaialble_IDs_HLA_typing$r_dpb1_typ1), as.character(avaialble_IDs_HLA_typing$r_dpb1_typ2))


plot_typing_summary(HLA_A, "HLA-A Typing")
# plot_typing_summary(HLA_B, "HLA-B", "#D55E00")
# plot_typing_summary(HLA_C, "HLA-C", "#0072B2")
# plot_typing_summary(HLA_DRB1, "HLA-DRB1", "#CC79A7")
# plot_typing_summary(HLA_DQB1, "HLA-DQB1", "#56B4E9")
# plot_typing_summary(HLA_DPB1, "HLA-DPB1", "#E69B00")
plot_typing_summary(HLA_B, "HLA-B Typing")
plot_typing_summary(HLA_C, "HLA-C Typing")
plot_typing_summary(HLA_DRB1, "HLA-DRB1 Typing")
plot_typing_summary(HLA_DQB1, "HLA-DQB1 Typing")
plot_typing_summary(HLA_DPB1, "HLA-DPB1 Typing")


## Broad Race
avaialble_IDs_HLA_typing$rid.broad.race <- as.character(avaialble_IDs_HLA_typing$rid.broad.race)
avaialble_IDs_HLA_typing$did.broad.race <- as.character(avaialble_IDs_HLA_typing$did.broad.race)

Broad_pair <- sapply(1:num_groups, function(x) paste0(avaialble_IDs_HLA_typing[x, c("rid.broad.race","did.broad.race")], collapse = "-"))
Group_Broad_pair <- data.frame(Typing = Broad_pair,
                               GroupType = group_type,
                               stringsAsFactors = F)
plot_typing_summary (Group_Broad_pair, "Broad Race Pairs (Recipient - Donor)", fill_color = "#56B4E9")


### Sex
avaialble_IDs_HLA_typing$rid_sex <- as.character(avaialble_IDs_HLA_typing$rid_sex)
avaialble_IDs_HLA_typing$donor_sex <- as.character(avaialble_IDs_HLA_typing$donor_sex)

# Sex <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$rid_sex), as.character(avaialble_IDs_HLA_typing$donor_sex)),
#                   GroupType = rep(group_type, each = 2),
#                   stringsAsFactors = F)

Sex_pair <- sapply(1:num_groups, function(x) paste0(avaialble_IDs_HLA_typing[x, c("donor_sex", "rid_sex")], collapse = " -> "))
Group_Gender_pair <- data.frame(Typing = Sex_pair,
                               GroupType = group_type,
                               stringsAsFactors = F)
plot_typing_summary (Group_Gender_pair, "Sex (Recipient - Donor)", fill_color = "#56B4E9", Angle = 0, H = 0.5)

ggplot(Group_Gender_pair, aes(Typing, fill = GroupType)) + 
  geom_bar(position="dodge") +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) 

tab_sex_pair <- as.data.frame(table(Group_Gender_pair))
mat_sex_pair <- cbind(tab_sex_pair$Freq[1:4], tab_sex_pair$Freq[5:8])
colnames(mat_sex_pair) <- c("aGVHD", "non-aGVHD")
rownames(mat_sex_pair) <- as.character(tab_sex_pair$Typing)[1:4]

## overall 
chisq.test(mat_sex_pair, correct = T)

## sex-match vs outcome
chisq.test(mat_sex_pair[c(1,4),], correct = T)  # p-value = 0.2576
fisher.test(mat_sex_pair[c(1,4),]) # p-value = 0.3021

## sex-mismatch vs outcome
chisq.test(mat_sex_pair[c(2,3),], correct = T)  # p = 0.001587
fisher.test(mat_sex_pair[c(2,3),]) # p = 0.001059

## sex mismatch/match vs outcome
matched <- colSums(mat_sex_pair[c(1,4),]) 
mismatched <- colSums(mat_sex_pair[c(2,3),])

new_mat <- rbind(matched, mismatched)
rownames(new_mat) <- c("matched", "mismatched")

chisq.test(new_mat, correct = T) # p = 0.8398
fisher.test(new_mat) # p = 0.8883

## given recipient male, donor vs. outcome
chisq.test(mat_sex_pair[c(2,4),], correct = T)  # p = 0.03522
fisher.test(mat_sex_pair[c(2,4),]) # p = 0.02473

## given recipient female, donor vs. outcome
chisq.test(mat_sex_pair[c(1, 3),], correct = T)  # p = 1
fisher.test(mat_sex_pair[c(1,3),]) # p = 1

## recipient sex vs outcome
rec_male <- colSums(mat_sex_pair[c(2,4),]) 
rec_female <- colSums(mat_sex_pair[c(1,3),])

rec_sex <- rbind(rec_male, rec_female)
rownames(rec_sex) <- c("recMale", "recFemale")
chisq.test(rec_sex, correct = T)  # p = 0.00935
fisher.test(rec_sex) # p = 0.00771

### donor sex vs outcome
donor_male <- colSums(mat_sex_pair[c(3,4),]) 
donor_female <- colSums(mat_sex_pair[c(1,2),])

donor_sex <- rbind(donor_male, donor_female)
rownames(donor_sex) <- c("donorMale", "donorFemale")
chisq.test(donor_sex, correct = T)  # p = 0.2011
fisher.test(donor_sex) # p = 0.1673

#########################################
## summary of ambiguous typings
#########################################
source("util.R")

library(ggplot2)
# library(xlsx)
# HLA_typings <- read.xlsx2(file = "../HLI_hla_mg_v3.xlsx", sheetIndex = 1)
HLA_typings <- read.csv(file = "../HLI_hla_mg_v3.csv")

# avialable IDs
load("../Data/ID_table.RData")

R_aa <- sort(c(as.character(HLA_typings$r_a_typ1), as.character(HLA_typings$r_a_typ2)))
D_aa <- sort(c(as.character(HLA_typings$d_a_typ1), as.character(HLA_typings$d_a_typ2)))

num_pairs <- dim(HLA_typings)[1]
unmatched <- list()
HLA_typing <- data.frame(ID = numeric(dim(HLA_typings)[1]*6),
                         typ1 = character(dim(HLA_typings)[1]*6),
                         typ2 = character(dim(HLA_typings)[1]*6), 
                         Ambiguous = numeric(dim(HLA_typings)[1]*6),
                         Group = character(dim(HLA_typings)[1]*6),
                         stringsAsFactors = F)

available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]

available_ID_index <- which(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID)

for(id in 1:num_pairs){
  
  HLA_typing$ID[(6*(id-1)+1):(6*id)] <- id
  
  if(id %in% available_ID_index){
    groupType <- available_IDs$Group[available_IDs$R_D_ID %in% HLA_typings$nmdp_rid[id]]
    if(groupType == "n") HLA_typing$Group[(6*(id-1)+1):(6*id)] <- "non-aGVHD" else if(groupType == "a") HLA_typing$Group[(6*(id-1)+1):(6*id)] <- "aGVHD"
  }
  
  ### HLA-A
  HLA_A_d <- sort(c(as.character(HLA_typings$d_a_typ1[id]), as.character(HLA_typings$d_a_typ2[id])))
  HLA_A_r <- sort(c(as.character(HLA_typings$r_a_typ1[id]), as.character(HLA_typings$r_a_typ2[id])))
  if(identical(HLA_A_d, HLA_A_r)){ # matching
    HLA_typing[6*(id-1)+1, c(2,3)] <- HLA_A_d
    AA <- strsplit(HLA_A_d, "\\*")
    HLA_typing$Ambiguous[6*(id-1)+1] <-length(which(sapply(AA, function(x) grepl("[[:alpha:]]", x[2]))))
  }else { # not matching
    if(HLA_A_d[1] != HLA_A_r[1]){ # First allele, if not matching
      # remove false positive mismatches
      HLA_A_d_field <- unlist(strsplit(HLA_A_d[1], ":"))
      HLA_A_r_field <- unlist(strsplit(HLA_A_r[1], ":"))
      if(HLA_A_d_field[1] == HLA_A_r_field[1] & HLA_A_d_field[2] == HLA_A_r_field[2]){
        HLA_typing$typ1[6*(id-1)+1] <-  paste0(HLA_A_d[1], "~" ,HLA_A_r[1]) # false positive
      }else HLA_typing$typ1[6*(id-1)+1] <- paste0(HLA_A_d[1], "<>" ,HLA_A_r[1]) 
      
    }else{ # if the first allele is matching
      HLA_typing$typ1[6*(id-1)+1] <- HLA_A_d[1]
    }
    
    if(HLA_A_d[2] != HLA_A_r[2]){ # second allele
      # remove false positive mismatches
      HLA_A_d_field <- unlist(strsplit(HLA_A_d[2], ":"))
      HLA_A_r_field <- unlist(strsplit(HLA_A_r[2], ":"))
      if(HLA_A_d_field[1] == HLA_A_r_field[1] & HLA_A_d_field[2] == HLA_A_r_field[2]){
        HLA_typing$typ2[6*(id-1)+1] <-  paste0(HLA_A_d[2], "~" ,HLA_A_r[2]) # false positive
      }else HLA_typing$typ2[6*(id-1)+1] <- paste0(HLA_A_d[2], "<>" ,HLA_A_r[2])
  
    }else{# if the second allele is matching
      HLA_typing$typ2[6*(id-1)+1] <- HLA_A_d[2]
    }
    
    AA <- strsplit(unique(c(HLA_A_d, HLA_A_r)), "\\*")
    HLA_typing$Ambiguous[6*(id-1)+1] <-length(which(sapply(AA, function(x) grepl("[[:alpha:]]", x[2]))))
  }
  
  ### HLA-B
  HLA_B_d <- sort(c(as.character(HLA_typings$d_b_typ1[id]), as.character(HLA_typings$d_b_typ2[id])))
  HLA_B_r <- sort(c(as.character(HLA_typings$r_b_typ1[id]), as.character(HLA_typings$r_b_typ2[id])))
  if(identical(HLA_B_d, HLA_B_r)){
    HLA_typing[6*(id-1)+2, c(2,3)] <- HLA_B_d
    BB <- strsplit(HLA_B_d, "\\*")
    HLA_typing$Ambiguous[6*(id-1)+2] <-length(which(sapply(BB, function(x) grepl("[[:alpha:]]", x[2]))))
  }else{
    # HLA_typing[6*(id-1)+2, c(2,3)] <- "No Matching"
    if(HLA_B_d[1] != HLA_B_r[1]){
      # remove false positive mismatches
      HLA_B_d_field <- unlist(strsplit(HLA_B_d[1], ":"))
      HLA_B_r_field <- unlist(strsplit(HLA_B_r[1], ":"))
      if(HLA_B_d_field[1] == HLA_B_r_field[1] & HLA_B_d_field[2] == HLA_B_r_field[2]){
        HLA_typing$typ1[6*(id-1)+2] <-  paste0(HLA_B_d[1], "~" ,HLA_B_r[1]) # false positive
      }else HLA_typing$typ1[6*(id-1)+2] <- paste0(HLA_B_d[1], "<>" ,HLA_B_r[1])
      
    }else{
      HLA_typing$typ1[6*(id-1)+2] <- HLA_B_d[1]
    }
    if(HLA_B_d[2] != HLA_B_r[2]){ # second allele
      # remove false positive mismatches
      HLA_B_d_field <- unlist(strsplit(HLA_B_d[2], ":"))
      HLA_B_r_field <- unlist(strsplit(HLA_B_r[2], ":"))
      if(HLA_B_d_field[1] == HLA_B_r_field[1] & HLA_B_d_field[2] == HLA_B_r_field[2]){
        HLA_typing$typ2[6*(id-1)+2] <-  paste0(HLA_B_d[2], "~" ,HLA_B_r[2]) # false positive
      }else HLA_typing$typ2[6*(id-1)+2] <- paste0(HLA_B_d[2], "<>" ,HLA_B_r[2])
      
    }else{
      HLA_typing$typ2[6*(id-1)+2] <- HLA_B_d[2]
    }

    BB <- strsplit(unique(c(HLA_B_d, HLA_B_r)), "\\*")
    HLA_typing$Ambiguous[6*(id-1)+2] <-length(which(sapply(BB, function(x) grepl("[[:alpha:]]", x[2]))))
  }
  
  ## HLA-C
  HLA_C_d <- sort(c(as.character(HLA_typings$d_c_typ1[id]), as.character(HLA_typings$d_c_typ2[id])))
  HLA_C_r <- sort(c(as.character(HLA_typings$r_c_typ1[id]), as.character(HLA_typings$r_c_typ2[id])))
  if(identical(HLA_C_d, HLA_C_r)){
    HLA_typing[6*(id-1)+3, c(2,3)] <- HLA_C_d
    CC <- strsplit(HLA_C_d, "\\*")
    HLA_typing$Ambiguous[6*(id-1)+3] <- length(which(sapply(CC, function(x) grepl("[[:alpha:]]", x[2]))))
  }else{
    # HLA_typing[6*(id-1)+3, c(2,3)] <- "No Matching"
    if(HLA_C_d[1] != HLA_C_r[1]){
      # remove false positive mismatches
      HLA_C_d_field <- unlist(strsplit(HLA_C_d[1], ":"))
      HLA_C_r_field <- unlist(strsplit(HLA_C_r[1], ":"))
      if(HLA_C_d_field[1] == HLA_C_r_field[1] & HLA_C_d_field[2] == HLA_C_r_field[2]){
        HLA_typing$typ1[6*(id-1)+3] <-  paste0(HLA_C_d[1], "~" ,HLA_C_r[1]) # false positive
      }else HLA_typing$typ1[6*(id-1)+3] <- paste0(HLA_C_d[1], "<>" ,HLA_C_r[1])
      
      
    }else{
      HLA_typing$typ1[6*(id-1)+3] <- HLA_C_d[1]
    }
    
    if(HLA_C_d[2] != HLA_C_r[2]){
      # remove false positive mismatches
      HLA_C_d_field <- unlist(strsplit(HLA_C_d[2], ":"))
      HLA_C_r_field <- unlist(strsplit(HLA_C_r[2], ":"))
      if(HLA_C_d_field[1] == HLA_C_r_field[1] & HLA_C_d_field[2] == HLA_C_r_field[2]){
        HLA_typing$typ2[6*(id-1)+3] <-  paste0(HLA_C_d[2], "~" ,HLA_C_r[2]) # false positive
      }else HLA_typing$typ2[6*(id-1)+3] <- paste0(HLA_C_d[2], "<>" ,HLA_C_r[2])
      
     
    }else{
      HLA_typing$typ2[6*(id-1)+3] <- HLA_C_d[2]
    }
    
    
    CC <- strsplit(unique(c(HLA_C_d, HLA_C_r)), "\\*")
    HLA_typing$Ambiguous[6*(id-1)+3] <-length(which(sapply(CC, function(x) grepl("[[:alpha:]]", x[2]))))
  }
  
  ##### HLA-DRB1
  HLA_DRB1_d <- sort(c(as.character(HLA_typings$d_drb1_typ1[id]), as.character(HLA_typings$d_drb1_typ2[id])))
  HLA_DRB1_r <- sort(c(as.character(HLA_typings$r_drb1_typ1[id]), as.character(HLA_typings$r_drb1_typ2[id])))
  if(identical(HLA_DRB1_d, HLA_DRB1_r)){
    HLA_typing[6*(id-1)+4, c(2,3)] <- HLA_DRB1_d
    DRBB <- strsplit(HLA_DRB1_d, "\\*")
    HLA_typing$Ambiguous[6*(id-1)+4] <- length(which(sapply(DRBB, function(x) grepl("[[:alpha:]]", x[2]))))
  }else{
    # HLA_typing[6*(id-1)+4, c(2,3)] <- "No Matching"
    
    if(HLA_DRB1_d[1] != HLA_DRB1_r[1]){
      # remove false positive mismatches
      HLA_DRB1_d_field <- unlist(strsplit(HLA_DRB1_d[1], ":"))
      HLA_DRB1_r_field <- unlist(strsplit(HLA_DRB1_r[1], ":"))
      if(HLA_DRB1_d_field[1] == HLA_DRB1_r_field[1] & HLA_DRB1_d_field[2] == HLA_DRB1_r_field[2]){
        HLA_typing$typ1[6*(id-1)+4] <-  paste0(HLA_DRB1_d[1], "~" ,HLA_DRB1_r[1]) # false positive
      }else HLA_typing$typ1[6*(id-1)+4] <- paste0(HLA_DRB1_d[1], "<>" ,HLA_DRB1_r[1])
      
    }else{
      HLA_typing$typ1[6*(id-1)+4] <- HLA_DRB1_d[1]
    }
    if(HLA_DRB1_d[2] != HLA_DRB1_r[2]){
      # remove false positive mismatches
      HLA_DRB1_d_field <- unlist(strsplit(HLA_DRB1_d[2], ":"))
      HLA_DRB1_r_field <- unlist(strsplit(HLA_DRB1_r[2], ":"))
      if(HLA_DRB1_d_field[1] == HLA_DRB1_r_field[1] & HLA_DRB1_d_field[2] == HLA_DRB1_r_field[2]){
        HLA_typing$typ2[6*(id-1)+4] <-  paste0(HLA_DRB1_d[2], "~" ,HLA_DRB1_r[2]) # false positive
      }else HLA_typing$typ2[6*(id-1)+4] <- paste0(HLA_DRB1_d[2], "<>" ,HLA_DRB1_r[2])
      
      
    }else{
      HLA_typing$typ2[6*(id-1)+4] <- HLA_DRB1_d[2]
    }
    
    
    DRBB <- strsplit(unique(c(HLA_DRB1_d, HLA_DRB1_r)), "\\*")
    HLA_typing$Ambiguous[6*(id-1)+4] <-length(which(sapply(DRBB, function(x) grepl("[[:alpha:]]", x[2]))))
  }
  
  ####### HLA-DQB1
  HLA_DQB1_d <- sort(c(as.character(HLA_typings$d_dqb1_typ1[id]), as.character(HLA_typings$d_dqb1_typ2[id])))
  HLA_DQB1_r <- sort(c(as.character(HLA_typings$r_dqb1_typ1[id]), as.character(HLA_typings$r_dqb1_typ2[id])))
  if(identical(HLA_DQB1_d, HLA_DQB1_r)){
    HLA_typing[6*(id-1)+5, c(2,3)] <- HLA_DQB1_d
    DQBB <- strsplit(HLA_DQB1_d, "\\*")
    HLA_typing$Ambiguous[6*(id-1)+5] <- length(which(sapply(DQBB, function(x) grepl("[[:alpha:]]", x[2]))))
  }else{
    # HLA_typing[6*(id-1)+5, c(2,3)] <- "No Matching"
    
    if(HLA_DQB1_d[1] != HLA_DQB1_r[1]){ # first allele
      # remove false positive mismatches
      HLA_DQB1_d_field <- unlist(strsplit(HLA_DQB1_d[1], ":"))
      HLA_DQB1_r_field <- unlist(strsplit(HLA_DQB1_r[1], ":"))
      if(HLA_DQB1_d_field[1] == HLA_DQB1_r_field[1] & HLA_DQB1_d_field[2] == HLA_DQB1_r_field[2]){
        HLA_typing$typ1[6*(id-1)+5] <-  paste0(HLA_DQB1_d[1], "~" ,HLA_DQB1_r[1]) # false positive
      }else HLA_typing$typ1[6*(id-1)+5] <- paste0(HLA_DQB1_d[1], "<>" ,HLA_DQB1_r[1])
      
      
    }else{
      HLA_typing$typ1[6*(id-1)+5] <- HLA_DQB1_d[1]
    }
    if(HLA_DQB1_d[2] != HLA_DQB1_r[2]){
      # remove false positive mismatches
      HLA_DQB1_d_field <- unlist(strsplit(HLA_DQB1_d[2], ":"))
      HLA_DQB1_r_field <- unlist(strsplit(HLA_DQB1_r[2], ":"))
      if(HLA_DQB1_d_field[1] == HLA_DQB1_r_field[1] & HLA_DQB1_d_field[2] == HLA_DQB1_r_field[2]){
        HLA_typing$typ2[6*(id-1)+5] <-  paste0(HLA_DQB1_d[2], "~" ,HLA_DQB1_r[2]) # false positive
      }else HLA_typing$typ2[6*(id-1)+5] <- paste0(HLA_DQB1_d[2], "<>" ,HLA_DQB1_r[2])
      
      
    }else{
      HLA_typing$typ2[6*(id-1)+5] <- HLA_DQB1_d[2]
    }
    
    DQBB <- strsplit(unique(c(HLA_DQB1_d, HLA_DQB1_r)), "\\*")
    HLA_typing$Ambiguous[6*(id-1)+5] <-length(which(sapply(DQBB, function(x) grepl("[[:alpha:]]", x[2]))))
  }
  
  ##### HLA-DPB1
  HLA_DPB1_d <- sort(c(as.character(HLA_typings$d_dpb1_typ1[id]), as.character(HLA_typings$d_dpb1_typ2[id])))
  HLA_DPB1_r <- sort(c(as.character(HLA_typings$r_dpb1_typ1[id]), as.character(HLA_typings$r_dpb1_typ2[id])))
  if(identical(HLA_DPB1_d, HLA_DPB1_r)){
    HLA_typing[6*id, c(2,3)] <- HLA_DPB1_d
    DPBB <- strsplit(HLA_DPB1_d, "\\*")
    HLA_typing$Ambiguous[6*id] <- length(which(sapply(DPBB, function(x) grepl("[[:alpha:]]", x[2]))))
  }else{
    # HLA_typing[6*id, c(2,3)] <- "No Matching"
    
    if(HLA_DPB1_d[1] != HLA_DPB1_r[1]){
      # remove false positive mismatches
      HLA_DPB1_d_field <- unlist(strsplit(HLA_DPB1_d[1], ":"))
      HLA_DPB1_r_field <- unlist(strsplit(HLA_DPB1_r[1], ":"))
      
      if(length(HLA_DPB1_r_field) >0 & length(HLA_DPB1_d_field) >0 
         & HLA_DPB1_d_field[1] == HLA_DPB1_r_field[1] & HLA_DPB1_d_field[2] == HLA_DPB1_r_field[2]){
        HLA_typing$typ1[6*id] <-  paste0(HLA_DPB1_d[1], "~" ,HLA_DPB1_r[1]) # false positive
      }else HLA_typing$typ1[6*id] <- paste0(HLA_DPB1_d[1], "<>" ,HLA_DPB1_r[1])
      
     
    }else{
      HLA_typing$typ1[6*id] <- HLA_DPB1_d[1]
    }
    if(HLA_DPB1_d[2] != HLA_DPB1_r[2]){
      # remove false positive mismatches
      HLA_DPB1_d_field <- unlist(strsplit(HLA_DPB1_d[2], ":"))
      HLA_DPB1_r_field <- unlist(strsplit(HLA_DPB1_r[2], ":"))
      if(length(HLA_DPB1_r_field) >0 & length(HLA_DPB1_d_field) >0 
         & HLA_DPB1_d_field[1] == HLA_DPB1_r_field[1] & HLA_DPB1_d_field[2] == HLA_DPB1_r_field[2]){
        HLA_typing$typ2[6*id] <-  paste0(HLA_DPB1_d[2], "~" ,HLA_DPB1_r[2]) # false positive
      }else HLA_typing$typ2[6*id] <- paste0(HLA_DPB1_d[2], "<>" ,HLA_DPB1_r[2])
      
      
    }else{
      HLA_typing$typ2[6*id] <- HLA_DPB1_d[2]
    }
    
    DPBB <- strsplit(unique(c(HLA_DPB1_d, HLA_DPB1_r)), "\\*")
    HLA_typing$Ambiguous[6*id] <-length(which(sapply(DPBB, function(x) grepl("[[:alpha:]]", x[2]))))
  }
}

HLA_typing_stats_ten <- HLA_typing[-6*(1:num_pairs), ]
sum(HLA_typing_stats_ten$Ambiguous) ### 676 typings are ambiguous out of 1255 ~  0.5386454 
sum(HLA_typing$Ambiguous)  ### 1039 typings are ambiguous out of 1506 ~ 0.689907  (6 loci)

length(which(HLA_typing$Group!=""))/6  # 205 groups have WGS data

available_HLA_typing <- HLA_typing[HLA_typing$Group!="", ]
available_HLA_typing_ten <- available_HLA_typing[-6*(1:205), ]
sum(available_HLA_typing$Ambiguous)  ### 912 typings are ambiguous, out of 1230 ~ 0.7414634 (6 loci)
sum(available_HLA_typing_ten$Ambiguous)  ### 606 typings are ambiguous, out 1025 ~ 0.5912195

available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]

availble_IDs_HLA_typing <- HLA_typings[(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID),]

Ambiguous_groups <- which(sapply(1:205, function(x) sum(available_HLA_typing_ten$Ambiguous[(5*(x-1)+1):(5*x)]))==0)
##  70 groups out of 205 groups ambiguous typings 

table(available_HLA_typing_ten$Group[Ambiguous_groups])
# #    aGVHD non-aGVHD 
# #      28        42 

# No matching
typ1_mismatch <- which( grepl("<>",available_HLA_typing_ten$typ1))
typ2_mismatch <- which( grepl("<>",available_HLA_typing_ten$typ2))
length(unique(c(typ1_mismatch, typ2_mismatch)))  ## 85
table(available_HLA_typing_ten$Group[unique(c(typ1_mismatch, typ2_mismatch))])
# aGVHD non-aGVHD 
# 40        45 
# remove false positive

## different resolution
typ1_res <- which( grepl("~",available_HLA_typing_ten$typ1))
typ2_res <- which( grepl("~",available_HLA_typing_ten$typ2))
length(unique(c(typ1_res, typ2_res)))  ## 57
table(available_HLA_typing_ten$Group[unique(c(typ1_res, typ2_res))])
# aGVHD non-aGVHD 
# 28        29 
# remove false positive

length(unique(c(typ1_mismatch, typ2_mismatch, typ1_res, typ2_res)))

#########################################################################
HLA_typings <- read.csv(file = "../HLI_hla_mg_v3.csv")
load("../Data/ID_table.RData")

available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]

available_IDs_HLA_typing <- HLA_typings[(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID),]
num_groups <- dim(available_IDs_HLA_typing)[1]
group_type <- sapply(1:num_groups, function(x) available_IDs$Group[available_IDs$R_D_ID %in% available_IDs_HLA_typing$nmdp_rid[x]])
group_type[group_type == "a"] <- "aGVHD"
group_type[group_type == "n"] <- "non-GVHD"

reformated_HLA_typing_list <- available_IDs_HLA_typing[, 1:3]
reformated_HLA_typing_list$RRace <- available_IDs_HLA_typing$rid.broad.race
reformated_HLA_typing_list$DRace <- available_IDs_HLA_typing$did.broad.race
reformated_HLA_typing_list$RSex <- available_IDs_HLA_typing$rid_sex
reformated_HLA_typing_list$DSex <- available_IDs_HLA_typing$donor_sex

reformated_HLA_typing_list$HLA_A1 <- sapply(as.character(available_IDs_HLA_typing$r_a_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_A2 <- sapply(as.character(available_IDs_HLA_typing$r_a_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_B1 <- sapply(as.character(available_IDs_HLA_typing$r_b_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_B2 <- sapply(as.character(available_IDs_HLA_typing$r_b_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_C1 <- sapply(as.character(available_IDs_HLA_typing$r_c_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_C2 <- sapply(as.character(available_IDs_HLA_typing$r_c_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_DRB11 <- sapply(as.character(available_IDs_HLA_typing$r_drb1_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_DRB12 <- sapply(as.character(available_IDs_HLA_typing$r_drb1_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_DQB11 <- sapply(as.character(available_IDs_HLA_typing$r_dqb1_typ1.gl), function(x) unlist(strsplit(x, "/"))[1])
reformated_HLA_typing_list$HLA_DQB12 <- sapply(as.character(available_IDs_HLA_typing$r_dqb1_typ2.gl), function(x) unlist(strsplit(x, "/"))[1])

get_two_field <- function(HLA_typing){
  
  reformat_typing <- unlist(strsplit(HLA_typing, ":"))
  if(length(reformat_typing) > 2){
    
    two_fields <- paste0(reformat_typing[c(1,2)], collapse = ":")
    
  }else two_fields <- HLA_typing
  
  return(two_fields)
}

reformated_HLA_typing_list$HLA_A1 <- sapply(reformated_HLA_typing_list$HLA_A1, get_two_field)
reformated_HLA_typing_list$HLA_A2 <- sapply(reformated_HLA_typing_list$HLA_A2, get_two_field)
reformated_HLA_typing_list$HLA_B1 <- sapply(reformated_HLA_typing_list$HLA_B1, get_two_field)
reformated_HLA_typing_list$HLA_B2 <- sapply(reformated_HLA_typing_list$HLA_B2, get_two_field)
reformated_HLA_typing_list$HLA_C1 <- sapply(reformated_HLA_typing_list$HLA_C1, get_two_field)
reformated_HLA_typing_list$HLA_C2 <- sapply(reformated_HLA_typing_list$HLA_C2, get_two_field)
reformated_HLA_typing_list$HLA_DRB11 <- sapply(reformated_HLA_typing_list$HLA_DRB11, get_two_field)
reformated_HLA_typing_list$HLA_DRB12 <- sapply(reformated_HLA_typing_list$HLA_DRB12, get_two_field)
reformated_HLA_typing_list$HLA_DQB11 <- sapply(reformated_HLA_typing_list$HLA_DQB11, get_two_field)
reformated_HLA_typing_list$HLA_DQB12 <- sapply(reformated_HLA_typing_list$HLA_DQB12, get_two_field)

write.csv(reformated_HLA_typing_list, file = "../ClinVar/reformated_HLA_typing.cvs")

HLA_A <- data.frame(Typing = c(reformated_HLA_typing_list$HLA_A1, reformated_HLA_typing_list$HLA_A2),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_A_redundant_ID <- seq(from = 1, to = dim(HLA_A)[1], by = 2)[sapply(seq(from = 1, to = dim(HLA_A)[1], by = 2), function(x) if(HLA_A$Typing[x] == HLA_A$Typing[x+1]) TRUE else FALSE)]
HLA_A_non_redundant <- HLA_A[-HLA_A_redundant_ID, ]

HLA_B <- data.frame(Typing = c(reformated_HLA_typing_list$HLA_B1, reformated_HLA_typing_list$HLA_B2),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_B_redundant_ID <- seq(from = 1, to = dim(HLA_B)[1], by = 2)[sapply(seq(from = 1, to = dim(HLA_B)[1], by = 2), function(x) if(HLA_B$Typing[x] == HLA_B$Typing[x+1]) TRUE else FALSE)]
HLA_B_non_redundant <- HLA_B[-HLA_B_redundant_ID, ]

HLA_C <- data.frame(Typing = c(reformated_HLA_typing_list$HLA_C1, reformated_HLA_typing_list$HLA_C2),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_C_redundant_ID <- seq(from = 1, to = dim(HLA_C)[1], by = 2)[sapply(seq(from = 1, to = dim(HLA_C)[1], by = 2), function(x) if(HLA_C$Typing[x] == HLA_C$Typing[x+1]) TRUE else FALSE)]
HLA_C_non_redundant <- HLA_C[-HLA_C_redundant_ID, ]

HLA_DRB1 <- data.frame(Typing = c(reformated_HLA_typing_list$HLA_DRB11, reformated_HLA_typing_list$HLA_DRB12),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_DRB1_redundant_ID <- seq(from = 1, to = dim(HLA_DRB1)[1], by = 2)[sapply(seq(from = 1, to = dim(HLA_DRB1)[1], by = 2), function(x) if(HLA_DRB1$Typing[x] == HLA_DRB1$Typing[x+1]) TRUE else FALSE)]
HLA_DRB1_non_redundant <- HLA_DRB1[-HLA_DRB1_redundant_ID, ]

HLA_DQB1 <- data.frame(Typing = c(reformated_HLA_typing_list$HLA_DQB11, reformated_HLA_typing_list$HLA_DQB12),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_DQB1_redundant_ID <- seq(from = 1, to = dim(HLA_DQB1)[1], by = 2)[sapply(seq(from = 1, to = dim(HLA_DQB1)[1], by = 2), function(x) if(HLA_DQB1$Typing[x] == HLA_DQB1$Typing[x+1]) TRUE else FALSE)]
HLA_DQB1_non_redundant <- HLA_DQB1[-HLA_DQB1_redundant_ID, ]

plot_typing_summary(HLA_A_non_redundant, "HLA-A Typing")
# plot_typing_summary(HLA_A, "HLA-A Typing")

plot_typing_summary(HLA_B_non_redundant, "HLA-B Typing")
plot_typing_summary(HLA_C_non_redundant, "HLA-C Typing")
plot_typing_summary(HLA_DRB1_non_redundant, "HLA-DRB1 Typing")
plot_typing_summary(HLA_DQB1_non_redundant, "HLA-DQB1 Typing")


unique_HLA_A <- sapply(cbind(reformated_HLA_typing_list$HLA_A1, reformated_HLA_typing_list$HLA_A2), function(x) unique(x))






A_0201_ID <- sort(union(which(highest_IDs_table$reformated_HLA_typing_list %in% "02:01"), which(reformated_HLA_typing_list$HLA_A2 %in% "02:01")))
table(highest_IDs_table$HLA_B1[A_0201_ID])

A_0201_B_0702_ID2 <- sort(union(which(highest_IDs_table$HLA_B1[A_0201_ID] %in% "07:02"), which(highest_IDs_table$HLA_B2[A_0201_ID] %in% "07:02")))
A_0201_B_0702_ID <- A_0201_ID[A_0201_B_0702_ID2]

table(highest_IDs_table$HLA_C1[A_0201_B_0702_ID])
