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

