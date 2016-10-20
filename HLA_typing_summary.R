source("util.R")

library(ggplot2)
library(xlsx)
HLA_typings <- read.xlsx2(file = "../HLI_hla_mg_v3.xlsx", sheetIndex = 1)

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

Sex_pair <- sapply(1:num_groups, function(x) paste0(avaialble_IDs_HLA_typing[x, c("rid_sex","donor_sex")], collapse = "-"))
Group_Gender_pair <- data.frame(Typing = Sex_pair,
                               GroupType = group_type,
                               stringsAsFactors = F)
plot_typing_summary (Group_Gender_pair, "Sex (Recipient - Donor)", fill_color = "#56B4E9", Angle = 0, H = 0.5)

ggplot(Group_Gender_pair, aes(Typing, fill = GroupType)) + 
  geom_bar(position="dodge") +
  scale_fill_manual(values = c("#D55E00", "#0072B2")) + 
  scale_x_discrete(labels = c("F -> F", "M -> F", "F -> M", "M -> M"))
