DPB1_table <- read.csv("../ClinVar/DPB1-MiHA/MiHA-DPB1.01-29.csv")

IDs <- unique(DPB1_table$PAIR_ID)
num_IDs <- length(IDs)

highest_IDs_table <- data.frame(pairID = numeric(num_IDs),
                                RID = numeric(num_IDs),
                                DID = numeric(num_IDs),
                                RRACE.DRACE = character(num_IDs),
                                RSEX.DSEX = character(num_IDs),
                                OUTCOME = character(num_IDs),
                                RDP_MIHA = character(num_IDs),
                                DDP_MIHA = character(num_IDs),
                                # Match.of.10 = character(num_IDs),
                                EXPRESSION = character(num_IDs),
                                TCE = character(num_IDs),
                                FD = character(num_IDs),
                                DP_MATCH = character(num_IDs),
                                X04_MATCH = character(num_IDs),
                                PROB = character(num_IDs),
                                stringsAsFactors = F)

fid <- sapply(DPB1_table, is.factor)
DPB1_table2 <- DPB1_table
DPB1_table2[fid] <- lapply(DPB1_table2[fid], as.character)
for(id in 1:num_IDs){
  
  index <- which(DPB1_table2$PAIR_ID %in% IDs[id])
  max_id <- which.max(DPB1_table2$PROB[index])
  
  highest_IDs_table[id, ] <- DPB1_table2[index[max_id], ]
  
}

table(highest_IDs_table$TCE, highest_IDs_table$OUTCOME, highest_IDs_table$RDP_MIHA)
highest_IDs_table$SexMatch <- gsub(":", " <- ", highest_IDs_table$RSEX.DSEX)

##########
# barplot
##########
library(ggplot2)

#sexmismatch vs outcome vs RDP_MiHA
tb1 <- as.data.frame(table(highest_IDs_table$SexMatch, highest_IDs_table$OUTCOME, highest_IDs_table$RDP_MIHA))
colnames(tb1) <- c("SexMatch", "Outcome", "RDP_MIHA", "count")

# ggplot(tb1, aes(x = factor(RDP_MIHA), y = count, fill= Outcome , color= Outcome)) + 
#   geom_bar(stat = "identity", position=position_dodge()) +
#   geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
#   # geom_text(aes(y=count, ymax=count, label=count), position= position_dodge(width=0.9), vjust=-.5, color="black") +
#   scale_y_continuous("Count of cases",limits=c(0, 35), breaks=seq(0, 35, 5)) + 
#   scale_x_discrete("MiHA peptide in Recipient") # +
  # scale_fill_discrete(name ="Outcome", labels=c("Acute GvHD", "Non-GvHD"))

g1 <- ggplot(tb1, aes(x = factor(RDP_MIHA), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ SexMatch) +
  scale_y_continuous("Count of cases",limits=c(0, 35), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("MiHA peptide in Recipient") + 
  labs(title = "Fig.1(a) - MiHA peptide in Recipient") 

## 
#sexmismatch vs outcome vs RDP_MiHA
tb2 <- as.data.frame(table(highest_IDs_table$SexMatch, highest_IDs_table$OUTCOME, highest_IDs_table$EXPRESSION))
colnames(tb2) <- c("SexMatch", "Outcome", "EXPRESSION", "count")
# ggplot(tb2, aes(x = factor(EXPRESSION), y = count, fill= Outcome , color= Outcome)) + 
#   geom_bar(stat = "identity", position=position_dodge()) +
#   geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
#   # geom_text(aes(y=count, ymax=count, label=count), position= position_dodge(width=0.9), vjust=-.5, color="black") +
#   scale_y_continuous("Count of cases",limits=c(0, 35), breaks=seq(0, 35, 5)) + 
#   scale_x_discrete("Expression")
g2 <- ggplot(tb2, aes(x = factor(EXPRESSION), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ SexMatch) +
  scale_y_continuous("Count of cases",limits=c(0, 35), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("Expression") +
  labs(title = "Fig.2 - DPB1 Expression Risk") 

# ft.tb2 <- as.data.frame(table(highest_IDs_table$EXPRESSION, highest_IDs_table$OUTCOME, highest_IDs_table$RDP_MIHA))
# fisher.test(ft.tb2)

### TCE - permissibility
tb3 <- as.data.frame(table(highest_IDs_table$SexMatch, highest_IDs_table$OUTCOME, highest_IDs_table$TCE))
colnames(tb3) <- c("SexMatch", "Outcome", "TCE", "count")
g3 <- ggplot(tb3, aes(x = factor(TCE), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ SexMatch) +
  theme(axis.text.x=element_text(angle= 0, hjust= 0.5)) + 
  labs(title = "Fig.3 - T-Cell Epitope permissibility") +
  scale_y_continuous("Count of cases",limits=c(0, 30), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("TCE DPB1 Permissibility")

## FD - 
tb4 <- as.data.frame(table(highest_IDs_table$SexMatch, highest_IDs_table$OUTCOME, highest_IDs_table$FD))
colnames(tb4) <- c("SexMatch", "Outcome", "FD", "count")
g4 <- ggplot(tb4, aes(x = factor(FD), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ SexMatch) +
  # theme(axis.text.x=element_text(angle= 15, hjust= 0.5)) + 
  labs(title = "Fig.4 - Functional Distance permissibility") +
  scale_y_continuous("Count of cases",limits=c(0, 35), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("FD")

#### DP_MATCH
tb5 <- as.data.frame(table(highest_IDs_table$SexMatch, highest_IDs_table$OUTCOME, highest_IDs_table$DP_MATCH))
colnames(tb5) <- c("SexMatch", "Outcome", "DP_MATCH", "count")
g5 <- ggplot(tb5, aes(x = factor(DP_MATCH), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ SexMatch) +
  theme(axis.text.x=element_text(angle= 0, hjust= .5, vjust = .5)) +
  labs(title = "Fig.5 - HLA-DPB1 Match-mismatch") +
  scale_y_continuous("Count of cases",limits=c(0, 20), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("DPB1 Matching")


##### X04 MActh
tb6 <- as.data.frame(table(highest_IDs_table$SexMatch, highest_IDs_table$OUTCOME, highest_IDs_table$X04_MATCH))
colnames(tb6) <- c("SexMatch", "Outcome", "X04_MATCH", "count")
g6 <- ggplot(tb6, aes(x = factor(X04_MATCH), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ SexMatch) +
  theme(axis.text.x=element_text(angle= 90, hjust= 1, vjust = .5)) +
  labs(title = "Fig.6 - HLA-DPB1*04 Match") +
  scale_y_continuous("Count of cases",limits=c(0, 10), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("HLA-DPB1*04")

#sexmismatch vs outcome vs RDP_MiHA
tb7 <- as.data.frame(table(highest_IDs_table$SexMatch, highest_IDs_table$OUTCOME, highest_IDs_table$DDP_MIHA))
colnames(tb7) <- c("SexMatch", "Outcome", "DDP_MIHA", "count")

g7 <- ggplot(tb7, aes(x = factor(DDP_MIHA), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ SexMatch) +
  scale_y_continuous("Count of cases",limits=c(0, 35), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("MiHA peptide in Donor") + 
  labs(title = "Fig.1(b) - MiHA SNP allelic variable in Donor") 

#DPB1*04 vs outcome vs RDP_MiHA
tb8 <- as.data.frame(table(highest_IDs_table$DP_MATCH, highest_IDs_table$OUTCOME, highest_IDs_table$RDP_MIHA))
colnames(tb8) <- c("DP_MATCH", "Outcome", "RDP_MIHA", "count")

g8 <- ggplot(tb8, aes(x = factor(RDP_MIHA), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ DP_MATCH) +
  scale_y_continuous("Count of cases",limits=c(0, 35), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("MiHA peptide in Recipient") + 
  labs(title = "Fig.7 - MiHA in Recipient vs DPB1 match") 

#DPB1*04 vs outcome vs RDP_MiHA
tb9 <- as.data.frame(table(highest_IDs_table$X04_MATCH, highest_IDs_table$OUTCOME, highest_IDs_table$RDP_MIHA))
colnames(tb9) <- c("X04_MATCH", "Outcome", "RDP_MIHA", "count")

g9 <- ggplot(tb9, aes(x = factor(RDP_MIHA), y = count, fill= Outcome , color= Outcome)) +
  geom_bar(stat = "identity", position=position_dodge()) + 
  geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
  facet_grid(. ~ X04_MATCH) +
  scale_y_continuous("Count of cases",limits=c(0, 35), breaks=seq(0, 35, 5)) + 
  scale_x_discrete("MiHA peptide in Recipient") + 
  labs(title = "Fig.8 - MiHA in Recipient vs DPB1*04 match") 


######## DPB1*04 vs outcome vs RDP_MiHA vs Sex-mismatch
library(gridExtra)
SexMatch <- unique(highest_IDs_table$SexMatch)
sexMatch_types_num <- length(SexMatch)
gg10 <- list() 
for(id in 1:sexMatch_types_num){
  temp_table <- highest_IDs_table[which(highest_IDs_table$SexMatch == SexMatch[id]),]
  tbtb10 <- as.data.frame(table(temp_table$DP_MATCH, temp_table$OUTCOME, temp_table$RDP_MIHA))
  colnames(tbtb10) <- c("DP_MATCH", "Outcome", "RDP_MIHA", "count")
  
  temp_gg <- ggplot(tbtb10, aes(x = factor(RDP_MIHA), y = count, fill= Outcome , color= Outcome)) +
    geom_bar(stat = "identity", position=position_dodge()) + 
    geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
    facet_grid(. ~ DP_MATCH) +
    scale_y_continuous("Count of cases",limits=c(0, 20), breaks=seq(0, 20, 5)) + 
    scale_x_discrete("MiHA peptide in Recipient") + 
    labs(title = sprintf("%s", SexMatch[id])) 
  
  gg10[[id]] <-  temp_gg
  
}

g10 <- grid.arrange(gg10[[1]],gg10[[2]], gg10[[3]], gg10[[4]],  ncol = 2, nrow = 2, top = "Fig.9 - MiHA in Recipient vs DPB1 match")

######## DPB1*04 vs outcome vs RDP_MiHA vs Sex-mismatch
# library(gridExtra)
SexMatch <- unique(highest_IDs_table$SexMatch)
sexMatch_types_num <- length(SexMatch)
gg11 <- list() 
for(id in 1:sexMatch_types_num){
  temp_table <- highest_IDs_table[which(highest_IDs_table$SexMatch == SexMatch[id]),]
  tbtb11 <- as.data.frame(table(temp_table$X04_MATCH, temp_table$OUTCOME, temp_table$RDP_MIHA))
  colnames(tbtb11) <- c("X04_MATCH", "Outcome", "RDP_MIHA", "count")
  
  temp_gg <- ggplot(tbtb11, aes(x = factor(RDP_MIHA), y = count, fill= Outcome , color= Outcome)) +
    geom_bar(stat = "identity", position=position_dodge()) + 
    geom_text(aes(y=count, ymax=count, label=count),position= position_dodge(width=0.9), vjust=-.5) +
    facet_grid(. ~ X04_MATCH) +
    scale_y_continuous("Count of cases",limits=c(0, 20), breaks=seq(0, 20, 5)) + 
    scale_x_discrete("MiHA peptide in Recipient") + 
    labs(title = sprintf("%s", SexMatch[id])) 
  
  gg11[[id]] <-  temp_gg
  
}
g11 <- grid.arrange(gg11[[1]],gg11[[2]], gg11[[3]], gg11[[4]],  ncol = 2, nrow = 2, top = "Fig.10 - MiHA in Recipient vs DPB1*04 match")

########
HLA_typings <- read.csv(file = "../HLI_hla_mg_v3.csv")
load("../Data/ID_table.RData")

available_IDs <-ID_table[ID_table$GroupID %in% ID_table$GroupID[duplicated(ID_table$GroupID)],]

avaialble_IDs_HLA_typing <- HLA_typings[(HLA_typings$nmdp_rid %in% available_IDs$R_D_ID),]
num_groups <- dim(avaialble_IDs_HLA_typing)[1]
group_type <- sapply(1:num_groups, function(x) available_IDs$Group[available_IDs$R_D_ID %in% avaialble_IDs_HLA_typing$nmdp_rid[x]])
group_type[group_type == "a"] <- "aGVHD"
group_type[group_type == "n"] <- "non-GVHD"

highest_IDs_table$HLA_A1 <- avaialble_IDs_HLA_typing$r_a_typ1.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_A2 <- avaialble_IDs_HLA_typing$r_a_typ2.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_B1 <- avaialble_IDs_HLA_typing$r_b_typ1.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_B2 <- avaialble_IDs_HLA_typing$r_b_typ2.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_C1 <- avaialble_IDs_HLA_typing$r_c_typ1.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_C2 <- avaialble_IDs_HLA_typing$r_c_typ2.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_DRB11 <- avaialble_IDs_HLA_typing$r_drb1_typ1.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_DRB12 <- avaialble_IDs_HLA_typing$r_drb1_typ2.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_DQB11 <- avaialble_IDs_HLA_typing$r_dqb1_typ1.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]
highest_IDs_table$HLA_DQB12 <- avaialble_IDs_HLA_typing$r_dqb1_typ2.gl[which(avaialble_IDs_HLA_typing$nmdp_rid %in% highest_IDs_table$RID)]

highest_IDs_table$HLA_A1 <- sapply(as.character(highest_IDs_table$HLA_A1), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_A2 <- sapply(as.character(highest_IDs_table$HLA_A2), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_B1 <- sapply(as.character(highest_IDs_table$HLA_B1), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_B2 <- sapply(as.character(highest_IDs_table$HLA_B2), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_C1 <- sapply(as.character(highest_IDs_table$HLA_C1), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_C2 <- sapply(as.character(highest_IDs_table$HLA_C2), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_DRB11 <- sapply(as.character(highest_IDs_table$HLA_DRB11), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_DRB12 <- sapply(as.character(highest_IDs_table$HLA_DRB12), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_DQB11 <- sapply(as.character(highest_IDs_table$HLA_DQB11), function(x) unlist(strsplit(x, "/"))[1])
highest_IDs_table$HLA_DQB12 <- sapply(as.character(highest_IDs_table$HLA_DQB12), function(x) unlist(strsplit(x, "/"))[1])

get_two_field <- function(HLA_typing){
  
  reformat_typing <- unlist(strsplit(HLA_typing, ":"))
  if(length(reformat_typing) > 2){
    
    two_fields <- paste0(reformat_typing[c(1,2)], collapse = ":")
    
  }else two_fields <- HLA_typing
  
  return(two_fields)
}

highest_IDs_table$HLA_A1 <- sapply(highest_IDs_table$HLA_A1, get_two_field)
highest_IDs_table$HLA_A2 <- sapply(highest_IDs_table$HLA_A2, get_two_field)
highest_IDs_table$HLA_B1 <- sapply(highest_IDs_table$HLA_B1, get_two_field)
highest_IDs_table$HLA_B2 <- sapply(highest_IDs_table$HLA_B2, get_two_field)
highest_IDs_table$HLA_C1 <- sapply(highest_IDs_table$HLA_C1, get_two_field)
highest_IDs_table$HLA_C2 <- sapply(highest_IDs_table$HLA_C2, get_two_field)
highest_IDs_table$HLA_DRB11 <- sapply(highest_IDs_table$HLA_DRB11, get_two_field)
highest_IDs_table$HLA_DRB12 <- sapply(highest_IDs_table$HLA_DRB12, get_two_field)
highest_IDs_table$HLA_DQB11 <- sapply(highest_IDs_table$HLA_DQB11, get_two_field)
highest_IDs_table$HLA_DQB12 <- sapply(highest_IDs_table$HLA_DQB12, get_two_field)


A_0201_ID <- sort(union(which(highest_IDs_table$HLA_A1 %in% "02:01"), which(highest_IDs_table$HLA_A2 %in% "02:01")))
table(highest_IDs_table$HLA_B1[A_0201_ID])

A_0201_B_0702_ID2 <- sort(union(which(highest_IDs_table$HLA_B1[A_0201_ID] %in% "07:02"), which(highest_IDs_table$HLA_B2[A_0201_ID] %in% "07:02")))
A_0201_B_0702_ID <- A_0201_ID[A_0201_B_0702_ID2]

table(highest_IDs_table$HLA_C1[A_0201_B_0702_ID])

#### 
# ggsave(filename = "../ClinVar/DPB1-MiHA/summary_plots.pdf", plot = list(g1, g2, g3, g4, g5, g6), device = "pdf")
pdf("../ClinVar/DPB1-MiHA/DPB1_MiHA_summary_plots_01_29.pdf", width= 20, height=7)
print(g1)
print(g7)
print(g2)
print(g3)
print(g4)
print(g5)
print(g6)
print(g8)
print(g9)
plot(g10)
plot(g11)
dev.off()

#####################
# MiHA table - DPB1 permissibility
#####################
MiHA_table <- read.csv("../ClinVar/")


