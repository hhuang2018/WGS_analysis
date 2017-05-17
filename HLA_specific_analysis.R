##### HLA restricted MiHAs
source("util.R")

library(ggplot2)
known_MiHA_table <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T)

#known_miha_freq <- read.csv("../WW_MiHA/HLA_restricted_knownMiHAs.csv", header = F)
known_miha_freq <- read.table("../WW_MiHA/Restricted_known_MiHAs.txt", header = F)
colnames(known_miha_freq) <- c("GroupType", "GroupID", "HLA_type",  "SNP", "CHROM", "REF", "ALT")

table(known_miha_freq[c("HLA_type","SNP")])

# known_miha_freq$HLA_SNP <- sapply(1:dim(known_miha_freq)[1], 
#                                   function(x) paste0(known_miha_freq$HLA_type[x], "-", known_miha_freq$SNP[x]))

table_index <- sapply(1:dim(known_miha_freq)[1], function(x) which(known_MiHA_table$SNPs %in% known_miha_freq$SNP[x]))


known_MiHA_stats <- cbind(known_miha_freq, known_MiHA_table[table_index, ])
known_MiHA_stats$MiHA_HLA_SNP_Gene <- sapply(1:dim(known_MiHA_stats)[1],
                                  function(x) paste0(known_MiHA_stats$MiHAs[x], "<>", known_MiHA_stats$HLA_type[x], "<>", known_MiHA_stats$SNP[x],
                                                     "<>", known_MiHA_stats$Gene[x]))
tabs <- as.data.frame(table(known_MiHA_stats[, c("GroupType", "MiHA_HLA_SNP_Gene")]))
tabs$GroupType <- as.character(tabs$GroupType)
tabs$GroupType[which(tabs$GroupType == "n")] <- "non-aGVHD"
tabs$GroupType[which(tabs$GroupType == "a")] <- "aGVHD"
# tabs$GroupType <- as.factor(tabs$GroupType)
tabs$Freq[which(tabs$GroupType == "non-aGVHD")] <- -tabs$Freq[which(tabs$GroupType == "non-aGVHD")]
# 
# aa <- sapply(1:74, function(x) strsplit(as.character(tabs$MiHA_HLA_SNP_Gene[x]), "<>"))
# aa <- do.call("rbind", aa)
ggplot(tabs, aes(x = HLA_SNP_Gene, ymax = Freq, ymin = 0, color = GroupType)) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  ggtitle("HLA Restricted MiHAs") +
  theme_bw() +
  xlab("") +
  # theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip() 

##### Log Odds
unique_MiHAs <- unique(tabs$MiHA_HLA_SNP_Gene)
num_unique_MiHAs <- length(unique_MiHAs)

aGVHD_index <- which(tabs$GroupType == "aGVHD")
nGVHD_index <- which(tabs$GroupType == "non-aGVHD")
# tabs$HLA_SNP_Gene <- as.character(tabs$HLA_SNP_Gene)
# n_index <- integer(length(nGVHD_index))
counter <- 0

MiHA_LLR <- data.frame(MiHA_HLA_SNP_Gene = character(num_unique_MiHAs),
                       aGVHD_count = numeric(num_unique_MiHAs),
                       nGVHD_count = numeric(num_unique_MiHAs),
                       LLR = numeric(num_unique_MiHAs),
                       stringsAsFactors = F)

for(id in 1:num_unique_MiHAs){
  # counter <- counter +1
  # n_index <- nGVHD_index[which(tabs$HLA_SNP_Gene[nGVHD_index] %in% tabs$HLA_SNP_Gene[id])]
  # LLR_odds$HLA_SNP_Gene[counter] <- as.character(tabs$HLA_SNP_Gene[id])
  
  index <- which(tabs$MiHA_HLA_SNP_Gene %in% unique_MiHAs[id])
  
  # count_summary <- as.data.frame(table(known_miha_freq$GroupType[index]))
  
  MiHA_LLR$MiHA_HLA_SNP_Gene[id] <- as.character(unique_MiHAs[id])
  MiHA_LLR$aGVHD_count[id] <- abs(tabs$Freq[index[which(tabs$GroupType[index] == "aGVHD")]])
  MiHA_LLR$nGVHD_count[id] <- abs(tabs$Freq[index[which(tabs$GroupType[index] == "non-aGVHD")]])
  # MiHA_LLR$LLR[id] <- log10(MiHA_LLR$aGVHD_count[id] / MiHA_LLR$nGVHD_count[id])
  
  if(MiHA_LLR$nGVHD_count[id] == 0){
    cat("aGVHD Only: ", as.character(MiHA_LLR$MiHA_HLA_SNP_Gene[id]), "\n") 
    MiHA_LLR$LLR[id] <- NA
    
  }else if(MiHA_LLR$aGVHD_count[id] == 0){
    cat("non-aGVHD Only: ", as.character(MiHA_LLR$MiHA_HLA_SNP_Gene[id]), "\n") 
    MiHA_LLR$LLR[id] <- NA
    
  }else{
    MiHA_LLR$LLR[id] <- log10(MiHA_LLR$aGVHD_count[id]/MiHA_LLR$nGVHD_count[id])
  }
}

write.csv(MiHA_LLR, file = "../WW_MiHA/Restricted_Known_MiHA_counts_LOD.csv", row.names = F)

MiHA_LLR <- MiHA_LLR[-which(is.na(MiHA_LLR$LLR)), ]
mc_LLR <- MiHA_LLR[order(MiHA_LLR$LLR, decreasing = TRUE), ]
mc_LLR <- within(mc_LLR, MiHA_HLA_SNP_Gene <- factor(MiHA_HLA_SNP_Gene, levels=factor(mc_LLR$MiHA_HLA_SNP_Gene)))

ggplot(mc_LLR, aes(x = MiHA_HLA_SNP_Gene, ymax = LLR, ymin = 0, color = "#D55E00")) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  ggtitle("LOD of HLA restricted MiHAs [log10(aGVHD/nGVHD)]") +
  theme_bw() +
  xlab("") +
  # scale_colour(values = "#D55E00") +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  coord_flip() 


tabs[c(25, 39, 65, 71, 73),]
aa <- sapply(1:74, function(x) strsplit(as.character(tabs$HLA_SNP_Gene[x]), "<>"))
aa <- do.call("rbind", aa)
colnames(aa) <- c("HLA_typing", "SNPs", "Gene", "MiHAs")
new_tabs <- cbind(tabs, aa)

write.csv(new_tabs, file = "../Data/HLA_restricted_MiHAs_counts.csv", row.names = F)

# pp + scale_linetype_discrete(name="GroupType")
##############
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
group_type[group_type == "n"] <- "non-aGVHD"
# num_groups <- dim(avaialble_IDs_HLA_typing)[1]

HLA_A <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_a_typ1.gl), as.character(avaialble_IDs_HLA_typing$r_a_typ2.gl)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_B <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_b_typ1.gl), as.character(avaialble_IDs_HLA_typing$r_b_typ2.gl)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_C <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_c_typ1.gl), as.character(avaialble_IDs_HLA_typing$r_c_typ2.gl)),
                    GroupType = rep(group_type, each = 2),
                    stringsAsFactors = F)
HLA_DRB1 <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_drb1_typ1.gl), as.character(avaialble_IDs_HLA_typing$r_drb1_typ2.gl)),
                       GroupType = rep(group_type, each = 2),
                       stringsAsFactors = F)
HLA_DQB1 <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_dqb1_typ1.gl), as.character(avaialble_IDs_HLA_typing$r_dqb1_typ2.gl)),
                       GroupType = rep(group_type, each = 2),
                       stringsAsFactors = F)
HLA_DPB1 <- data.frame(Typing = c(as.character(avaialble_IDs_HLA_typing$r_dpb1_typ1.gl), as.character(avaialble_IDs_HLA_typing$r_dpb1_typ2.gl)),
                       GroupType = rep(group_type, each = 2),
                       stringsAsFactors = F)

