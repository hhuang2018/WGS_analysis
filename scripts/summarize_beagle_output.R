source("util.R")

#### Summarize Beagles output

beagles_output_dir <- "../Output/Beagle_output/"

IBD_file <- "aGVHD_all_pairs_IBD.RData"

load(file = paste0(beagles_output_dir, IBD_file))
aGVHD_table <- new_IBD_table

aSampleID1 <- as.data.frame(t(as.data.frame(strsplit(aGVHD_table$SampleID1, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)
aSampleID2 <- as.data.frame(t(as.data.frame(strsplit(aGVHD_table$SampleID2, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)

aMatched_pair_ID <- which(aSampleID1$GroupID == aSampleID2$GroupID)
aMatched_pair_ID <- aMatched_pair_ID[which(sapply(1:length(aMatched_pair_ID), function(x) aSampleID1$R_D[aMatched_pair_ID[x]]!=aSampleID2$R_D[aMatched_pair_ID[x]]))]

aRandom_pair_ID <- which(aSampleID1$GroupID != aSampleID2$GroupID)
aRandom_pair_ID <- aRandom_pair_ID[which(sapply(1:length(aRandom_pair_ID), function(x) aSampleID1$R_D[aRandom_pair_ID[x]]!=aSampleID2$R_D[aRandom_pair_ID[x]]))]

IBD_file <- "nGVHD_all_pairs_IBD.RData"

load(file = paste0(beagles_output_dir, IBD_file))
nGVHD_table <- new_IBD_table

nSampleID1 <- as.data.frame(t(as.data.frame(strsplit(nGVHD_table$SampleID1, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)
nSampleID2 <- as.data.frame(t(as.data.frame(strsplit(nGVHD_table$SampleID2, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)

nMatched_pair_ID <- which(aSampleID1$GroupID == aSampleID2$GroupID)
nMatched_pair_ID <- nMatched_pair_ID[which(sapply(1:length(nMatched_pair_ID), function(x) nSampleID1$R_D[nMatched_pair_ID[x]]!=nSampleID2$R_D[nMatched_pair_ID[x]]))]

nRandom_pair_ID <- which(aSampleID1$GroupID != aSampleID2$GroupID)
nRandom_pair_ID <- nRandom_pair_ID[which(sapply(1:length(nRandom_pair_ID), function(x) nSampleID1$R_D[nRandom_pair_ID[x]]!=nSampleID2$R_D[nRandom_pair_ID[x]]))]

# aGVHD vs. nGVHD
plot_LOD_boxplot(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID,
                 nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                 Region = "all")

plot_LOD_histogram(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID, 
                   nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                   Region = "all")
# GVHD_LOD_mat <- data.frame(LOD = c(aGVHD_table$LOD[aMatched_pair_ID],nGVHD_table$LOD[nMatched_pair_ID]), 
#                            Group = c(rep("aGVHD", length(aMatched_pair_ID)), 
#                                      rep("non-GVHD", length(nMatched_pair_ID))),
#                            stringsAsFactors = F)
# 
# pp <- ggplot(GVHD_LOD_mat, aes(x = factor(Group), y = LOD))
# pp + geom_boxplot(aes(fill = Group)) +
#   ggtitle("IBD score (LOD) distribution on Chr6") +
#   labs(x="Group", y = "LOD")
# 
# # Matching R-D pairs vs. random R-D pairs
# LOD_score_mat <- data.frame(LOD = c(aGVHD_table$LOD[aMatched_pair_ID],nGVHD_table$LOD[nMatched_pair_ID], 
#                                     aGVHD_table$LOD[aRandom_pair_ID], nGVHD_table$LOD[nRandom_pair_ID]), 
#                             Group = c(rep("Matching Recipient-Donor Pair", length(aMatched_pair_ID)+length(nMatched_pair_ID)), 
#                                       rep("Random Recipient-Donor Pair", length(aRandom_pair_ID)+length(nRandom_pair_ID))),
#                             stringsAsFactors = F)
# colnames(LOD_score_mat)
# library(ggplot2)
# 
# p <- ggplot(LOD_score_mat, aes(x = factor(Group), y = LOD))
# p + geom_boxplot(aes(fill = Group)) + 
#   ggtitle("IBD score (LOD) distribution on Chr6") +
#   # geom_density(stat = "density") +
#   labs(x="Group", y = "LOD")

## MHC regions
plot_LOD_boxplot(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID,
                 nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                 Region = "HLA-A")
plot_LOD_boxplot(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID,
                 nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                 Region = "HLA-B")
plot_LOD_boxplot(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID,
                 nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                 Region = "HLA-C")
plot_LOD_boxplot(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID,
                 nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                 Region = "HLA-DRB1")
plot_LOD_boxplot(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID,
                 nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                 Region = "HLA-DQB1")
plot_LOD_boxplot(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID,
                 nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                 Region = "HLA-DPB1")

plot_LOD_histogram(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID, 
                   nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                   Region = "HLA-A")
plot_LOD_histogram(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID, 
                   nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                   Region = "HLA-B")
plot_LOD_histogram(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID, 
                   nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                   Region = "HLA-C")
plot_LOD_histogram(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID, 
                   nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                   Region = "HLA-DRB1")
plot_LOD_histogram(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID, 
                   nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                   Region = "HLA-DQB1")
plot_LOD_histogram(aGVHD_table, aMatched_pair_ID, aRandom_pair_ID, 
                   nGVHD_table, nMatched_pair_ID, nRandom_pair_ID, 
                   Region = "HLA-DPB1")

# a_HLA_A_IDs <- which(grepl("HLA-A", aGVHD_table$GeneNames[aMatched_pair_ID]))
# n_HLA_A_IDs <- which(grepl("HLA-A", nGVHD_table$GeneNames[nMatched_pair_ID]))
# 
# a_random_HLA_A_IDs <- which(grepl("HLA-A", aGVHD_table$GeneNames[aRandom_pair_ID]))
# n_random_HLA_A_IDs <- which(grepl("HLA-A", nGVHD_table$GeneNames[nRandom_pair_ID]))
# 
# HLA_LOD_score_mat <- data.frame(LOD = c(aGVHD_table$LOD[aMatched_pair_ID[a_HLA_A_IDs]],nGVHD_table$LOD[nMatched_pair_ID[n_HLA_A_IDs]], 
#                                     aGVHD_table$LOD[aRandom_pair_ID[a_random_HLA_A_IDs[1:length(a_HLA_A_IDs)]]], nGVHD_table$LOD[nRandom_pair_ID[1:length(n_HLA_A_IDs)]]), 
#                             Group = c(rep("Matching R-D Pair", length(a_HLA_A_IDs)+length(n_HLA_A_IDs)), 
#                                       rep("Random R-D Pair", length(a_HLA_A_IDs)+length(n_HLA_A_IDs))),
#                             stringsAsFactors = F)
# 
# p1 <- ggplot(HLA_LOD_score_mat, aes(x = factor(Group), y = LOD))
# p1 + geom_boxplot(aes(fill = Group)) + 
#   ggtitle("IBD score (LOD) distribution on HLA-A") +
#   # geom_density(stat = "density") +
#   labs(x="Group", y = "LOD")
