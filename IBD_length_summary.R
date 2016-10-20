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

# # IBD1 segment length
# aGVHD_sample1 <- unique(aGVHD_table$SampleID1)
# aGVHD_sample2 <- unique(aGVHD_table$SampleID2)
# aGVHD_num_sample1 <- length(aGVHD_sample1)
# aGVHD_num_sample2 <- length(aGVHD_sample2)
# 
# IBD_matrix <- matrix(data = 0, nrow = aGVHD_num_sample1, ncol = aGVHD_num_sample1)
# rownames(IBD_matrix) <- aGVHD_sample1
# colnames(IBD_matrix) <- aGVHD_sample1
# 
# for(id in 1:aGVHD_num_sample1){
#   
#   samp1 <- aGVHD_sample1[id]
#   index1 <- which(aGVHD_table$SampleID1 == samp1)
#   
#   for(jd in 1:aGVHD_num_sample2){
#     
#     samp2 <- aGVHD_sample2[jd]
#     index2 <- which(aGVHD_table$SampleID2[index1] == samp2)
#     
#     IBD_matrix[samp1, samp2] <- sum(aGVHD_table$EndID[index1[index2]] - aGVHD_table$StartID[index1[index2]] + 1) 
#   }
#   
# }

chrLength <- 170805979 # chr6

aGVHD_IBD <- IBD_segment_matrix(aGVHD_table, chrLength)

a_col_index <- sapply(1:dim(aGVHD_IBD)[1], function(x) which.max(aGVHD_IBD[x, ]))
aGVHD_maxIBD_pairs <- data.frame(Sample1 = rownames(aGVHD_IBD),
                                 Sample2 = colnames(aGVHD_IBD)[a_col_index],
                                 stringsAsFactors = FALSE)

aGVHD_same_group <- which(sapply(1:length(a_col_index), function(x) unlist(strsplit(aGVHD_maxIBD_pairs$Sample1[x], "\\."))[2] == unlist(strsplit(aGVHD_maxIBD_pairs$Sample2[x], "\\."))[2]))
length(aGVHD_same_group)


nGVHD_IBD <- IBD_segment_matrix(nGVHD_table, chrLength)
n_col_index <- sapply(1:dim(nGVHD_IBD)[1], function(x) which.max(nGVHD_IBD[x, ]))
nGVHD_maxIBD_pairs <- data.frame(Sample1 <- rownames(nGVHD_IBD),
                                 Sample2 <- colnames(nGVHD_IBD)[n_col_index],
                                 stringsAsFactors = FALSE)
nGVHD_same_group <- which(sapply(1:length(n_col_index), function(x) unlist(strsplit(nGVHD_maxIBD_pairs$Sample1[x], "\\."))[2] == unlist(strsplit(nGVHD_maxIBD_pairs$Sample2[x], "\\."))[2]))
length(nGVHD_same_group)


###########################
# number of IBD segments
###########################
