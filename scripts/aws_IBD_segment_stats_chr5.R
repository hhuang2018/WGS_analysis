# process HBD/IBD result files from BEAGLE
library(Matrix)
source('utils/util.R', echo = FALSE)

#load("../Data/ID_table.RData")

IBD_reformat_dir <- "/home/hhuang/data/IBD/IBD_seq_output/"
IBD_segment_count_output_dir <-  "/home/hhuang/data/IBD/segment_stats_winsize_1000/"

chr <- 5
window_size <- 1000

IBD_files <- paste0("ibdseq_output_all_chr",chr,"_IBD.RData")
###################
# Overall IBD segments on genome
###################
chrLengths <- c(248956422, # chr1
242193529, # chr2
198295559, # chr3
190214555, # chr4
181538259, # chr5
170805979, # chr6
159345973, # chr7
145138636, # chr8
138394717, # chr9
133797422, # chr10
135086622, # chr11
133275309, # chr12
114364328, # chr13
107043718, # chr14
101991189, # chr15
90338345,  # chr16
83257441,  # chr17
80373285,  # chr18
58617616,  # chr19
64444167,  # chr20
46709983,  # chr21
50818468   # chr22
)

ptm <- proc.time()

# IBD_table <- read.table(file = paste0(IBD_file_dir, IBD_files[fid]))
# colnames(IBD_table) <- c("SampleID1", "HapID1", "SampleID2", "HapID2", "Chr", "StartInd", "EndInd", "LOD")

load(paste0(IBD_reformat_dir, IBD_files))
# new_IBD_table
# $SampleID1
# $HapID1
# $SampleID2
# $HapID2
# $Chr
# $StartID
# $EndID
# $LOD
# $GeneNames
num_rows <- dim(new_IBD_table)[1]
new_IBD_table$Group1 <- sapply(1:num_rows, function(x) paste0(unlist(strsplit(new_IBD_table$SampleID1[x], "\\."))[1:2], collapse = '.'))
new_IBD_table$Group2 <- sapply(1:num_rows, function(x) paste0(unlist(strsplit(new_IBD_table$SampleID2[x], "\\."))[1:2], collapse = '.'))

new_IBD_table$agvhd1 <- sapply(1:num_rows, function(x) if(unlist(strsplit(new_IBD_table$SampleID1[x], "\\."))[1] == 'a') TRUE else FALSE)
new_IBD_table$agvhd2 <- sapply(1:num_rows, function(x) if(unlist(strsplit(new_IBD_table$SampleID2[x], "\\."))[1] == 'a') TRUE else FALSE)

new_IBD_table$DR1 <- sapply(1:num_rows, function(x) if(unlist(strsplit(new_IBD_table$SampleID1[x], "\\."))[3] == 'R') TRUE else FALSE)
new_IBD_table$DR2 <- sapply(1:num_rows, function(x) if(unlist(strsplit(new_IBD_table$SampleID2[x], "\\."))[3] == 'R') TRUE else FALSE)

new_IBD_table$GroupPair <- sapply(1:num_rows, function(x) paste0(sort(c(new_IBD_table$Group1[x], new_IBD_table$Group2[x])), collapse = '-'))

# $DR - Donor:FALSE; Recipient:TRUE
# $agvhd - aGVHD:TRUE ; non-GVHD:FALSE
# $group

### Matched Pairs
MatchedPairIDs_index <- which(new_IBD_table$Group1 == new_IBD_table$Group2)
MatchedPairIDs <- unique(new_IBD_table$Group1[MatchedPairIDs_index])
num_matchedPairs <- length(MatchedPairIDs)

MatchedPair_IBD <- Matrix(0, nrow = num_matchedPairs, ncol = chrLengths[chr]/window_size, sparse = T)
#MatchedPair_IBD <- as.data.frame(MatchedPair_IBD, row.names = MatchedPairIDs)
num_rows <- length(MatchedPairIDs_index)

for(id in 1:num_rows){
    
    start_index <- floor(new_IBD_table$StartID[MatchedPairIDs_index[id]]/window_size)
    end_index <- ceiling(new_IBD_table$EndID[MatchedPairIDs_index[id]]/window_size)
    
    row_index <- which(MatchedPairIDs %in% new_IBD_table$Group1[MatchedPairIDs_index[id]])
    
    MatchedPair_IBD[row_index, start_index:end_index] <- 1
}

save(MatchedPair_IBD, MatchedPairIDs, file = paste0(IBD_segment_count_output_dir, 'MatchedPair_IBD_chr', chr,'.RData'))
rm(MatchedPair_IBD)

### Random DD pairs
Random_DD_index <- which((new_IBD_table$Group1 != new_IBD_table$Group2) & (!new_IBD_table$DR1) & (!new_IBD_table$DR2))
Random_DD_IDs <- unique(new_IBD_table$GroupPair[Random_DD_index])
num_RandomDDPairs <- length(Random_DD_IDs)

Random_DD_IBD <- Matrix(0, nrow = num_RandomDDPairs, ncol = chrLengths[chr]/window_size, sparse = T)

for(id in 1:num_rows){
    
    start_index <- floor(new_IBD_table$StartID[Random_DD_index[id]]/window_size)
    end_index <- ceiling(new_IBD_table$EndID[Random_DD_index[id]]/window_size)
    
    row_index <- which(Random_DD_IDs %in% new_IBD_table$GroupPair[Random_DD_index[id]])
    
    Random_DD_IBD[row_index, start_index:end_index] <- 1
    
}

save(Random_DD_IBD, Random_DD_IDs, file = paste0(IBD_segment_count_output_dir, 'Random_DD_IBD_chr', chr,'.RData'))
rm(Random_DD_IBD)

### Random RD pairs
Random_RD_index <- which((new_IBD_table$Group1 != new_IBD_table$Group2) &
(((!new_IBD_table$DR1) & new_IBD_table$DR2)|(!new_IBD_table$DR2) & new_IBD_table$DR1))
Random_RD_IDs <- unique(new_IBD_table$GroupPair[Random_RD_index])
num_RandomRDPairs <- length(Random_RD_IDs)

Random_RD_IBD <- Matrix(0, nrow = num_RandomRDPairs, ncol = chrLengths[chr]/window_size, sparse = T)

for(id in 1:num_rows){
    
    start_index <- floor(new_IBD_table$StartID[Random_RD_index[id]]/window_size)
    end_index <- ceiling(new_IBD_table$EndID[Random_RD_index[id]]/window_size)
    
    row_index <- which(Random_RD_IDs %in% new_IBD_table$GroupPair[Random_RD_index[id]])
    
    Random_RD_IBD[row_index, start_index:end_index] <- 1
    
}

save(Random_RD_IBD, Random_RD_IDs, file = paste0(IBD_segment_count_output_dir, 'Random_RD_IBD_chr', chr,'.RData'))
rm(Random_RD_IBD)

### Random RR pairs
Random_RR_index <- which((new_IBD_table$Group1 != new_IBD_table$Group2) & (new_IBD_table$DR1 & new_IBD_table$DR2))
Random_RR_IDs <- unique(new_IBD_table$GroupPair[Random_RR_index])
num_RandomRR_Pairs <- length(Random_RR_IDs)

Random_RR_IBD <- Matrix(0, nrow = num_RandomRR_Pairs, ncol = chrLengths[chr]/window_size, sparse = T)

for(id in 1:num_rows){
    
    start_index <- floor(new_IBD_table$StartID[Random_RR_index[id]]/window_size)
    end_index <- ceiling(new_IBD_table$EndID[Random_RR_index[id]]/window_size)
    
    row_index <- which(Random_RR_IDs %in% new_IBD_table$GroupPair[Random_RR_index[id]])
    
    Random_RR_IBD[row_index, start_index:end_index] <- 1
    
}

save(Random_RR_IBD, Random_RR_IDs, file = paste0(IBD_segment_count_output_dir, 'Random_RR_IBD_chr', chr,'.RData'))
rm(Random_RR_IBD)
