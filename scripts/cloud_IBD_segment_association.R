IBD_reformatted_fp <- '/Users/hhuang2 (Deleted)/Documents/NGSProject/HLI/Output/IBDseq/R_reformated/'

chrom <- 'chr1'
load(paste0(IBD_reformatted_fp, 'ibdseq_output_all_', chrom,'_IBD.RData')) # new_IBD_table

new_IBD_table$length <- new_IBD_table$EndID - new_IBD_table$StartID + 1

SampleID1 <- as.data.frame(t(as.data.frame(strsplit(new_IBD_table$SampleID1, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)
SampleID2 <- as.data.frame(t(as.data.frame(strsplit(new_IBD_table$SampleID2, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)

Matched_pair_ID <- which(SampleID1$GroupID == SampleID2$GroupID)
Matched_pair_ID <- Matched_pair_ID[which(sapply(1:length(Matched_pair_ID), function(x) SampleID1$R_D[Matched_pair_ID[x]]!=SampleID2$R_D[Matched_pair_ID[x]]))]

Random_pair_ID <- which(SampleID1$GroupID != SampleID2$GroupID)
Random_pair_ID <- Random_pair_ID[which(sapply(1:length(Random_pair_ID), function(x) SampleID1$R_D[Random_pair_ID[x]]!=SampleID2$R_D[Random_pair_ID[x]]))]

###
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
##
library(reshape)
IBD_count_fp <- '/Users/hhuang2 (Deleted)/Documents/NGSProject/2018WGS/Data/HLI/ibdseq_output_all_chr22_segment_count.RData'
load(IBD_count_fp)

chr <- 22
df_IBD <- data.frame(POS=1:chrLengths[chr], 
                     aGVHD=vector(mode = "numeric", length = chrLengths[chr]),
                     nGVHD=vector(mode = "numeric", length = chrLengths[chr]),
                     random=vector(mode = "numeric", length = chrLengths[chr]))

#aa <- melt(df_IBD, id='POS')

MatchedPair_IBD_segment_counter <- vector(mode = "numeric", length = chrLengths[chr])
  RandomPair_IBD_segment_counter <- vector(mode = "numeric", length = chrLengths[chr])
  aGVHD_IBD_segment_counter <- vector(mode = "numeric", length = chrLengths[chr])
  nGVHD_IBD_segment_counter <- vector(mode = "numeric", length = chrLengths[chr])

# 