source("util.R")
library(BSgenome.Hsapiens.UCSC.hg38)

chr.len = seqlengths(BSgenome.Hsapiens.UCSC.hg38)  # get chromosome lengths
# remove X,Y,M and random chromosomes
chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]

##############################
chr <- (1:22)[-11]
IBD_reformated_dir <- "../Output/IBDseq/cloud_Rformat/"
stats_table <- data.frame(CHROM = character(length(chr)*2),
                          MeanIBD = numeric(length(chr)*2),
                          MedianIBD = numeric(length(chr)*2),
                          ChrLength = numeric(length(chr)*2),
                          Group = character(length(chr)*2),
                          stringsAsFactors = F)
# stats_table <- cbind(stats_table)
# stats_table <- list()
counter <- 0

load(paste0(IBD_reformated_dir, "summaryIBDseq_summary_chr_1.RData"))
all_random_Percent <- data.frame(SampleID1 = Random_high_pert$SampleID1, 
                                 SampleID2 = Random_high_pert$SampleID2,
                                 IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 IBDPercent = numeric(length(Random_high_pert$SampleID2)),
                                 stringsAsFactors = F)
all_matched_Percent <- data.frame(SampleID1 = Matched_high_pert$SampleID1, 
                                  SampleID2 = Matched_high_pert$SampleID2,
                                  IBDLength = numeric(length(Matched_high_pert$SampleID1)),
                                  IBDPercent = numeric(length(Matched_high_pert$SampleID1)),
                                  stringsAsFactors = F)
for (chr_id in chr) {
  
  load(paste0(IBD_reformated_dir, "summaryIBDseq_summary_chr_", chr_id, ".RData"))
  
  # temp_stats <- data.frame(CHROM = character(2),
  #                          MeanIBD = numeric(2),
  #                          MedianIBD = numeric(2),
  #                          stringsAsFactors = F)
  # 
  counter <- counter + 1
  stats_table$CHROM[(2*(counter-1)+1) : (2*counter)] <- paste0("chr", chr_id)
  stats_table$ChrLength[(2*(counter-1)+1) : (2*counter)] <- chr.len[paste0("chr", chr_id)]
  
  stats_table$MeanIBD[2*(counter-1)+1] <- mean(Random_high_pert$Percent)
  stats_table$MedianIBD[2*(counter-1)+1] <- median(Random_high_pert$Percent)
  stats_table$Group[2*(counter-1)+1] <- "RandomPair"
  
  stats_table$MeanIBD[2*counter] <- mean(Matched_high_pert$Percent)
  stats_table$MedianIBD[2*counter] <- median(Matched_high_pert$Percent)
  stats_table$Group[2*counter] <- "HLAMatched"
  
  ######## all length
  # random pairs
  existing_random_pair_IDs <- paste0(all_random_Percent$SampleID1, all_random_Percent$SampleID2)
  new_random_pair_IDs <- paste0(Random_high_pert$SampleID1, Random_high_pert$SampleID2)
  
  IDs_in_existing <- which(existing_random_pair_IDs %in% new_random_pair_IDs)
  IDs_in_new <- which(new_random_pair_IDs %in% existing_random_pair_IDs)
  
  ordered_IDs_in_new <- sapply(1:length(IDs_in_existing), function(x) which(new_random_pair_IDs[IDs_in_new] == existing_random_pair_IDs[IDs_in_existing[x]]))
  
  all_random_Percent$IBDLength[IDs_in_existing] <- all_random_Percent$IBDLength[IDs_in_existing] +
    Random_high_pert$Length[IDs_in_new[ordered_IDs_in_new]]
  
  if(length(IDs_in_new) < length(new_random_pair_IDs)){
    
    new_segments <- Random_high_pert[-IDs_in_new, 1:3]
    colnames(new_segments) <- c("SampleID1", "SampleID2", "IBDLength")
    new_segments$IBDPercent <- 0
    
    all_random_Percent <- rbind(all_random_Percent, new_segments)
    
  }
  
  ##### HLA mathced pairs
  existing_matched_pair_IDs <- paste0(all_matched_Percent$SampleID1, all_matched_Percent$SampleID2)
  new_matched_pair_IDs <- paste0(Matched_high_pert$SampleID1, Matched_high_pert$SampleID2)
  
  IDs_in_existing2 <- which(existing_matched_pair_IDs %in% new_matched_pair_IDs)
  IDs_in_new2 <- which(new_matched_pair_IDs %in% existing_matched_pair_IDs)
  
  ordered_IDs_in_new2 <- sapply(1:length(IDs_in_existing2), function(x) which(new_matched_pair_IDs[IDs_in_new2] == existing_matched_pair_IDs[IDs_in_existing2[x]]))
  
  all_matched_Percent$IBDLength[IDs_in_existing2] <- all_matched_Percent$IBDLength[IDs_in_existing2] +
    Matched_high_pert$Length[IDs_in_new2[ordered_IDs_in_new2]]
  
  if(length(IDs_in_new2) < length(new_matched_pair_IDs)){
    new_segments2 <- Matched_high_pert[-IDs_in_new2, 1:3]
    colnames(new_segments2) <- c("SampleID1", "SampleID2", "IBDLength")
    new_segments2$IBDPercent <- 0
    
    all_matched_Percent <- rbind(all_matched_Percent, new_segments2)
    
  }
  
  
}

total_IBD <- data.frame(CHROM = c("Total", "Total"),
                        MeanIBD = mean(stats_table) )

