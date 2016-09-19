source("util.R")
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

#### Summarize Beagles output

# beagles_output_dir <- "../Output/IBDseq/R_reformated/"
beagles_output_dir <- "/mnt/cloudbiodata_nfs_2/hhuang/IBD/IBD_seq_output/"
IBDseq_summary_output <- "/mnt/cloudbiodata_nfs_2/hhuang/IBD/IBD_seq_output/summary"
#chr <- 10 
chr_list <- c(1, 2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22)

for(chr in chr_list){
  
  IBD_file <- paste0("ibdseq_output_all_chr",chr,"_IBD.RData")
  
  load(file = paste0(beagles_output_dir, IBD_file))
  
  SampleID1 <- as.data.frame(t(as.data.frame(strsplit(new_IBD_table$SampleID1, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)
  SampleID2 <- as.data.frame(t(as.data.frame(strsplit(new_IBD_table$SampleID2, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)
  
  Matched_pair_ID <- which(SampleID1$GroupID == SampleID2$GroupID)
  Matched_pair_ID <- Matched_pair_ID[which(sapply(1:length(Matched_pair_ID), function(x) SampleID1$R_D[Matched_pair_ID[x]]!=SampleID2$R_D[Matched_pair_ID[x]]))]
  
  Random_pair_ID <- which(SampleID1$GroupID != SampleID2$GroupID)
  Random_pair_ID <- Random_pair_ID[which(sapply(1:length(Random_pair_ID), function(x) SampleID1$R_D[Random_pair_ID[x]]!=SampleID2$R_D[Random_pair_ID[x]]))]
  
  # Donor_random_pair_ID <- which(grepl(".D", SampleID1$R_D) && grepl(".D",SampleID2$R_D))
  
  ## Number of segments
  # aa <- unique(new_IBD_table[Matched_pair_ID, c(1, 3)])
  Seg_summary <- aggregate(numdup ~., data=transform(new_IBD_table[Matched_pair_ID, c(1, 3)], numdup=1), length)
  
  aGVHD_group_segments <- Seg_summary[which(grepl("a.",Seg_summary[, 1])),3]
  nGVHD_group_segments <- Seg_summary[which(grepl("n.",Seg_summary[, 1])),3]
  all_matched_segments <- Seg_summary[, 3]
  
  aa <- unique(new_IBD_table[Random_pair_ID, c(1, 3)])
  random_Seg_summary <- aggregate(numdup ~., data=transform(new_IBD_table[Random_pair_ID, c(1, 3)], numdup=1), length)
  
  random_pair_segments <- random_Seg_summary[, 3]
  
  All_segs_summary <- aggregate(numdup ~., data=transform(new_IBD_table[, c(1, 3)], numdup=1), length)
  all_pairs_segments <- All_segs_summary[, 3]
  
  library(ggplot2)
  
  mc <- data.frame(NumSegments = c(aGVHD_group_segments, nGVHD_group_segments, all_matched_segments, 
                                   random_pair_segments, all_pairs_segments),
                   Group = c(rep("aGVHD", length(aGVHD_group_segments)), rep("non-GVHD", length(nGVHD_group_segments)),
                             rep("Matched", length(all_matched_segments)), rep("Random Pairs", length(random_pair_segments)),
                             rep("All Pairs", length(all_pairs_segments))),
                   stringsAsFactors = F)
  mc$Group <- factor(mc$Group, levels = c("aGVHD", "non-GVHD", "Matched", "Random Pairs", "All Pairs"))
  
  p1 <- ggplot(mc, aes(x = Group , y = NumSegments, fill = Group)) + 
    geom_boxplot() +
    guides(fill=FALSE) + 
    coord_flip() +
    ggtitle(paste0("Number of IBD segments \non Chromosome ", chr))
  
  # histogram
  # ggplot(mc[1:(length(aGVHD_group_segments)+length(nGVHD_group_segments)),], aes(x = NumSegments, fill = Group)) +
  #   geom_histogram(binwidth = 1, alpha=.6, position="identity") + 
  #   ggtitle(paste0("Number of IBD segments distribution \non Chromosome", chr)) +
  #   # geom_density(stat = "density") +
  #   labs(x="Number of IBD segments", y = "Count")
  
  
  ## IBD proportion over each chromosome
  # chrLength <- 133797422 # chr10 #170805979 # chr6
  
  GVHD_IBD <- IBD_segment_matrix(new_IBD_table, chrLengths[chr])
  
  Matched_group_length <- IBD_segment_data_frame(new_IBD_table[Matched_pair_ID, ], chrLengths[chr])
  
  aGVHD_group_length <- Matched_group_length[grepl("a.", Matched_group_length$SampleID1), ]
  nGVHD_group_length <- Matched_group_length[grepl("n.", Matched_group_length$SampleID1), ]
  
  Random_group_length <- IBD_segment_data_frame(new_IBD_table[Random_pair_ID, ], chrLengths[chr])
  
  mc_length <- data.frame(Percent = c(aGVHD_group_length$Percent, nGVHD_group_length$Percent, 
                                      Matched_group_length$Percent, Random_group_length$Percent),
                          Group = c(rep("aGVHD", length(aGVHD_group_length$Percent)), rep("non-GVHD", length(nGVHD_group_length$Percent)),
                                    rep("Matched", length(Matched_group_length$Percent)), rep("Random Pairs", length( Random_group_length$Percent))),
                          stringsAsFactors = F)
  mc_length$Group <- factor(mc_length$Group, levels = c("aGVHD", "non-GVHD", "Matched", "Random Pairs"))
  
  p2 <- ggplot(mc_length, aes(x = Group , y = Percent, fill = Group)) + 
    geom_boxplot() +
    guides(fill=FALSE) + 
    coord_flip() +
    ggtitle(paste0("Length of IBD segments \non Chromosome ", chr))
  
  # ordered
  ordered_Rand <- Random_group_length[order(Random_group_length$Percent, decreasing = T), ] 
  Random_high_pert <- ordered_Rand[which(ordered_Rand$Percent > 0), ]
  ordered_Matched <- Matched_group_length[order(Matched_group_length$Percent, decreasing = T), ]
  Matched_high_pert <- ordered_Matched[which(ordered_Matched$Percent > 0), ]
  # }
  save(mc, mc_length, p1, p2, file = paste0(IBDseq_summary_output, "IBDseq_summary_chr_", chr, ".RData"))
}
