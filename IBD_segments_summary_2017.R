source("util.R")
library(BSgenome.Hsapiens.UCSC.hg38)
# chr.len <- c(248956422, # chr1
#                 242193529, # chr2
#                 198295559, # chr3
#                 190214555, # chr4
#                 181538259, # chr5
#                 170805979, # chr6
#                 159345973, # chr7
#                 145138636, # chr8
#                 138394717, # chr9
#                 133797422, # chr10
#                 135086622, # chr11
#                 133275309, # chr12
#                 114364328, # chr13
#                 107043718, # chr14
#                 101991189, # chr15
#                 90338345,  # chr16
#                 83257441,  # chr17
#                 80373285,  # chr18
#                 58617616,  # chr19
#                 64444167,  # chr20
#                 46709983,  # chr21
#                 50818468   # chr22
# )

chr.len = seqlengths(BSgenome.Hsapiens.UCSC.hg38)  # get chromosome lengths
# remove X,Y,M and random chromosomes
chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]
#### 

##############################
chr <- (1:22)
IBD_reformated_dir <- "../Output/IBDseq/cloud_Rformat/"
# IBD_reformated_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/summary/"

stats_table <- data.frame(CHROM = character(length(chr)*2),
                          MeanIBD = numeric(length(chr)*2),
                          MedianIBD = numeric(length(chr)*2),
                          ChrLength = numeric(length(chr)*2),
                          Group = character(length(chr)*2),
                          stringsAsFactors = F)
# stats_table <- cbind(stats_table)
# stats_table <- list()


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
counter <- 0
for (chr_id in chr) {
  
  load(paste0(IBD_reformated_dir, "summaryIBDseq_summary_chr_", chr_id, ".RData"))
  
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
                        MeanIBD = mean(stats_table))

################# MHC region 
#     HLA-A: chr6:29,941,260-29,945,884; 
#     HLA-B: chr6:31,353,872-31,357,187
#     HLA-C: chr6:31,268,749-31,272,086
#  HLA-DRB1: chr6:32,578,769-32,589,848 
#  HLA-DQB1: chr6:32,659,467-32,668,383
MHC_region_index <- data.frame(HLA = c("A", "B", "C", "DRB1", "DQB1"),
                               startID = c(29941260, 31353872, 31268749, 32578769, 32659467),
                               endID = c(29945884, 31357187, 31272086, 32589848, 32668383),
                               stringsAsFactors = F)
MHC_region_index <- MHC_region_index[order(MHC_region_index$startID), ]
MHC_region_index_vector <- sort(c(MHC_region_index$startID, MHC_region_index$endID), decreasing = F)

MHC_ARS_index <- data.frame(HLA_ARS = c("A2", "A3", "C2", "C3", "B2", "B3", "DRB1", "DQB1"),
                            startID = c(29942757, 29943268, 31271599, 31271073, 31356688,31356167, 32584109, 32664798),
                            endID = c(29943026, 29943543, 31271868, 31271348, 31356957, 31356442, 32584378, 32665067), 
                            stringsAsFactors = F)
MHC_ARS_index <- MHC_ARS_index[order(MHC_ARS_index$startID),]
MHC_ARS_index_vector <- sort(c(MHC_ARS_index$startID, MHC_ARS_index$endID))

MHC_random_Percent <- data.frame(SampleID1 = Random_high_pert$SampleID1, 
                                 SampleID2 = Random_high_pert$SampleID2,
                                 HLA.A.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 HLA.B.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 HLA.C.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 HLA.DRB1.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 HLA.DQB1.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 # IBDPercent = numeric(length(Random_high_pert$SampleID2)),
                                 stringsAsFactors = F)
MHC_matched_Percent <- data.frame(SampleID1 = Matched_high_pert$SampleID1, 
                                  SampleID2 = Matched_high_pert$SampleID2,
                                  HLA.A.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  HLA.B.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  HLA.C.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  HLA.DRB1.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  HLA.DQB1.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  stringsAsFactors = F)
IBDseq_output_dir <- "../Output/IBDseq/R_reformated/"
# IBDseq_output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/"
IBD_file <- "ibdseq_output_all_chr6_IBD.RData"

load(file = paste0(IBDseq_output_dir, IBD_file))

num_region <- dim(MHC_region_index)[1]

SampleID1 <- as.data.frame(t(as.data.frame(strsplit(new_IBD_table$SampleID1, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)
SampleID2 <- as.data.frame(t(as.data.frame(strsplit(new_IBD_table$SampleID2, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)

Matched_pair_ID <- which(SampleID1$GroupID == SampleID2$GroupID)
Matched_pair_ID <- Matched_pair_ID[which(sapply(1:length(Matched_pair_ID), function(x) SampleID1$R_D[Matched_pair_ID[x]]!=SampleID2$R_D[Matched_pair_ID[x]]))]

Random_pair_ID <- which(SampleID1$GroupID != SampleID2$GroupID)
Random_pair_ID <- Random_pair_ID[which(sapply(1:length(Random_pair_ID), function(x) SampleID1$R_D[Random_pair_ID[x]]!=SampleID2$R_D[Random_pair_ID[x]]))]

Seg_summary <- aggregate(numdup ~., data=transform(new_IBD_table[Matched_pair_ID, c(1, 3)], numdup=1), length)

sorted_table <- new_IBD_table[order(new_IBD_table$StartID), ]
startID <- unique(sorted_table$StartID)
endID <- unique(sorted_table$EndID)

start_interval_flags <-findInterval(startID, MHC_region_index_vector)
end_interval_flags <- findInterval(endID, MHC_region_index_vector)

start_flag_id <- sapply(1:length(startID), function(x) which(new_IBD_table$StartID %in% startID[x]))
end_flag_id <- sapply(1:length(endID), function(x) which(new_IBD_table$EndID %in% endID[x]))

num_flag_ids <- length(startID)
new_IBD_table$StartID_flag <- NA
for(id in 1:num_flag_ids){
  
  new_IBD_table$StartID_flag[start_flag_id[[id]]] <- start_interval_flags[id]
  
}

num_flag_ids <- length(endID)
new_IBD_table$EndID_flag <- NA
for(id in 1:num_flag_ids){
  
  new_IBD_table$EndID_flag[end_flag_id[[id]]] <- end_interval_flags[id]
  
}

#### ARS region
ARS_start_interval_flags <- findInterval(startID, MHC_ARS_index_vector)
ARS_end_interval_flags <- findInterval(endID, MHC_ARS_index_vector)

num_flag_ids <- length(startID)
new_IBD_table$StartID_ARS_flag <- NA
for(id in 1:num_flag_ids){
  
  new_IBD_table$StartID_ARS_flag[start_flag_id[[id]]] <- ARS_start_interval_flags[id]
  
}

num_flag_ids <- length(endID)
new_IBD_table$EndID_ARS_flag <- NA
for(id in 1:num_flag_ids){
  
  new_IBD_table$EndID_ARS_flag[end_flag_id[[id]]] <- ARS_end_interval_flags[id]
  
}

######### 
load("../Data/IBD_chr6_wMHC_ARS_2.RData")

MHC_IBD_table <- new_IBD_table
MHC_IBD_matched_id <- which(((new_IBD_table$StartID_flag %% 2 == 1) | (new_IBD_table$EndID_flag %% 2 == 1)) | 
                              (floor((new_IBD_table$EndID_flag - new_IBD_table$EndID_flag)/2)>0) ) # multigenes

# MHC_IBD_table <- MHC_IBD_table[MHC_IBD_matched_id, ]

MHC_region_number <- length(MHC_IBD_matched_id)
MHC_IBD_table$MHC_IBD_length <- 0
MHC_IBD_table$ARS_IBD_length <- 0
for(id in 1: MHC_region_number){
  
  # MHC region IBD length
  record <- MHC_IBD_table[MHC_IBD_matched_id[id], ]
  if(record$StartID_flag == record$EndID_flag){ # within the range
    
    MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] <- min(MHC_region_index_vector[record$EndID_flag], record$EndID) - 
                                                               max(MHC_region_index_vector[record$StartID_flag], record$StartID) + 1
    
  }else{ # partly within the range
    cross_genes <- floor((record$EndID_flag - record$StartID_flag + 1) / 2) # number of crossed regions
    if(cross_genes == 1){ # only within one MHC gene
      
      if(record$StartID_flag %% 2 == 1){ # startID is within the range
        
        MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] <- MHC_region_index_vector[record$EndID_flag] - record$StartID + 1
        
      }else{ # endID is within the range
        
        MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] <- record$EndID - MHC_region_index_vector[record$StartID_flag+1] + 1
        
      }
      
    }else{ # cross multiple regions
      
      if(record$StartID_flag %% 2 == 1) {
        start_regions_id <- record$StartID_flag + 1
      }else start_regions_id <- record$StartID_flag 
      if(record$EndID_flag %% 2 == 1){
        end_regions_id <- record$EndID_flag
      }else end_regions_id <- record$EndID_flag - 1
      
      for(region_id in 1:cross_genes){
        if(region_id == 1){
          MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] <- MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] + 
            MHC_region_index_vector[start_regions_id+1] - MHC_region_index_vector[start_regions_id] + 1
        }else if(region_id == cross_genes){
          MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] <- MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] + 
            MHC_region_index_vector[end_regions_id] - MHC_region_index_vector[start_regions_id] + 1
        }
        MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] <- MHC_IBD_table$MHC_IBD_length[MHC_IBD_matched_id[id]] +
          MHC_IBD_matched_id
        
      }
      
    }
      
    
    
  }
  MHC_region_index
  
  
}

######
for(id in 1:num_region){

  aGVHD_group_segments <- Seg_summary[which(grepl("a.",Seg_summary[, 1])),3]
  nGVHD_group_segments <- Seg_summary[which(grepl("n.",Seg_summary[, 1])),3]
  all_matched_segments <- Seg_summary[, 3]
  
  aa <- unique(new_IBD_table[Random_pair_ID, c(1, 3)])
  random_Seg_summary <- aggregate(numdup ~., data=transform(new_IBD_table[Random_pair_ID, c(1, 3)], numdup=1), length)
  
  random_pair_segments <- random_Seg_summary[, 3]
  
  All_segs_summary <- aggregate(numdup ~., data=transform(new_IBD_table[, c(1, 3)], numdup=1), length)
  all_pairs_segments <- All_segs_summary[, 3]
  
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
  
  mc_length <- data.frame(Proportion = c(aGVHD_group_length$Percent, nGVHD_group_length$Percent, 
                                         Matched_group_length$Percent, Random_group_length$Percent),
                          Group = c(rep("aGVHD", length(aGVHD_group_length$Percent)), rep("non-GVHD", length(nGVHD_group_length$Percent)),
                                    rep("Matched", length(Matched_group_length$Percent)), rep("Random Pairs", length( Random_group_length$Percent))),
                          stringsAsFactors = F)
  mc_length$Group <- factor(mc_length$Group, levels = c("aGVHD", "non-GVHD", "Matched", "Random Pairs"))
  
  p2 <- ggplot(mc_length, aes(x = Group , y = Percent, fill = Group)) + 
    geom_boxplot() +
    guides(fill=FALSE) + 
    coord_flip() +
    ggtitle(paste0("Length (%) of IBD segments \non Chromosome ", chr))
  
  # ordered
  ordered_Rand <- Random_group_length[order(Random_group_length$Percent, decreasing = T), ] 
  Random_high_pert <- ordered_Rand[which(ordered_Rand$Percent > 0), ]
  ordered_Matched <- Matched_group_length[order(Matched_group_length$Percent, decreasing = T), ]
  Matched_high_pert <- ordered_Matched[which(ordered_Matched$Percent > 0), ]
  # }
  save(mc, mc_length, p1, p2, Random_high_pert, Matched_high_pert, file = paste0(IBDseq_summary_output, "IBDseq_summary_chr_", chr, ".RData"))
}


#################
# Outlier analysis -- most likely related pairs,
# 
library(ggplot2)
NumSeg_percentiles <- quantile(mc$NumSegments)
NumSeg_Q1 <- as.numeric(NumSeg_percentiles["25%"])
NumSeg_Q3 <- as.numeric(NumSeg_percentiles["75%"])
NumSeg_IQ <- NumSeg_Q3 - NumSeg_Q1
NumSeg_outlier_th <- NumSeg_Q3 + 3 * NumSeg_IQ

NumSeg_outliers_index <- which(mc$NumSegments >= NumSeg_outlier_th)

mc_new <- mc[-which(mc$NumSegments >= 200),]

p1 <- ggplot(mc_new, aes(x = Group , y = NumSegments, fill = Group)) + 
  geom_boxplot() +
  guides(fill=FALSE) + 
  coord_flip() +
  ggtitle(paste0("Number of IBD segments \non Chromosome ", chr))

t.test(mc$NumSegments[mc$Group == "Matched"], mc$NumSegments[mc$Group == "Random Pairs"])

