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

##################################### chr 1:22 all 
chr <- c(1:5, 7:22)
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
# chr_all_length <- 0
sink("../Output/Chr_tTest_wilcox.txt")
for (chr_id in chr) {
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

  # chr_all_length <- chr_all_length + as.numeric(chr.len[paste0("chr", chr_id)])
  # chr_all_length <- chr_all_length + chr.len[chr_id]
  
  load(paste0(IBD_reformated_dir, "summaryIBDseq_summary_chr_", chr_id, ".RData"))
  
  # counter <- counter + 1
  # stats_table$CHROM[(2*(counter-1)+1) : (2*counter)] <- paste0("chr", chr_id)
  # stats_table$ChrLength[(2*(counter-1)+1) : (2*counter)] <- chr.len[paste0("chr", chr_id)]
  # 
  # stats_table$MeanIBD[2*(counter-1)+1] <- mean(Random_high_pert$Percent)
  # stats_table$MedianIBD[2*(counter-1)+1] <- median(Random_high_pert$Percent)
  # stats_table$Group[2*(counter-1)+1] <- "RandomPair"
  # 
  # stats_table$MeanIBD[2*counter] <- mean(Matched_high_pert$Percent)
  # stats_table$MedianIBD[2*counter] <- median(Matched_high_pert$Percent)
  # stats_table$Group[2*counter] <- "HLAMatched"
  
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
  
  cat(paste0("Chromosome ", chr_id, ": \n"))
  print(t.test(all_matched_Percent$IBDLength, all_random_Percent$IBDLength))
  cat("\n-----------------------------------------\n")
  print(wilcox.test(all_matched_Percent$IBDLength, all_random_Percent$IBDLength, alternative = "two.sided"))
  
}
sink()

all_matched_Percent$IBDPercent <- all_matched_Percent$IBDLength / chr_all_length
all_random_Percent$IBDPercent <- all_random_Percent$IBDLength / chr_all_length

total_IBD <- data.frame(CHROM = c("Total", "Total"),
                        MeanIBD = mean(stats_table))

save(all_matched_Percent, all_random_Percent, file = "../Output/IBDseq/summary/all_chromosome_IBD_04062017.RData")

#########
load("../Output/IBDseq/summary/all_chromosome_IBD_04062017.RData")

ordered_Rand <- all_random_Percent[order(all_random_Percent$IBDLength, decreasing = T), ] 
all_Random_high_pert <- ordered_Rand[which(ordered_Rand$IBDPercent > 0), ]
ordered_Matched <- all_matched_Percent[order(all_matched_Percent$IBDPercent, decreasing = T), ]
all_Matched_high_pert <- ordered_Matched[which(ordered_Matched$IBDPercent > 0), ]
all_Random_high_pert[1:10, ]
all_Matched_high_pert[1:10,]

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

# Seg_summary <- aggregate(numdup ~., data=transform(new_IBD_table[Matched_pair_ID, c(1, 3)], numdup=1), length)

#### sort the position indice in an increasing order, remove duplicates
sorted_table <- new_IBD_table[order(new_IBD_table$StartID), ]
startID <- unique(sorted_table$StartID)
endID <- unique(sorted_table$EndID)

### region flags
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
MHC_ARS_index <- data.frame(HLA_ARS = c("A2", "A3", "C2", "C3", "B2", "B3", "DRB1", "DQB1"),
                            startID = c(29942757, 29943268, 31271599, 31271073, 31356688,31356167, 32584109, 32664798),
                            endID = c(29943026, 29943543, 31271868, 31271348, 31356957, 31356442, 32584378, 32665067), 
                            stringsAsFactors = F)
MHC_ARS_index <- MHC_ARS_index[order(MHC_ARS_index$startID),]
MHC_ARS_index_vector <- sort(c(MHC_ARS_index$startID, MHC_ARS_index$endID))

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
save(new_IBD_table, file = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_renamed_all_chr6/IBD/IBD_chr6_wMHC_ARS.RData")

save(new_IBD_table, file = "../Data/IBD_chr6_wMHC_ARS_2.RData")

######### 
load("../Data/IBD_chr6_wMHC_ARS_2.RData")

MHC_IBD_table <- new_IBD_table
# MHC_IBD_matched_id <- which(((new_IBD_table$StartID_flag %% 2 == 1) | (new_IBD_table$EndID_flag %% 2 == 1)) | 
#                               (floor((new_IBD_table$EndID_flag - new_IBD_table$EndID_flag)/2)>0) ) # multigenes
# rule_out_IDs <- which((MHC_IBD_table$StartID_flag %% 2 == 0)& #(MHC_IBD_table$EndID_flag %% 2 ==0) & 
#                         (MHC_IBD_table$EndID_flag == MHC_IBD_table$StartID_flag))
# MHC_IBD_table <- MHC_IBD_table[-rule_out_IDs, ]
# MHC_IBD_table <- MHC_IBD_table[MHC_IBD_matched_id, ]

MHC_region_number <- dim(MHC_IBD_table)[1]
MHC_IBD_table$MHC_IBD_length <- 0
MHC_IBD_table$ARS_IBD_length <- 0

MHC_IBD_table$MHC_A_IBD_length <- 0    # flag - 1
MHC_IBD_table$MHC_C_IBD_length <- 0    # 3
MHC_IBD_table$MHC_B_IBD_length <- 0    # 5
MHC_IBD_table$MHC_DRB1_IBD_length <- 0 # 7
MHC_IBD_table$MHC_DQB1_IBD_length <- 0 # 9

for(id in 1: MHC_region_number){
  
  # MHC region IBD length
  startID_even_odd <- MHC_IBD_table$StartID_flag[id] %% 2 ## 0: even; 1: odd
  endID_even_odd <- MHC_IBD_table$EndID_flag[id] %% 2
  
  if(startID_even_odd * endID_even_odd == 1){ # both are odd numbers
    region_number <- (MHC_IBD_table$EndID_flag[id] - MHC_IBD_table$StartID_flag[id])/2 + 1
  }else{
    region_number <- ceiling((MHC_IBD_table$EndID_flag[id] - MHC_IBD_table$StartID_flag[id])/2)
  }
  if(region_number>0){
    reg_start <- vector(mode = "numeric", length=region_number)
    reg_end <- vector(mode = "numeric", length=region_number)
    if(region_number>1){ # multiple regions
      
      ### reg_start - first ID
      if(startID_even_odd == 0) {# startID is even, meaning between regions
        
        reg_start[1] <- MHC_region_index_vector[MHC_IBD_table$StartID_flag[id]+1]
        reg_end[1] <- MHC_region_index_vector[MHC_IBD_table$StartID_flag[id]+2]
        
      }else { # startID is odd, meaning within regions.
        
        reg_start[1] <- MHC_IBD_table$StartID[id]
        reg_end[1] <- MHC_region_index_vector[MHC_IBD_table$StartID_flag[id]+1]
        
      }
      ### reg_end - last ID
      if(endID_even_odd == 0){# endID is even,meaning between regions 
        
        reg_end[region_number] <- MHC_region_index_vector[MHC_IBD_table$EndID_flag[id]]
        reg_start[region_number] <- MHC_region_index_vector[MHC_IBD_table$EndID_flag[id] - 1]
        
      }else{#  endID is odd, meaning within regions. 
        
        reg_end[region_number] <- MHC_IBD_table$EndID[id]
        reg_start[region_number] <- MHC_region_index_vector[MHC_IBD_table$EndID_flag[id]]
        
      }
      
      ####### if region_number > 2
      if(region_number > 2){
        
        for(mid_id in 2:(region_number-1)){
          
          mid_start <- which(MHC_region_index_vector %in% reg_end[mid_id - 1]) + 1
          # mid_end <- which(MHC_region_index_vector %in% reg_start[region_number - mid_id + 2]) - 1
          if(mid_start %% 2 == 0) stop("wrong index for startID.")
          reg_start[mid_id] <- MHC_region_index_vector[mid_start]
          reg_end[mid_id] <- MHC_region_index_vector[mid_start + 1]
          
        }
        
      }
      
    } else{ # only one region
      ### reg_start - first ID
      if(startID_even_odd == 0) {# startID is even, meaning between regions
        
        reg_start <- MHC_region_index_vector[MHC_IBD_table$StartID_flag[id]+1]
        # reg_end <- MHC_region_index_vector[MHC_IBD_table$StartID_flag[id]+2]
        
      } else { # startID is odd, meaning within regions.
        
        reg_start <- MHC_IBD_table$StartID[id]
        # reg_end <- MHC_region_index_vector[MHC_IBD_table$StartID_flag[id]+1]
        
      }
      ### reg_end - last ID
      if(endID_even_odd == 0){# endID is even,meaning between regions 
        
        reg_end <- MHC_region_index_vector[MHC_IBD_table$EndID_flag[id]]
        
      }else{#  endID is odd, meaning within regions. 
        
        reg_end <- MHC_IBD_table$EndID[id]
        
      }
    }
    
    MHC_IBD_table$MHC_IBD_length[id] <- sum(reg_end - reg_start) + region_number
    
    MHC_region_id <- floor(MHC_IBD_table$StartID_flag[id]/2) # A: 0; C:1 ; B:2; DRB1: 3; DQB1: 4
    MHC_IBD_table[id, (16 + MHC_region_id):(15 + MHC_region_id +region_number)] <- reg_end - reg_start + 1
  }
  
  ####################### ARS region IBD length
  startID_ARS_even_odd <- MHC_IBD_table$StartID_ARS_flag[id] %% 2 ## 0: even; 1: odd
  endID_ARS_even_odd <- MHC_IBD_table$EndID_ARS_flag[id] %% 2
  
  if(startID_ARS_even_odd * endID_ARS_even_odd == 1){ # both are odd numbers
    region_ARS_number <- (MHC_IBD_table$EndID_ARS_flag[id] - MHC_IBD_table$StartID_ARS_flag[id])/2 + 1
  }else{
    region_ARS_number <- ceiling((MHC_IBD_table$EndID_ARS_flag[id] - MHC_IBD_table$StartID_ARS_flag[id])/2)
  }
  if(region_ARS_number >0 ){
    reg_ARS_start <- vector(mode = "numeric", length=region_ARS_number)
    reg_ARS_end <- vector(mode = "numeric", length=region_ARS_number)
    if(region_ARS_number>1){ # multiple regions
      
      ### reg_ARS_start - first ID
      if(startID_ARS_even_odd == 0) {# startID is even, meaning between regions
        
        reg_ARS_start[1] <- MHC_ARS_index_vector[MHC_IBD_table$StartID_ARS_flag[id]+1]
        reg_ARS_end[1] <- MHC_ARS_index_vector[MHC_IBD_table$StartID_ARS_flag[id]+2]
        
      }else { # startID is odd, meaning within regions.
        
        reg_ARS_start[1] <- MHC_IBD_table$StartID[id]
        reg_ARS_end[1] <- MHC_ARS_index_vector[MHC_IBD_table$StartID_ARS_flag[id]+1]
        
      }
      ### reg_end - last ID
      if(endID_ARS_even_odd == 0){# endID is even,meaning between regions 
        
        reg_ARS_end[region_ARS_number] <- MHC_ARS_index_vector[MHC_IBD_table$EndID_ARS_flag[id]]
        reg_ARS_start[region_ARS_number] <- MHC_ARS_index_vector[MHC_IBD_table$EndID_ARS_flag[id] - 1]
        
      }else{#  endID is odd, meaning within regions. 
        
        reg_ARS_end[region_ARS_number] <- MHC_IBD_table$EndID[id]
        reg_ARS_start[region_ARS_number] <- MHC_ARS_index_vector[MHC_IBD_table$EndID_ARS_flag[id]]
        
      }
      
      ####### if region_number > 2
      if(region_ARS_number > 2){
        
        for(mid_id in 2:(region_ARS_number-1)){
          
          mid_start <- which(MHC_ARS_index_vector %in% reg_ARS_end[mid_id - 1]) + 1
          # mid_end <- which(MHC_region_index_vector %in% reg_start[region_number - mid_id + 2]) - 1
          if(mid_start %% 2 == 0) stop("wrong index for startID.")
          reg_ARS_start[mid_id] <- MHC_ARS_index_vector[mid_start]
          reg_ARS_end[mid_id] <- MHC_ARS_index_vector[mid_start + 1]
          
        }
        
      }
      
    }else{ # only one region
      ### reg_start - first ID
      if(startID_ARS_even_odd == 0) {# startID is even, meaning between regions
        
        reg_ARS_start <- MHC_ARS_index_vector[MHC_IBD_table$StartID_ARS_flag[id]+1]
        # reg_ARS_end <- MHC_ARS_index_vector[MHC_IBD_table$StartID_ARS_flag[id]+2]
        
      }else { # startID is odd, meaning within regions.
        
        reg_ARS_start <- MHC_IBD_table$StartID[id]
        # reg_ARS_end <- MHC_ARS_index_vector[MHC_IBD_table$StartID_ARS_flag[id]+1]
        
      }
      ### reg_end - last ID
      if(endID_ARS_even_odd == 0){# endID is even,meaning between regions 
        
        reg_ARS_end <- MHC_ARS_index_vector[MHC_IBD_table$EndID_ARS_flag[id]]
        # if(reg_end_check != reg_ARS_end) stop("Wrong/Inconsistent ARS_EndID!")
        # reg_start[region_number] <- MHC_region_index_vector[MHC_IBD_table$EndID_flag[id] - 1]
        
      }else{#  endID is odd, meaning within regions. 
        
        reg_ARS_end <- MHC_IBD_table$EndID[id]
        # if(reg_end_check != reg_ARS_end) stop("Wrong/Inconsistent ARS_EndID!")
        # reg_start[region_number] <- MHC_region_index_vector[MHC_IBD_table$EndID_flag[id]]
        
      }
    }
    
    MHC_IBD_table$ARS_IBD_length[id] <- sum(reg_ARS_end - reg_ARS_start) + region_ARS_number
  }
}

# save(MHC_IBD_table, file = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_renamed_all_chr6/IBD/IBD_info_Original_chr6_wMHC_0501.RData")

save(MHC_IBD_table, file = "../Data/IBD_info_chr6_wMHC_ARS_0406.RData")

###### Percentage MHC region
load("../Data/IBD_info_chr6_wMHC_ARS_0406.RData")

IBDseq_output_dir <-  "../Output/IBDseq/cloud_Rformat/"
# IBDseq_output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/"
IBD_file <- "summaryIBDseq_summary_chr_6.RData"  # chr 6

load(file = paste0(IBDseq_output_dir, IBD_file))

num_region <- dim(MHC_region_index)[1]

SampleID1 <- as.data.frame(t(as.data.frame(strsplit(MHC_IBD_table$SampleID1, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)
SampleID2 <- as.data.frame(t(as.data.frame(strsplit(MHC_IBD_table$SampleID2, "\\."), row.names = c("GVHD", "GroupID", "R_D"), stringsAsFactors = F)), stringsAsFactors = F)

Matched_pair_ID <- which(SampleID1$GroupID == SampleID2$GroupID)
Matched_pair_ID <- Matched_pair_ID[which(sapply(1:length(Matched_pair_ID), function(x) SampleID1$R_D[Matched_pair_ID[x]]!=SampleID2$R_D[Matched_pair_ID[x]]))]

## Samples from different group IDs
Random_pair_ID_idx <- which(SampleID1$GroupID != SampleID2$GroupID)
##  Make sure random samples are R-D pairs
Random_pair_ID <- Random_pair_ID_idx[which(sapply(1:length(Random_pair_ID_idx), function(x) SampleID1$R_D[Random_pair_ID_idx[x]]!=SampleID2$R_D[Random_pair_ID_idx[x]]))]

MHC_random_Percent <- data.frame(SampleID1 = Random_high_pert$SampleID1, 
                                 SampleID2 = Random_high_pert$SampleID2,
                                 HLA.A.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 HLA.B.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 HLA.C.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 HLA.DRB1.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 HLA.DQB1.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 MHC.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 ARS.IBDLength = numeric(length(Random_high_pert$SampleID2)),
                                 # IBDPercent = numeric(length(Random_high_pert$SampleID2)),
                                 stringsAsFactors = F)
MHC_matched_Percent <- data.frame(SampleID1 = Matched_high_pert$SampleID1, 
                                  SampleID2 = Matched_high_pert$SampleID2,
                                  HLA.A.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  HLA.B.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  HLA.C.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  HLA.DRB1.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  HLA.DQB1.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  MHC.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  ARS.IBDLength = numeric(length(Matched_high_pert$SampleID2)),
                                  stringsAsFactors = F)

# chrLength <- 133797422 # chr10 #170805979 # chr6
MHC_region_index <- data.frame(HLA = c("A", "B", "C", "DRB1", "DQB1"),
                               startID = c(29941260, 31353872, 31268749, 32578769, 32659467),
                               endID = c(29945884, 31357187, 31272086, 32589848, 32668383),
                               stringsAsFactors = F)
MHC_region_index <- MHC_region_index[order(MHC_region_index$startID), ]
MHC_region_index$Length <- MHC_region_index$endID - MHC_region_index$startID + 1

GVHD_IBD <- IBD_MHC_segment_matrix(MHC_IBD_table, RegionName = "MHC", sum(MHC_region_index$Length))

Matched_group_length <- IBD_MHC_segment_data_frame(MHC_IBD_table[Matched_pair_ID, ], RegionName = "MHC", sum(MHC_region_index$Length))
Random_group_length <- IBD_MHC_segment_data_frame(MHC_IBD_table[Random_pair_ID, ], RegionName = "MHC", sum(MHC_region_index$Length))

Matched_ARS_length <- IBD_MHC_segment_data_frame(MHC_IBD_table[Matched_pair_ID, ], RegionName = "ARS", sum(MHC_region_index$Length))
Random_ARS_length <- IBD_MHC_segment_data_frame(MHC_IBD_table[Random_pair_ID, ], RegionName = "ARS", sum(MHC_region_index$Length))

save(Matched_ARS_length, Random_ARS_length, Matched_group_length, Random_group_length,
     Matched_high_pert, Random_high_pert,
     file = "/mnt/cloudbiodata_nfs_2/users/hhuang/hli_renamed_all_chr6/IBD/ARS_MHC_chr6_WGS_IBD_0501.RData")

save(Matched_ARS_length, Random_ARS_length, Matched_group_length, Random_group_length,
     Matched_high_pert, Random_high_pert, all_matched_Percent, all_random_Percent,
     file = "../Output/IBDseq/summary/ARS_MHC_chr6_WGS_IBD_04062017.RData")


############ Plot
library(ggplot2)
load("../Output/IBDseq/summary/ARS_MHC_chr6_WGS_IBD_04062017.RData")
load("../Output/IBDseq/OriginalVCFIBD_chr6/ARS_MHC_chr6_WGS_IBD_0501.RData")
# load("../Output/IBDseq/summary/all_chromosome_stats_1030.RData")
# load(paste0(IBD_reformated_dir, "summaryIBDseq_summary_chr_6.RData"))

shared_DNA_theory <- data.frame(percent = c(0.5, 0.25, 0.125),
                                groups = c("parent-child/siblings", "uncle/aunt-niece/nephew", "1st Cousin"))

mc_length <- data.frame(Proportion = c(Matched_ARS_length$Percent, Random_ARS_length$Percent,     # ARS
                                       Matched_group_length$Percent, Random_group_length$Percent, # MHC
                                       Matched_high_pert$Percent, Random_high_pert$Percent), # Chr 6
                                       # all_matched_Percent$IBDPercent, all_random_Percent$IBDPercent), # WGS
                        Group = c(rep("Matched ARS", length(Matched_ARS_length$Percent)), rep("Random ARS", length(Random_ARS_length$Percent)),
                                  rep("Matched MHC", length(Matched_group_length$Percent)), rep("Random MHC", length(Random_group_length$Percent)),
                                  rep("Matched Chr6", length(Matched_high_pert$Percent)), rep("Random Chr6", length(Random_high_pert$Percent))),
                                  # rep("Matched WGS", length(all_matched_Percent$IBDPercent)), rep("Random WGS", length(all_random_Percent$IBDPercent))),
                        stringsAsFactors = F)

mc_length$Group <- factor(mc_length$Group, 
                          levels = c("Matched ARS", "Random ARS", 
                                     "Matched MHC", "Random MHC", 
                                     "Matched Chr6", "Random Chr6", 
                                     "Matched WGS", "Random WGS"))

mainplot <- ggplot(mc_length, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot(outlier.colour = "gray", outlier.size = 0) +
  guides(fill=FALSE)  +
  geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous("Percentage of IBD segments between pairs")
  # coord_flip() +
  # ggtitle(paste0("Length (%) of IBD segments \non Chromosome ", chr))
library(grid)

mc_length_chr6only <- data.frame(Proportion = c(Matched_high_pert$Percent, Random_high_pert$Percent),  # Chr 6
                                 Group = c(rep("Matched Chr6", length(Matched_high_pert$Percent)), rep("Random Chr6", length(Random_high_pert$Percent))),
                                 stringsAsFactors = F)
subplot_chr6 <- ggplot(mc_length_chr6only, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot() +
  guides(fill=FALSE) +
  # geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(element_blank()) + scale_x_discrete(element_blank()) 

vp_chr6 <- viewport(width = 0.2, height = 0.35, 
                   x = 0.73,
                   y = unit(0.75, "npc"), just = "right")

mc_length_WGSonly <- data.frame(Proportion = c(all_matched_Percent$IBDPercent, all_random_Percent$IBDPercent), # WGS
                        Group = c(rep("Matched WGS", length(all_matched_Percent$IBDPercent)), rep("Random WGS", length(all_random_Percent$IBDPercent))),
                        stringsAsFactors = F)

subplot_wgs <- ggplot(mc_length_WGSonly, aes(x = Group , y = Proportion, fill = Group)) + 
  geom_boxplot() +
  guides(fill=FALSE) +
  # geom_hline(yintercept = shared_DNA_theory$percent, color = "chocolate1", linetype="dashed") +
  scale_y_continuous(element_blank()) + scale_x_discrete(element_blank()) 

vp_wgs <- viewport(width = 0.2, height = 0.35, 
               x = 0.95,
               y = unit(0.75, "npc"), just = "right")

print(mainplot)

print(subplot_chr6, vp = vp_chr6)
print(subplot_wgs, vp = vp_wgs)
# theme_set(theme_bw())
