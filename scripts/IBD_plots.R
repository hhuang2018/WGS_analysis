##
library(Matrix)
library(ggplot2)
#library(reshape)
#IBD_count_fp <- '/Users/hhuang2 (Deleted)/Documents/NGSProject/2018WGS/Data/HLI/segment_stats_winsize_1000/'
IBD_count_fp <- '/home/hhuang/data/IBD/segment_stats_winsize_1000/'

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
#### 
# pdf(paste0(IBD_count_fp,'Plots/IBD_by_pairs.pdf'), onefile = T, width = 1500, height = 3.4, res=300)  # ppi 300
#### 
# chr <- 1
for(chr in 1:22){
  # Matched
  load(paste0(IBD_count_fp, 'MatchedPair_IBD_chr', chr,'.RData')) # MatchedPair_IBD and MatchedPairIDs  
  matched_points <- length(which(MatchedPair_IBD ==1))
  
  # Random donor-donor pairs
  load(paste0(IBD_count_fp, 'Random_DD_IBD_chr', chr,'.RData')) # Random_DD_IBD  and Random_DD_IDs  
  
  #Random_DD_points <- length(which(Random_DD_IBD == 1))
  
  if(chr == 1){
    RandomDD_205_pairs_index <- sample(1:length(Random_DD_IDs), 205, replace = F) # random donor-donor pair
    RandomDD_205_points <- length(which(Random_DD_IBD[RandomDD_205_pairs_index,] == 1))
    RandomDD_205_pairs <- Random_DD_IDs[RandomDD_205_pairs_index]
    save(RandomDD_205_pairs, file = paste0(IBD_count_fp, 'Random_DD_205paris.RData'))
  }else{
    load(paste0(IBD_count_fp, 'Random_DD_205paris.RData'))
    RandomDD_205_pairs_index <- which(Random_DD_IDs %in% RandomDD_205_pairs)
    RandomDD_205_points <- length(which(Random_DD_IBD[RandomDD_205_pairs_index,] == 1))
  }
  
  # Random donor-recipient pairs
  load(paste0(IBD_count_fp, 'Random_RD_IBD_chr', chr,'.RData')) # Random_RD_IBD  and Random_RD_IDs  
  
  #Random_RD_points <- length(which(Random_RD_IBD == 1))
  if(chr == 1){
    RandomRD_205_pairs_id <- sample(1:length(Random_RD_IDs), 205, replace = F) # random donor-recipient pair
    RandomRD_205_points <- length(which(Random_RD_IBD[RandomRD_205_pairs_id,] == 1))
    RandomRD_205_pairs <- Random_RD_IDs[RandomRD_205_pairs_id]
    save(RandomDD_205_pairs, file = paste0(IBD_count_fp, 'Random_RD_205paris.RData'))
  }else{
    load(paste0(IBD_count_fp, 'Random_RD_205paris.RData'))
    RandomRD_205_pairs_id <- which(Random_RD_IDs %in% RandomRD_205_pairs)
    RandomRD_205_points <- length(which(Random_RD_IBD[RandomRD_205_pairs_id,] == 1))
  }
  
  # Random recipient-recipient pairs
  load(paste0(IBD_count_fp, 'Random_RR_IBD_chr', chr,'.RData')) # Random_RR_IBD  and Random_RR_IDs  
  
  if(chr==1){
    RandomRR_205_pairs_id <- sample(1:length(Random_RR_IDs), 205, replace = F) # random recipient-recipient pair
    RandomRR_205_points <- length(which(Random_RR_IBD[RandomRR_205_pairs_id,] == 1))
    RandomRR_205_pairs <- Random_RR_IDs[RandomRR_205_pairs_id]
    save(RandomRR_205_pairs, file = paste0(IBD_count_fp, 'Random_RR_205paris.RData'))
  }else{
    load(paste0(IBD_count_fp, 'Random_RR_205paris.RData'))
    RandomRR_205_pairs_id <- which(Random_RR_IDs %in% RandomRR_205_pairs)
    RandomRR_205_points <- length(which(Random_RR_IBD[RandomRR_205_pairs_id,] == 1))
  }
  
  ###### Total
  total_points <- matched_points + RandomDD_205_points + RandomRD_205_points + RandomRR_205_points
  
  #total_points <- length(which(MatchedPair_IBD ==1))
  
  nrow <- dim(MatchedPair_IBD)[1] 
  IBD_dataFrame <- data.frame(PairID = integer(total_points),
                              POS = integer(total_points),
                              Type = character(total_points),
                              stringsAsFactors = F)
  
  counter <- 0
  for(id in 1:nrow){
    
    ## Matched
    POS <- which(MatchedPair_IBD[id,] == 1)
    
    if (length(POS)>0){
      index <- (counter+1):(counter+length(POS))
      
      IBD_dataFrame$PairID[index] <- id
      IBD_dataFrame$POS[index] <- POS
      
      temp_type <- unlist(strsplit(MatchedPairIDs[id], '\\.'))[1]
      
      if(temp_type == 'a'){
        Type <- 'aGVHD'
      }else{
        Type <- 'nonGVHD'
      }
      
      IBD_dataFrame$Type[index] <- Type
      
      counter <- counter+length(POS)
    }
    ### Random donor-donor pairs nrow*1 + id  
    POS <- which(Random_DD_IBD[RandomDD_205_pairs_index[id],] == 1)
    
    if (length(POS)>0){
      index <- (counter+1):(counter+length(POS))
      
      IBD_dataFrame$PairID[index] <- id + nrow * 1
      IBD_dataFrame$POS[index] <- POS
      
      IBD_dataFrame$Type[index] <- 'Random D-D'
      
      counter <- counter+length(POS)
    }
    
    ### Random donor-recipient pairs: nrow * 2 + id  
    POS <- which(Random_RD_IBD[RandomRD_205_pairs_id[id],] == 1)
    
    if (length(POS)>0){
      index <- (counter+1):(counter+length(POS))
      
      IBD_dataFrame$PairID[index] <- id + nrow * 2
      IBD_dataFrame$POS[index] <- POS
      
      IBD_dataFrame$Type[index] <- 'Random D-R'
      
      counter <- counter+length(POS)
    }
    ### Random recipient-recipient pairs: nrow * 3 + id  
    POS <- which(Random_RR_IBD[RandomRR_205_pairs_id[id],] == 1)
    
    if (length(POS)>0){
      index <- (counter+1):(counter+length(POS))
      
      IBD_dataFrame$PairID[index] <- id + nrow * 3
      IBD_dataFrame$POS[index] <- POS
      
      IBD_dataFrame$Type[index] <- 'Random R-R'
      
      counter <- counter+length(POS)
    }
  }
  
  Levels <- c('aGVHD', 'nonGVHD', 'Random D-R', 'Random D-D', 'Random R-R')
  IBD_dataFrame$Type <- factor(IBD_dataFrame$Type, levels = Levels)
  
  # plot(IBD_dataFrame$POS, IBD_dataFrame$Group)
  p1 <- ggplot(data = IBD_dataFrame, aes(x = POS, y = PairID, color=Type)) + 
    #xlim(0, chrLengths[chr]) + 
    geom_point(alpha = 1/100) +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) + # + geom_rug()
    ggtitle(paste0('aGVHD vs nonGVHD vs Random 205 pairs - Chr ', chr))
  #object.size(p1)
  #
  #ggsave(filename = paste0(IBD_count_fp,'Plots/Matched_chr',chr,'.png'), plot = p1, width = 6, height = 3, dpi = 300)
  png(filename = paste0(IBD_count_fp,'Plots/plots_IBD_by_pairs_chr_',chr,'.png'),
      width = 1500, height = 1000, units = "px")
  print(p1)
  save(p1, file = paste0(IBD_count_fp, 'Plots/plots_IBD_by_pairs_chr_', chr, '.RData'))
  dev.off()
  
  rm(p1)
}
