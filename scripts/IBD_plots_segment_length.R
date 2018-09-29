##
library(Matrix)
library(ggplot2)
#library(reshape)
IBD_count_fp <- '/Users/hhuang2 (Deleted)/Documents/NGSProject/2018WGS/Data/HLI/segment_stats_winsize_1000/'
# IBD_count_fp <- '/home/hhuang/data/IBD/segment_stats_winsize_1000/'

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
Matched_pair_n <- 205

for(chr in 1:22){
  # Matched
  load(paste0(IBD_count_fp, 'MatchedPair_IBD_chr', chr,'.RData')) # MatchedPair_IBD and MatchedPairIDs 
  
  # matched_points <- length(which(MatchedPair_IBD ==1))
  matched_IBD_perc <- colSums(MatchedPair_IBD)/Matched_pair_n *100
  
  GVHD_type <- sapply(1:length(MatchedPairIDs), function(x) unlist(strsplit(MatchedPairIDs[x], '\\.'))[1])
  aGVHD_index <- which(GVHD_type == "a")
  nonGVHD_index <- which(GVHD_type == "n")
  
  aGVHD_IBD_perc <- colSums(MatchedPair_IBD[aGVHD_index, ])/length(aGVHD_index) * 100
  nonGVHD_IBD_perc <- colSums(MatchedPair_IBD[nonGVHD_index, ])/length(nonGVHD_index) * 100
  
  # Random donor-donor pairs
  load(paste0(IBD_count_fp, 'Random_DD_IBD_chr', chr,'.RData')) # Random_DD_IBD  and Random_DD_IDs  
  
  Random_DD_pair_n <- length(Random_DD_IDs)
  Random_DD_IBD_perc <- colSums(Random_DD_IBD)/Random_DD_pair_n *100

  # Random donor-recipient pairs
  load(paste0(IBD_count_fp, 'Random_RD_IBD_chr', chr,'.RData')) # Random_RD_IBD  and Random_RD_IDs  
  
  Random_RD_pair_n <- length(Random_RD_IDs)
  Random_RD_IBD_perc <- colSums(Random_RD_IBD)/Random_RD_pair_n *100

  # Random recipient-recipient pairs
  load(paste0(IBD_count_fp, 'Random_RR_IBD_chr', chr,'.RData')) # Random_RR_IBD  and Random_RR_IDs  
  
  Random_RR_pair_n <- length(Random_RR_IDs)
  Random_RR_IBD_perc <- colSums(Random_RR_IBD)/Random_RR_pair_n *100
  
  ###### Total
  #total_points <- length(which(MatchedPair_IBD ==1))
  num_POS <- length(Random_RD_IBD_perc)
  total_points <- num_POS * 6
  POS <- 1:num_POS
  IBD_dataFrame <- data.frame(POS = integer(total_points),
                              Percent = double(total_points),
                              Type = character(total_points),
                              stringsAsFactors = F)
  # aGVHD
  INDEX <- 1:num_POS
  IBD_dataFrame$POS[INDEX] <- POS
  IBD_dataFrame$Percent[INDEX] <- aGVHD_IBD_perc
  IBD_dataFrame$Type[INDEX] <- 'aGVHD'
  
  # nonGVHD
  INDEX <- (num_POS+1):(num_POS*2)
  IBD_dataFrame$POS[INDEX] <- POS
  IBD_dataFrame$Percent[INDEX] <- nonGVHD_IBD_perc 
  IBD_dataFrame$Type[INDEX] <- 'nonGVHD'
  
  # Matched Pairs
  INDEX <- (num_POS*2+1):(num_POS*3)
  IBD_dataFrame$POS[INDEX] <- POS
  IBD_dataFrame$Percent[INDEX] <- matched_IBD_perc 
  IBD_dataFrame$Type[INDEX] <- 'Matched Pairs'
  
  # Random D-R
  INDEX <- (num_POS*3+1):(num_POS*4)
  IBD_dataFrame$POS[INDEX] <- POS
  IBD_dataFrame$Percent[INDEX] <- Random_RD_IBD_perc 
  IBD_dataFrame$Type[INDEX] <- 'Random D-R'
  
  # Random D-D
  INDEX <- (num_POS*4+1):(num_POS*5)
  IBD_dataFrame$POS[INDEX] <- POS
  IBD_dataFrame$Percent[INDEX] <- Random_DD_IBD_perc 
  IBD_dataFrame$Type[INDEX] <- 'Random D-D'
  
  # Random R-R
  INDEX <- (num_POS*5+1):(num_POS*6)
  IBD_dataFrame$POS[INDEX] <- POS
  IBD_dataFrame$Percent[INDEX] <- Random_RR_IBD_perc 
  IBD_dataFrame$Type[INDEX] <- 'Random R-R'
  
  Levels <- c('aGVHD', 'nonGVHD', 'Matched Pairs', 'Random D-R', 'Random D-D', 'Random R-R')
  IBD_dataFrame$Type <- factor(IBD_dataFrame$Type, levels = Levels)
  
  # plot(IBD_dataFrame$POS, IBD_dataFrame$Group)
  p1 <- ggplot(data = IBD_dataFrame, aes(x = POS, y = Percent, color=Type)) + 
    #xlim(0, chrLengths[chr]) + 
    geom_point(alpha = 1/50) + 
    guides(colour = guide_legend(override.aes = list(alpha = 1))) + # + geom_rug()
    ggtitle(paste0('Freq IBD seg - Chr ', chr)) +
    labs(y='Frequency (%)')
  #object.size(p1)
  #
  #ggsave(filename = paste0(IBD_count_fp,'Plots/Matched_chr',chr,'.png'), plot = p1, width = 6, height = 3, dpi = 300)
  png(filename = paste0(IBD_count_fp,'Plots/IBD_percents_by_pairs_chr_',chr,'.png'),
      width = 1500, height = 1000, units = "px")
  print(p1)
  save(IBD_dataFrame, p1, file = paste0(IBD_count_fp, 'Plots/IBD_segment_freq_chr_', chr, '.RData'))
  dev.off()
  
  rm(p1)
}
