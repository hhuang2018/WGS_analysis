source("utils/util.R")

############ 
IBD_reformated_dir <- "../Output/IBDseq/R_reformated//"

load(paste0(IBD_reformated_dir, "ibdseq_output_all_chr6_IBD.RData"))

total_sample_IDs <- unique(new_IBD_table$SampleID1)

aGVHD_samples <- total_sample_IDs[grep("a", total_sample_IDs)]

aGVHD_donor <- aGVHD_samples[grep(".D", aGVHD_samples)]     # n = 102
aGVHD_recipient <- aGVHD_samples[grep(".R", aGVHD_samples)] # n = 102

nGVHD_samples <-  total_sample_IDs[grep("n.", total_sample_IDs)]

nGVHD_donor <-  total_sample_IDs[grep(".D", nGVHD_samples)]     # n = 102
nGVHD_recipient <-  total_sample_IDs[grep(".R", nGVHD_samples)] # n = 103

#############

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
chr_len <- 170805979

IBD_table <- matrix(0, nrow = 4, ncol = chr_len)
IBD_table <- as.data.frame(IBD_table)

num_pairs <- dim(new_IBD_table)[1]

matched_DR_pairs <- vector("character", length = 0)
matched_DR_pairs_counter <- 0

random_DR_pairs <- vector("character", length = 0)
random_DR_pairs_counter <- 0

random_DD_pairs <- vector("character", length = 0)
random_DD_pairs_counter <- 0

random_RR_pairs <- vector("character", length = 0)
random_RR_pairs_counter <- 0

for(id in 1:num_pairs){
  
  group1 <- paste0(strsplit(new_IBD_table$SampleID1[id],'\\.')[[1]][1:2], collapse = '.')
  group2 <- paste0(strsplit(new_IBD_table$SampleID2[id],'\\.')[[1]][1:2], collapse = '.')
  
  if(group1 == group2){ # matching donor-recipient pairs - row 1
    IBD_table[1, new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- IBD_table[1, new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
    if(length(matched_DR_pairs)==0){
      matched_DR_pairs <- c(matched_DR_pairs, group1)
      matched_DR_pairs_counter <- matched_DR_pairs_counter + 1
    }else if(!(group1 %in% matched_DR_pairs)){
      
      matched_DR_pairs <- c(matched_DR_pairs, group1)
      matched_DR_pairs_counter <- matched_DR_pairs_counter + 1
      
    }
  }else{ # random pairs
    
    groupIDs <- paste0(sort(c(group1, group2)), collapse = '-')
    
    type1 <- strsplit(new_IBD_table$SampleID1[id],'\\.')[[1]][3]
    type2 <- strsplit(new_IBD_table$SampleID2[id],'\\.')[[1]][3]
    
    if(type1 != type2){ # random donor-recipient pars  - row 2
      IBD_table[2, new_IBD_table$StartID:new_IBD_table$EndID] <- IBD_table[2, new_IBD_table$StartID:new_IBD_table$EndID] + 1
      
      if(length(random_DR_pairs)==0){
        random_DR_pairs <- c(random_DR_pairs, groupIDs)
        random_DR_pairs_counter <- random_DR_pairs_counter + 1
      }else if(!(groupIDs %in% random_DR_pairs)){
        random_DR_pairs <- c(random_DR_pairs, groupIDs)
        random_DR_pairs_counter <- random_DR_pairs_counter + 1
      }
      
    }else{
      
      if(type1 == 'D'){  # random donor-donor - row 3
        IBD_table[3, new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- IBD_table[3, new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
        
        if(length(random_DD_pairs)==0){
          random_DD_pairs <- c(random_DD_pairs, groupIDs)
          random_DD_pairs_counter <- random_DD_pairs_counter + 1
        }else if(!(groupIDs %in% random_DR_pairs)){
          random_DD_pairs <- c(random_DD_pairs, groupIDs)
          random_DD_pairs_counter <- random_DD_pairs_counter + 1
        }
        
      }else{ # random recipient-recipient pairs - row 4
        
        IBD_table[4, new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- IBD_table[4, new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
        
        if(length(random_RR_pairs)==0){
          random_RR_pairs <- c(random_RR_pairs, groupIDs)
          random_RR_pairs_counter <- random_RR_pairs_counter + 1
        }else if(!(groupIDs %in% random_RR_pairs)){
          random_RR_pairs <- c(random_RR_pairs, groupIDs)
          random_RR_pairs_counter <- random_RR_pairs_counter + 1
        }
      }

    }
  }
    
}
  
#plot(IBD_table)
