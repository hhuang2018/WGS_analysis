# process HBD/IBD result files from BEAGLE
source('utils/util.R', echo = FALSE)

#load("../Data/ID_table.RData")

IBD_reformat_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/"
IBD_segment_count_output_dir <-  "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/segment_stats_plots/"

chr <- 4
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
  
df_IBD <- data.frame(POS = 1:chrLengths[chr],
                     MatchedPairs = vector(mode = "numeric", length = chrLengths[chr]),
                     RandomPairs = vector(mode = "numeric", length = chrLengths[chr]),
                     aGVHD = vector(mode = "numeric", length = chrLengths[chr]),
                     nGVHD = vector(mode = "numeric", length = chrLengths[chr]),
                     rand_DD = vector(mode = "numeric", length = chrLengths[chr]),
                     rand_RR = vector(mode = "numeric", length = chrLengths[chr]),
                     rand_DR = vector(mode = "numeric", length = chrLengths[chr])
                     )

num_rows <- dim(new_IBD_table)[1]

MatchedPairIDs <- vector(mode='character', length = 0)
RandomPairs <- vector(mode='character', length = 0)
aGVHDIDs <- vector(mode='character', length = 0)
nGVHDs <- vector(mode='character', length = 0)
MatchedPairIDs <- vector(mode='character', length = 0)
MatchedPairIDs <- vector(mode='character', length = 0)

for (id in 1: num_rows){
  
  # IBD_segment_counter[IBD_table$StartInd[id]:IBD_table$EndInd[id]] <- IBD_segment_counter[IBD_table$StartInd[id]:IBD_table$EndInd[id]] + 1
  Group_SampleID1 <- unlist(strsplit(new_IBD_table$SampleID1[id], "\\."))
  Group_SampleID2 <- unlist(strsplit(new_IBD_table$SampleID2[id], "\\."))
  
  if(paste0(Group_SampleID1[1:2], collapse = "") == paste0(Group_SampleID2[1:2], collapse = "")){ ## Matched pairs
    
    df_IBD$MatchedPairs[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- df_IBD$MatchedPairs[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
    
    if(Group_SampleID1[1] == "a"){ # aGVHD pairs
      
      df_IBD$aGVHD[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- df_IBD$aGVHD[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
      
    }else{ # non-aGVHD pairs
      
      df_IBD$nGVHD[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- df_IBD$nGVHD[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
      
    }
    
  }else { ## random Pairs
    
    df_IBD$RandomPairs[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- df_IBD$RandomPairs[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
    
    if(Group_SampleID1[3] == "D"){ # Donor 
      
      if(Group_SampleID2[3] == "D"){ # random DD
        
        df_IBD$rand_DD[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- df_IBD$rand_DD[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
    
      }else{ # Random DR
        
        df_IBD$rand_DR[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- df_IBD$rand_DR[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
    
      }
      
    }else{ # Recipient
      
      if(Group_SampleID2[3] == "D"){ # random DR
        
        df_IBD$rand_DR[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- df_IBD$rand_DR[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
    
      }else{ # Random RR
        
        df_IBD$rand_RR[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- df_IBD$rand_RR[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
    
      }      
    }

  }
  
}

file.out.RData <- gsub("_IBD.RData", "_segment_distribution.RData", IBD_files)
save(df_IBD, 
     file = paste0(IBD_segment_count_output_dir, file.out.RData))

cat(paste0("IBD table for ", IBD_files," has been counted and saved!\n"))
print(proc.time() - ptm)
cat("\n")

# plot
#library(ggplot2)
library(reshape)
file.out.reshaped.RData <- gsub("_IBD.RData", "_stats_reshaped_plot_ready.RData", IBD_files)
df_IBD_reshape <- melt(df_IBD, id='POS', variable_name = 'pair_counts')
save(df_IBD_reshape, file = paste0(IBD_segment_count_output_dir, file.out.RData))

cat(paste0("Reshaped IBD table for ", IBD_files," has been created and saved!\n"))
cat('----------------------------------')
