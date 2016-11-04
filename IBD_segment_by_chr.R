# process HBD/IBD result files from BEAGLE
source('util.R', echo = FALSE)

load("../Data/ID_table.RData")
# load("../Data/GRCh38_gene_list.RData")

# IBD_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/IBD_output/IBD_seq_output/"
# IBD_reformat_dir <- "../Output/IBDseq/cloud_Rformat/"
IBD_reformat_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/"
IBD_segment_count_output_dir <-  "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/segment_counts/"

# IBD_files <- list.files(IBD_reformat_dir, pattern = "\\.RData$")
# num_files <- length(IBD_files)
IBD_files <- "ibdseq_output_all_chr11.gt.ibd.RData"
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


# for(fid in 1:num_files){
 fid <- 1
  ptm <- proc.time()
  
  # IBD_table <- read.table(file = paste0(IBD_file_dir, IBD_files[fid]))
  # colnames(IBD_table) <- c("SampleID1", "HapID1", "SampleID2", "HapID2", "Chr", "StartInd", "EndInd", "LOD")
  
  load(paste0(IBD_reformat_dir, IBD_files[fid]))
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
  
  chr <- as.numeric(gsub("chr", "", new_IBD_table$Chr[1]))
  MatchedPair_IBD_segment_counter <- vector(mode = "numeric", length = chrLengths[chr])
  RandomPair_IBD_segment_counter <- vector(mode = "numeric", length = chrLengths[chr])
  aGVHD_IBD_segment_counter <- vector(mode = "numeric", length = chrLengths[chr])
  nGVHD_IBD_segment_counter <- vector(mode = "numeric", length = chrLengths[chr])
  
  num_rows <- dim(new_IBD_table)[1]
  
  for (id in 1: num_rows){
    
    # IBD_segment_counter[IBD_table$StartInd[id]:IBD_table$EndInd[id]] <- IBD_segment_counter[IBD_table$StartInd[id]:IBD_table$EndInd[id]] + 1
    Group_SampleID1 <- unlist(strsplit(new_IBD_table$SampleID1[id], "\\."))
    Group_SampleID2 <- unlist(strsplit(new_IBD_table$SampleID2[id], "\\."))
    
    if(paste0(Group_SampleID1[1:2], collapse = "") == paste0(Group_SampleID2[1:2], collapse = "")){ ## Matched pairs

      MatchedPair_IBD_segment_counter[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- MatchedPair_IBD_segment_counter[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
      
      if(Group_SampleID1[1] == "a"){ # aGVHD pairs
        
        aGVHD_IBD_segment_counter[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- aGVHD_IBD_segment_counter[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
        
      }else{ # non-aGVHD pairs
        
        nGVHD_IBD_segment_counter[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- nGVHD_IBD_segment_counter[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
        
      }
      
    }else { ## random Pairs
      
      RandomPair_IBD_segment_counter[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] <- RandomPair_IBD_segment_counter[new_IBD_table$StartID[id]:new_IBD_table$EndID[id]] + 1
      
    }
    
  }
  
  MatchedPair_IBD_dataframe <- IBD_count_2_data_frame(MatchedPair_IBD_segment_counter, new_IBD_table$Chr[1])
  RandomPair_IBD_dataframe <- IBD_count_2_data_frame(RandomPair_IBD_segment_counter, new_IBD_table$Chr[1])
  aGVHD_IBD_dataframe <- IBD_count_2_data_frame(aGVHD_IBD_segment_counter, new_IBD_table$Chr[1])
  nGVHD_IBD_dataframe <- IBD_count_2_data_frame(nGVHD_IBD_segment_counter, new_IBD_table$Chr[1])

  file.out.RData <- gsub("_IBD.RData", "_segment_count.RData", IBD_files[fid])
  save(MatchedPair_IBD_dataframe, RandomPair_IBD_dataframe, 
       aGVHD_IBD_dataframe, nGVHD_IBD_dataframe, 
       file = paste0(IBD_segment_count_output_dir, file.out.RData))
  
  cat(paste0("IBD table for ", IBD_files[fid]," has been counted and saved!\n"))
  print(proc.time() - ptm)
  cat("\n")
# }

