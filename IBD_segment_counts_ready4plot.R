# process HBD/IBD result files from BEAGLE
source('util.R', echo = FALSE)

# IBD_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/IBD_output/IBD_seq_output/"
# IBD_reformat_dir <- "../Output/IBDseq/cloud_Rformat/"
# IBD_reformat_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/"
IBD_segment_count_output_dir <-  "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/segment_counts/"

IBD_files <- list.files(IBD_segment_count_output_dir, pattern = "\\.RData$")
num_files <- length(IBD_files)

Matched_all_IBD <- data.frame(CHROM = character(0),
                              POS = numeric(0),
                              NumPairs = numeric(0),
                              Proportion = numeric(0),
                              stringsAsFactors = F)
Random_all_IBD <- Matched_all_IBD
aGVHD_all_IBD <- Matched_all_IBD
nGVHD_all_IBD <- Matched_all_IBD

for(fid in 1:num_files){
  
  load(paste0(IBD_segment_count_output_dir, IBD_files[fid]))
  
  Matched_all_IBD <- rbind(Matched_all_IBD, MatchedPair_IBD_dataframe)
  Random_all_IBD <- rbind(Random_all_IBD, RandomPair_IBD_dataframe)
  aGVHD_all_IBD <- rbind(aGVHD_all_IBD, aGVHD_IBD_dataframe)
  nGVHD_all_IBD <- rbind(nGVHD_all_IBD, nGVHD_IBD_dataframe)
 
  rm(MatchedPair_IBD_dataframe, RandomPair_IBD_dataframe, aGVHD_IBD_dataframe, nGVHD_IBD_dataframe) 
}

save(Matched_all_IBD, Random_all_IBD, file = paste0(IBD_segment_count_output_dir, "Matched_v_RandomPairs_IBD.RData"))
save(aGVHD_all_IBD, nGVHD_all_IBD, file = paste0(IBD_segment_count_output_dir, "aGVHD_v_nGVHD_IBD.RData"))

