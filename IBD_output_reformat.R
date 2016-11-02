# process HBD/IBD result files from BEAGLE
source('util.R', echo = FALSE)

load("../Data/ID_table.RData")
load("../Data/GRCh38_gene_list.RData")

IBD_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/IBD_output/IBD_seq_output/"
IBD_reformat_output <- "/mnt/cloudbiodata_nfs_2/users/hhuang/IBD/IBD_seq_output/"

IBD_files <- list.files(IBD_file_dir, pattern = "\\.gt.ibd$")
num_files <- length(IBD_files)

for(fid in 1:num_files){
  ptm <- proc.time() 
  
  IBD_table <- read.table(file = paste0(IBD_file_dir, IBD_files[fid]))
  
  colnames(IBD_table) <- c("SampleID1", "HapID1", "SampleID2", "HapID2", "Chr", "StartInd", "EndInd", "LOD")
  
  num_rows <- dim(IBD_table)[1]
  
  new_IBD_table <- data.frame(SampleID1 = character(num_rows),
                              HapID1 = numeric(num_rows),
                              SampleID2 = character(num_rows),
                              HapID2 = numeric(num_rows),
                              Chr = character(num_rows),
                              StartID = numeric(num_rows),
                              EndID = numeric(num_rows),
                              LOD = numeric(num_rows),
                              GeneNames = character(num_rows),
                              stringsAsFactors = FALSE)
  
  for (id in 1: num_rows){
    
    SampleID <- paste0(ID_table[which(ID_table$SeqID == IBD_table$SampleID1[id]), 1:3], collapse = ".")
    new_IBD_table$SampleID1[id] <- SampleID
    new_IBD_table$HapID1[id] <- IBD_table$HapID1[id]
    
    SampleID <- paste0(ID_table[which(ID_table$SeqID == IBD_table$SampleID2[id]), 1:3], collapse = ".")
    new_IBD_table$SampleID2[id] <- SampleID
    new_IBD_table$HapID2[id] <- IBD_table$HapID2[id]
    
    new_IBD_table$Chr[id] <- as.character(IBD_table$Chr[id])
    new_IBD_table$StartID[id] <- IBD_table$StartInd[id]
    new_IBD_table$EndID[id] <- IBD_table$EndInd[id]
    new_IBD_table$LOD[id] <- IBD_table$LOD[id]
    
    new_IBD_table$GeneNames[id] <- get.GeneNames(chrom = new_IBD_table$Chr[id], 
                                                 startPos = new_IBD_table$StartID[id],
                                                 endPos = new_IBD_table$EndID[id],
                                                 GeneList = GRCh38_gene_list)
    
    
    
  }
  
  file.out.RData <- gsub(".gt.ibd", "_IBD.RData", IBD_files[fid])
  file.out.csv <- gsub(".gt.ibd", "_IBD.csv", IBD_files[fid])
  save(new_IBD_table, file = paste0(IBD_reformat_output, file.out.RData))
  write.csv(new_IBD_table, file = paste0(IBD_reformat_output, file.out.csv), row.names = FALSE)
  
  cat(paste0("IBD table for ", IBD_files[fid]," has been converted and saved!\n"))
  print(proc.time() - ptm)
  cat("\n")
}


###################
# HBD
###################
# rm(list = ls())
IBD_files <- list.files(IBD_file_dir, pattern = "\\.gt.hbd$")
num_files <- length(HBD_files)

for(fid in 1:num_files){
  ptm <- proc.time() 
  IBD_table <- read.table(file = paste0(IBD_file_dir, IBD_files[fid]))
  colnames(IBD_table) <- c("SampleID1", "HapID1", "SampleID2", "HapID2", "Chr", "StartInd", "EndInd", "LOD")
  
  num_rows <- dim(IBD_table)[1]
  
  new_IBD_table <- data.frame(SampleID1 = character(num_rows),
                              HapID1 = numeric(num_rows),
                              SampleID2 = character(num_rows),
                              HapID2 = numeric(num_rows),
                              Chr = character(num_rows),
                              StartID = numeric(num_rows),
                              EndID = numeric(num_rows),
                              LOD = numeric(num_rows),
                              GeneNames = character(num_rows),
                              stringsAsFactors = FALSE)
  
  for (id in 1: num_rows){
    
    SampleID <- paste0(ID_table[which(ID_table$SeqID == IBD_table$SampleID1[id]), 1:3], collapse = ".")
    new_IBD_table$SampleID1[id] <- SampleID
    new_IBD_table$HapID1[id] <- IBD_table$HapID1[id]
    
    SampleID <- paste0(ID_table[which(ID_table$SeqID == IBD_table$SampleID2[id]), 1:3], collapse = ".")
    new_IBD_table$SampleID2[id] <- SampleID
    new_IBD_table$HapID2[id] <- IBD_table$HapID2[id]
    
    new_IBD_table$Chr[id] <- as.character(IBD_table$Chr[id])
    new_IBD_table$StartID[id] <- IBD_table$StartInd[id]
    new_IBD_table$EndID[id] <- IBD_table$EndInd[id]
    new_IBD_table$LOD[id] <- IBD_table$LOD[id]
    
    new_IBD_table$GeneNames[id] <- get.GeneNames(chrom = new_IBD_table$Chr[id], 
                                                 startPos = new_IBD_table$StartID[id],
                                                 endPos = new_IBD_table$EndID[id],
                                                 GeneList = GRCh38_gene_list)
    
  }
  
  file.out.RData <- gsub(".gt.hbd", "_HBD.RData", IBD_files[fid])
  file.out.csv <- gsub(".gt.hbd", "_HBD.csv", IBD_files[fid])
  save(new_IBD_table, file = paste0(IBD_reformat_output, file.out.RData))
  write.csv(new_IBD_table, file = paste0(IBD_reformat_output, file.out.csv), row.names = FALSE)
  
  cat(paste0("HBD table for ", IBD_files[fid]," has been converted and saved!\n"))
  print(proc.time() - ptm)
  cat("\n")
}