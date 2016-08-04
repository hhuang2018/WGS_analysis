# process HBD/IBD result files from BEAGLE
source('util.R', echo = FALSE)

load("../Data/ID_table.RData")

IBD_table <- read.table(file = "../aGVHD_10_pairs_ref.gt.ibd.ibd")
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

load("../Data/GRCh38_gene_list.RData")

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

save(new_IBD_table, file = "../Output/aGVHD_10_pairs_IBD.RData")
write.csv(new_IBD_table, file = "../Output/aGVHD_10_pairs_IBD.csv", row.names = FALSE)
cat("IBD table has been converted and saved!")
###################
# HBD
###################
rm(list = ls())
source('util.R', echo = FALSE)

load("../Data/ID_table.RData")
load("../Data/GRCh38_gene_list.RData")

IBD_table <- read.table(file = "../aGVHD_10_pairs_ref.gt.ibd.hbd")
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

# load("../Data/GRCh38_gene_list.RData")

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

save(new_IBD_table, file = "../Output/aGVHD_10_pairs_HBD.RData")
write.csv(new_IBD_table, file = "../Output/aGVHD_10_pairs_HBD.csv", row.names = FALSE)
cat("HBD table has been converted and saved!")