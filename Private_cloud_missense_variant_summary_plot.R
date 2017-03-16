# library("ggbio")
# # load(system.file("data", "hg19IdeogramCyto.rda", package="biovizBase", mustWork=TRUE))
# library(BSgenome.Hsapiens.UCSC.hg38)
# # library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(GenomicRanges)
# # library(biovizBase)
# chr.len = seqlengths(BSgenome.Hsapiens.UCSC.hg38)  # get chromosome lengths
# # remove X,Y,M and random chromosomes
# chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]
load("../Data/hg38_known_gene_symbols_HGNC.RData")

Stats_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants_RefSeq_canonical/stats/"
# Stats_dir <- "../Output/Missense_variant_stats/RefSeq_canon_0228/"
load(file = paste0(Stats_dir, "Recipient_missesense_stats_RefSeq_canon_0228.RData"))
recipient_missense_stats <- donor_missense_stats
load(file = paste0(Stats_dir, "donor_missesense_stats_RefSeq_canon_0228.RData"))

# load("../Output/Missense_variant_stats/donor_missesense_stats_updated.RData")
# load("../Output/Missense_variant_stats/recipient_missesense_stats_updated.RData")

# write.csv(recipient_missense_stats, file = "../Output/Missense_variant_stats/all_recipient_missesense_stats_updated.csv", row.names = F)
# write.csv(donor_missense_stats, file = "../Output/Missense_variant_stats/all_donor_missesense_stats.csv_updated", row.names = F)

# load("../Output/Missense_variant_stats/Matched_donor_missesense_stats_updated.RData")
# load("../Output/Missense_variant_stats/Matched_recipient_missesense_stats_updated.RData")
# 
# write.csv(recipient_missense_stats, file = "../Output/Missense_variant_stats/Matched_recipient_missesense_stats_updated.csv", row.names = F)
# write.csv(donor_missense_stats, file = "../Output/Missense_variant_stats/Matched_donor_missesense_stats_updated.csv", row.names = F)

# recipient_missense_stats$NumDiff <- -recipient_missense_stats$NumDiff
recipient_missense_stats$PertDiff <- recipient_missense_stats$NumDiff/240
donor_missense_stats$PertDiff <- donor_missense_stats$NumDiff/216
# missense_all_cases <- rbind(donor_missense_stats, recipient_missense_stats)
# missense_all_cases$group <- c(rep("donor", dim(donor_missense_stats)[1]), 
#                               rep("recipient", dim(recipient_missense_stats)[1]))
# missense_all_cases$POS2 <- as.numeric(missense_all_cases$POS)
# missense_all_cases$POS <- as.numeric(missense_all_cases$POS)
# GR_missense_all_cases <- makeGRangesFromDataFrame(missense_all_cases, keep.extra.columns = T,
#                                                   ignore.strand = T, seqnames = "CHROM",
#                                                   start.field = "POS", end.field = "POS")
# seqlengths(GR_missense_all_cases) <- chr.len[names(seqlengths(GR_missense_all_cases))]


#### Sort donor v patient 
num_POS <- dim(unique(rbind(recipient_missense_stats[, c(1,2)], donor_missense_stats[, c(1,2)])))[1]
Missense_stats_LOD <- data.frame(CHROM = character(num_POS),
                                 POS = numeric(num_POS),
                                 LOD = numeric(num_POS),
                                 stringsAsFactors = F)

counter <- 0
for(chr_id in 1:22){
  CHROM <- paste0("chr", chr_id)
  CHR_recipient_missense <- recipient_missense_stats[which(recipient_missense_stats$CHROM %in% CHROM),]
  CHR_donor_missense <- donor_missense_stats[which(donor_missense_stats$CHROM %in% CHROM), ]
  all_POS <- union(CHR_donor_missense$POS, CHR_recipient_missense$POS)
  chr_num_POS <- length(all_POS)
  for(POS_id in 1:chr_num_POS){
    
    counter <- counter + 1
    Missense_stats_LOD$CHROM[counter] <- CHROM
    Missense_stats_LOD$POS[counter] <- POS_id
    
    recipient_index <- which(CHR_recipient_missense$POS %in% all_POS[POS_id])
    donor_index <- which(CHR_donor_missense$POS %in% all_POS[POS_id])
    
    if(length(recipient_index)*length(donor_index) != 0){
      
      Missense_stats_LOD$LOD[counter] <- log10(CHR_recipient_missense$NumDiff[recipient_index] / CHR_donor_missense$NumDiff[donor_index])
      
    }else if(length(donor_index) == 0 ){
      
      Missense_stats_LOD$LOD[counter] <- log10(CHR_recipient_missense$NumDiff[recipient_index] / 1e-8)
      
    }else if(length(recipient_index) == 0){
      
      Missense_stats_LOD$LOD[counter] <- log10(1e-8 / CHR_donor_missense$NumDiff[donor_index])
      
    }
    
  }
  
}

sorted_Missense_stats_LOD <- Missense_stats_LOD[order(Missense_stats_LOD$LOD, decreasing = T),]
sorted_Missense_stats_LOD$GeneRegion <- ""
sorted_Missense_stats_LOD$ProteinID <- ""
sorted_Missense_stats_LOD$alignID <- ""
sorted_Missense_stats_LOD$name <- ""
num_rows <- dim(sorted_Missense_stats_LOD)
for(id in 1:num_rows){
  
  chr_ids <- which(GRCh38_known_gene_list$Chromosome %in% sorted_Missense_stats_LOD$CHROM[id])
  
  IntervIndex <- sapply(chr_ids, 
                        function(x) findInterval(sorted_Missense_stats_LOD$POS[id], 
                                                 c(GRCh38_known_gene_list$txStart[x], GRCh38_known_gene_list$txEnd[x])))
  
  Within_index <- which(IntervIndex == 1)
  if(length(Within_index) > 0){
    
    sorted_Missense_stats_LOD$GeneRegion[id] <- GRCh38_known_gene_list$genesymbol[chr_ids[IntervIndex[Within_index]]]
    sorted_Missense_stats_LOD$ProteinID[id] <- GRCh38_known_gene_list$proteinID[chr_ids[IntervIndex[Within_index]]]
    sorted_Missense_stats_LOD$alignID[id] <- GRCh38_known_gene_list$alignID[chr_ids[IntervIndex[Within_index]]]
    sorted_Missense_stats_LOD$name[id] <- GRCh38_known_gene_list$name[chr_ids[IntervIndex[Within_index]]]
    
  }
  
}

save(sorted_Missense_stats_LOD, file = paste0(Stats_dir, "Sorted_Missense_genes_LODstats.RData"))