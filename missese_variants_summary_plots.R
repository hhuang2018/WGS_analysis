library("ggbio")
# load(system.file("data", "hg19IdeogramCyto.rda", package="biovizBase", mustWork=TRUE))
library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
# library(biovizBase)
chr.len = seqlengths(BSgenome.Hsapiens.UCSC.hg38)  # get chromosome lengths
# remove X,Y,M and random chromosomes
chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]

load("../Output/Missense_variant_stats/RefSeq_canon_0228/Recipient_missesense_stats_RefSeq_canon_0228.RData")
recipient_missense_stats <- donor_missense_stats
load("../Output/Missense_variant_stats/RefSeq_canon_0228/donor_missesense_stats_RefSeq_canon_0228.RData")

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
missense_all_cases <- rbind(donor_missense_stats, recipient_missense_stats)
missense_all_cases$group <- c(rep("donor", dim(donor_missense_stats)[1]), 
                              rep("recipient", dim(recipient_missense_stats)[1]))
missense_all_cases$POS2 <- as.numeric(missense_all_cases$POS)
missense_all_cases$POS <- as.numeric(missense_all_cases$POS)
GR_missense_all_cases <- makeGRangesFromDataFrame(missense_all_cases, keep.extra.columns = T,
                                                  ignore.strand = T, seqnames = "CHROM",
                                                  start.field = "POS", end.field = "POS")
seqlengths(GR_missense_all_cases) <- chr.len[names(seqlengths(GR_missense_all_cases))]


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

# known gene list
load("../Data/hg38_known_gene_symbols_HGNC.RData")
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

# hg38 cytoband ideogram
load("../Data/hg38_cytobandIdeogram.RData")
# MHC gene symbols and coordinates
load("../Data/hg38_MHC_region_gene_symbols_HGNC.RData")
# known gene list
load("../Data/hg38_known_gene_symbols_HGNC.RData")

# known MiHAs 
# known_MiHA_SNPs <- makeGRangesFromDataFrame(known_MiHA_table, keep.extra.columns = T,
#                                             ignore.strand = T, seqnames = "Chr",
#                                             start.field = "Pos", end.field = "Pos")
# seqlengths(known_MiHA_SNPs) <- chr.len[names(seqlengths(known_MiHA_SNPs))]
# pp <- autoplot(known_MiHA_SNPs, layout = "karyogram", geom = "rect") + scale_color_discrete("brown")## SNPs in the chromosome space

## ClinVar - AML related genes (31)
clinVar.AML <- read.delim("../Data/clinvar_result_acuteMyeloidLeukemia.txt")
AML_genes <- GRCh38_known_gene_list[grepl(paste0(unique(as.character(clinVar.AML$Gene.s.)), collapse = "|"), GRCh38_known_gene_list$genesymbol),]
AML_genes$Chromosome <- as.character(AML_genes$Chromosome)
AML_genes <- AML_genes[!grepl("_", AML_genes$Chromosome), ]
# GRCh38_known_gene_list$genesymbol[grepl("HLA-A|HLA-B|HLA-C|HLA-DRB|HLA-DQB", GRCh38_known_gene_list$genesymbol)]
GR_AML_genes <- makeGRangesFromDataFrame(AML_genes, seqnames.field = "Chromosome",
                                           start.field = "txStart", end.field = "txEnd", 
                                           strand.field = "strand",
                                           keep.extra.columns = T)
seqlengths(GR_AML_genes) <- chr.len[names(seqlengths(GR_AML_genes))]
save(GR_AML_genes, AML_genes, file = "../Data/ClinVar_AML_genes_acuteML.RData")

load("../Data/ClinVar_AML_genes.RData")

#######
hg38<- keepSeqlevels(hg38_cytoband, paste0("chr", c(1:22, "X", "Y")))  ## extract only chr1~22
pp <- autoplot(hg38, layout = "karyogram", cytoband = T) ## + scale_fill_giemsa()

pp + layout_karyogram(GR_AML_genes, layout = "karyogram", size = 1, color = "red") + #scale_color_discrete("red") +
  layout_karyogram(GR_MHC_regions, layout = "karyogram", size = .5, color = "#0072B2") +
  layout_karyogram(GR_missense_all_cases, aes(x = POS2, y = PertDiff, color = factor(group)),
                   ylim = c(10,40), rect.height = 10, geom="point", size = 0.3, alpha = 0.2) + 
  # scale_colour_manual(values = alpha(c("#D55E00", "#0072B2"), .1)) +
  # scale_fill_manual(values = alpha(c("#D55E00", "#0072B2"), .1)) +
  theme(#axis.title.x=element_blank(),
    legend.position = "none")
