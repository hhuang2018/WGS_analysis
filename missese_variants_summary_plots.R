library("ggbio")
# load(system.file("data", "hg19IdeogramCyto.rda", package="biovizBase", mustWork=TRUE))
library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
# library(biovizBase)
chr.len = seqlengths(BSgenome.Hsapiens.UCSC.hg38)  # get chromosome lengths
# remove X,Y,M and random chromosomes
chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]

# load("../Output/Missense_variant_stats/RefSeq_canon_0228/Recipient_missesense_stats_RefSeq_canon_0228.RData")
# load("/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants_wwVersion/AMLPatients_missesense_stats_wwVersion.RData")
load("../ClinVar/Missense_mutation_summary_wwVersion/AMLPatients_missesense_stats_wwVersion.RData")
recipient_missense_stats <- donor_missense_stats
# load("../Output/Missense_variant_stats/RefSeq_canon_0228/donor_missesense_stats_RefSeq_canon_0228.RData")
# load("/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants_wwVersion/Donor_missesense_stats_wwVersion.RData")
load("../ClinVar/Missense_mutation_summary_wwVersion/Donor_missesense_stats_wwVersion.RData")

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
recipient_missense_stats$PertDiff <- -recipient_missense_stats$NumDiff/196 ## AML patients: 196; all patients: 240
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
# load("../Data/hg38_known_gene_symbols_HGNC.RData")
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
GRCh38_known_gene_list$name <- as.character(GRCh38_known_gene_list$name)
GRCh38_known_gene_list$Chromosome <- as.character(GRCh38_known_gene_list$Chromosome)
GRCh38_known_gene_list$strand <- as.character(GRCh38_known_gene_list$strand)
GRCh38_known_gene_list$txStart <- as.character(GRCh38_known_gene_list$txStart)
GRCh38_known_gene_list$txEnd <- as.character(GRCh38_known_gene_list$txEnd)
GRCh38_known_gene_list$cdsStart <- as.character(GRCh38_known_gene_list$cdsStart)
GRCh38_known_gene_list$cdsEnd <- as.character(GRCh38_known_gene_list$cdsEnd)
GRCh38_known_gene_list$exonCount <- as.character(GRCh38_known_gene_list$exonCount)
GRCh38_known_gene_list$exonStarts <- as.character(GRCh38_known_gene_list$exonStarts)
GRCh38_known_gene_list$exonEnds <- as.character(GRCh38_known_gene_list$exonEnds)
GRCh38_known_gene_list$proteinID <- as.character(GRCh38_known_gene_list$proteinID)
GRCh38_known_gene_list$alignID <- as.character(GRCh38_known_gene_list$alignID)
GRCh38_known_gene_list$genesymbol <- as.character(GRCh38_known_gene_list$genesymbol)
# known MiHAs 
# known_MiHA_SNPs <- makeGRangesFromDataFrame(known_MiHA_table, keep.extra.columns = T,
#                                             ignore.strand = T, seqnames = "Chr",
#                                             start.field = "Pos", end.field = "Pos")
# seqlengths(known_MiHA_SNPs) <- chr.len[names(seqlengths(known_MiHA_SNPs))]
# pp <- autoplot(known_MiHA_SNPs, layout = "karyogram", geom = "rect") + scale_color_discrete("brown")## SNPs in the chromosome space

## ClinVar - AML related genes (31)
# clinVar.AML <- read.delim("../Data/clinvar_result_acuteMyeloidLeukemia.txt")
clinVar.AML <- read.delim("../Data//AcuteMyeloidLeukemia_relatedGenes_0314.txt", stringsAsFactors = F)
clinVar.AML <- read.delim("../ClinVar//AcuteMyeloidLeukemia_relatedGenes_0314.txt", stringsAsFactors = F)
AML_index <- which(grepl("acute myeloid leukemia", as.character(clinVar.AML$Condition.s.), ignore.case = T))
clinVar.AML <- clinVar.AML[AML_index, ] # 237 

# AML_genes <- GRCh38_known_gene_list[grepl(paste0(unique(as.character(clinVar.AML$Gene.s.)), collapse = "|"), GRCh38_known_gene_list$genesymbol),]
# AML_genes$Chromosome <- as.character(AML_genes$Chromosome)
# AML_genes <- AML_genes[!grepl("_", AML_genes$Chromosome), ]

### GRCh38_known_gene_list$genesymbol[grepl("HLA-A|HLA-B|HLA-C|HLA-DRB|HLA-DQB", GRCh38_known_gene_list$genesymbol)]

# GR_AML_genes <- makeGRangesFromDataFrame(AML_genes, seqnames.field = "Chromosome",
#                                            start.field = "txStart", end.field = "txEnd", 
#                                            strand.field = "strand",
#                                            keep.extra.columns = T)


no_range_locations <- which(as.character(clinVar.AML$Location) %in% "")
long_range_locations <- which(grepl("-", as.character(clinVar.AML$Location)))


no_range_genes <- GRCh38_known_gene_list[grepl(paste0(unique(as.character(clinVar.AML$Gene.s.[no_range_locations])), collapse = "|"), GRCh38_known_gene_list$genesymbol),]
no_range_genes$startIndex <- no_range_genes$txStart
no_range_genes$endIndex <- no_range_genes$txEnd

gene_list_point <- clinVar.AML[-c(long_range_locations, no_range_locations), ]
point_genes <- GRCh38_known_gene_list[grepl(paste0(unique(as.character(gene_list_point$Gene.s.)), 
                                                   collapse = "|"), GRCh38_known_gene_list$genesymbol), ]
point_genes_index <- unlist(sapply(1:dim(gene_list_point)[1], function(x) 
                      {iid <- which(as.character(point_genes$genesymbol) %in% as.character(gene_list_point$Gene.s.[x]))
                      if(length(iid) > 0) iid else 0})
                     )

# gene_list_point$Gene.s.[which(point_genes_index == 0) ]

num_point_genes <- dim(gene_list_point)[1]
point_all_genes <- data.frame(name = character(num_point_genes),
                              Chromosome = character(num_point_genes),
                              strand = character(num_point_genes),
                              txStart = numeric(num_point_genes),
                              txEnd = numeric(num_point_genes),
                              cdsStart = numeric(num_point_genes),
                              cdsEnd = numeric(num_point_genes),
                              exonCount = character(num_point_genes),
                              exonStarts = character(num_point_genes),
                              exonEnds = character(num_point_genes),
                              proteinID = character(num_point_genes),
                              alignID = character(num_point_genes),
                              genesymbol = character(num_point_genes),
                              # startIndex = as.numeric(num_point_genes),
                              # endIndex = as.numeric(num_point_genes),
                              stringsAsFactors = F)
for(id in 1:num_point_genes){
  
  if(point_genes_index[id] > 0){
    
    point_all_genes[id, ] <- point_genes[point_genes_index[id], ]
    # point_all_genes$startIndex[id] <- as.numeric(gene_list_point$Location[id])
    # point_all_genes$endIndex[id] <- as.numeric(gene_list_point$Location[id])
    
  }else {
    point_all_genes$genesymbol[id] <- gene_list_point$Gene.s.[id]
    point_all_genes$Chromosome[id] <- paste0("chr", gene_list_point$Chromosome[id])
    point_all_genes$strand[id] <- "+"
    # point_all_genes$startIndex[id] <- as.numeric(gene_list_point$Location[id])
    # point_all_genes$endIndex[id] <- as.numeric(gene_list_point$Location[id])
  }
  
}
point_all_genes$startIndex <- as.numeric(gene_list_point$Location)
point_all_genes$endIndex <- as.numeric(gene_list_point$Location)
# multi_genes <- GRCh38_known_gene_list[grepl(as.character(clinVar.AML$Gene.s.[long_range_locations[which(long_range_gene_index == 0)]]),
#                                             GRCh38_known_gene_list$genesymbol),]

long_range_genes <- GRCh38_known_gene_list[grepl(paste0(unique(clinVar.AML$Gene.s.[long_range_locations]), collapse = "|"), GRCh38_known_gene_list$genesymbol),]
# long_range_genes$name <- as.character(long_range_genes$name)
# long_range_genes$Chromosome <- as.character(long_range_genes$Chromosome)
# long_range_genes$strand <- as.character(long_range_genes$strand)
# long_range_genes$txStart <- as.character(long_range_genes$txStart)
# long_range_genes$txEnd <- as.character(long_range_genes$txEnd)
# long_range_genes$cdsStart <- as.character(long_range_genes$cdsStart)
# long_range_genes$cdsEnd <- as.character(long_range_genes$cdsEnd)
# long_range_genes$exonCount <- as.character(long_range_genes$exonCount)
# long_range_genes$exonStarts <- as.character(long_range_genes$exonStarts)
# long_range_genes$exonEnds <- as.character(long_range_genes$exonEnds)
# long_range_genes$proteinID <- as.character(long_range_genes$proteinID)
# long_range_genes$alignID <- as.character(long_range_genes$alignID)
# long_range_genes$genesymbol <- as.character(long_range_genes$genesymbol)

long_range_gene_index <- unlist(sapply(1:length(long_range_locations), function(x) 
{iid <- which(as.character(long_range_genes$genesymbol) %in% clinVar.AML$Gene.s.[long_range_locations[x]])

if(length(iid) > 0) iid else 0})
)

multi_genes <- GRCh38_known_gene_list[grepl(as.character(clinVar.AML$Gene.s.[long_range_locations[which(long_range_gene_index == 0)]]),
                             GRCh38_known_gene_list$genesymbol),]
long_range_genes <- rbind(long_range_genes, long_range_genes[dim(long_range_genes)[1],])
long_range_genes$name[dim(long_range_genes)[1]] <- paste0(multi_genes$name, collapse = "|")
long_range_genes$Chromosome[dim(long_range_genes)[1]] <- as.character(multi_genes$Chromosome[1])
long_range_genes$txStart[dim(long_range_genes)[1]] <- min(multi_genes$txStart)
long_range_genes$txEnd[dim(long_range_genes)[1]] <- max(multi_genes$txEnd)
long_range_genes$cdsStart[dim(long_range_genes)[1]] <- min(multi_genes$cdsStart)
long_range_genes$cdsEnd[dim(long_range_genes)[1]] <- max(multi_genes$cdsEnd)
long_range_genes$exonCount[dim(long_range_genes)[1]] <- paste0(multi_genes$exonCount, collapse = "|")
long_range_genes$exonStarts[dim(long_range_genes)[1]] <- paste0(multi_genes$exonStarts, collapse = "|")
long_range_genes$exonEnds[dim(long_range_genes)[1]] <- paste0(multi_genes$exonEnds, collapse = "|")
long_range_genes$proteinID[dim(long_range_genes)[1]] <- paste0(multi_genes$proteinID, collapse = "|")
long_range_genes$alignID[dim(long_range_genes)[1]] <- paste0(multi_genes$alignID, collapse = "|")
long_range_genes$genesymbol[dim(long_range_genes)[1]] <- paste0(multi_genes$genesymbol, collapse = "|")

long_range_gene_index[which(long_range_gene_index == 0)] <- dim(long_range_genes)[1]

long_range_all_genes <- long_range_genes[long_range_gene_index,]

long_range_indices <- do.call("rbind", strsplit(as.character(clinVar.AML$Location[long_range_locations]), " - "))
long_range_all_genes$startIndex <- as.numeric(long_range_indices[,1])
long_range_all_genes$endIndex <- as.numeric(long_range_indices[,2])

AML_genes <- rbind(no_range_genes, point_all_genes, long_range_all_genes)

GR_AML_genes <- makeGRangesFromDataFrame(AML_genes, seqnames.field = "Chromosome",
                                         start.field = "startIndex", end.field = "endIndex",
                                         strand.field = "strand",
                                         keep.extra.columns = T)

seqlengths(GR_AML_genes) <- chr.len[names(seqlengths(GR_AML_genes))]
save(GR_AML_genes, AML_genes, file = "../Data/AML_patientOnly_ClinVar_AML_genes_acuteML_0415.RData")

load("../Data/AML_patientOnly_ClinVar_AML_genes_acuteML_0415.RData")

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

##########################
# AML related genes LOD
#############
num_AML_genes <- dim(AML_genes)[1]

# for(id in 1:num_AML_genes){
#   
#   chrom <- AML_genes$Chromosome[id]
#   
#   recipient_chrom_id <- which(recipient_missense_stats$CHROM == chrom)
#   donor_chrom_id <- which(donor_missense_stats$CHROM %in% chrom)
#   
#   aa <- which(recipient_missense_stats$POS[recipient_chrom_id] >= as.numeric(AML_genes$txStart[id]) & 
#           recipient_missense_stats$POS[recipient_chrom_id] <= as.numeric(AML_genes$txEnd[id]))
#   bb <- which(donor_missense_stats$POS[donor_chrom_id] >= as.numeric(AML_genes$txStart[id]) & 
#                 donor_missense_stats$POS[donor_chrom_id] <= as.numeric(AML_genes$txEnd[id]))
#   
# }
num_AML_points <- dim(gene_list_point)[1]
# strsplit(strsplit(gene_list_point$Name[1],":")[[1]][2], " ")

AML_point_mutation_table <- data.frame(CHROM = paste0("chr", gene_list_point$Chromosome),
                                       POS = as.numeric(gene_list_point$Location),
                                       Gene_name = gene_list_point$Name,
                                       AMLpatient_count = numeric(num_AML_points),
                                       donor_count = numeric(num_AML_points),
                                       AMLpatient_pathogenic_count = numeric(num_AML_points),
                                       donor_pathogenic_count = numeric(num_AML_points),
                                       Condition = gene_list_point$Condition.s.,
                                       Clinical.Significance.Last.Reviewed = gene_list_point$Clinical.significance..Last.reviewed.,
                                       stringsAsFactors = F)
for(id in 1:num_AML_points){
  
  CHROM <- AML_point_mutation_table$CHROM[id]
    
  Mutation <- unlist(strsplit(unlist(strsplit(strsplit(strsplit(AML_point_mutation_table$Gene_name[id],":")[[1]][2], " ")[[1]][1], split = "[0-9]+"))[2], ">"))
  
  recipient_chrom_id <- which(recipient_missense_stats$CHROM == CHROM)
  donor_chrom_id <- which(donor_missense_stats$CHROM %in% CHROM)
  
  if(length(recipient_chrom_id)>0){
    
    r_index <- which(recipient_missense_stats$POS[recipient_chrom_id] == AML_point_mutation_table$POS[id])
    
    if(length(r_index) > 0){
      AML_point_mutation_table$AMLpatient_count[id] <- recipient_missense_stats$NumDiff[recipient_chrom_id[r_index]]
      
      switch(Mutation[2],
             "A" = AML_point_mutation_table$AMLpatient_pathogenic_count[id] <- recipient_missense_stats$ALT.A[recipient_chrom_id[r_index]], 
             "T" = AML_point_mutation_table$AMLpatient_pathogenic_count[id] <- recipient_missense_stats$ALT.T[recipient_chrom_id[r_index]], 
             "G" = AML_point_mutation_table$AMLpatient_pathogenic_count[id] <- recipient_missense_stats$ALT.G[recipient_chrom_id[r_index]], 
             "C" = AML_point_mutation_table$AMLpatient_pathogenic_count[id] <- recipient_missense_stats$ALT.C[recipient_chrom_id[r_index]]) 
    }
    
  }
  
  if(length(donor_chrom_id)>0){
    
    d_index <- which(donor_missense_stats$POS[donor_chrom_id] == gene_list_point$Location[id])
    
    if(length(d_index) > 0){
      AML_point_mutation_table$donor_count[id] <- donor_missense_stats$NumDiff[donor_chrom_id[d_index]]
      
      switch(Mutation[2],
             "A" = AML_point_mutation_table$donor_pathogenic_count[id] <- donor_missense_stats$ALT.A[donor_chrom_id[d_index]], 
             "T" = AML_point_mutation_table$donor_pathogenic_count[id] <- donor_missense_stats$ALT.T[donor_chrom_id[d_index]], 
             "G" = AML_point_mutation_table$donor_pathogenic_count[id] <- donor_missense_stats$ALT.G[donor_chrom_id[d_index]], 
             "C" = AML_point_mutation_table$donor_pathogenic_count[id] <- donor_missense_stats$ALT.C[donor_chrom_id[d_index]]) 
      
    }
    
  }
  
}

write.csv(AML_point_mutation_table, file = "../ClinVar/AML_patient_donor_point_Mutation_counts.csv", row.names = F)


Unique_AML_genes <- unique(AML_genes)
num_AML_genes <- dim(Unique_AML_genes)[1]
AML_genes_mutation_table <- data.frame(CHROM = Unique_AML_genes$Chromosome,
                                       strand = Unique_AML_genes$strand,
                                       txStart = as.numeric(Unique_AML_genes$txStart),
                                       txEnd = as.numeric(Unique_AML_genes$txEnd),
                                       startIndex = as.numeric(Unique_AML_genes$startIndex),
                                       endIndex = as.numeric(Unique_AML_genes$endIndex),
                                       exonCount = Unique_AML_genes$exonCount,
                                       Gene_name = Unique_AML_genes$genesymbol,
                                       AMLpatient_count = numeric(num_AML_genes), 
                                       AMLPatient_mutation_POS_count = numeric(num_AML_genes),
                                       donor_count = numeric(num_AML_genes),
                                       donor_mutation_POS_count = numeric(num_AML_genes),
                                       # AMLpatient_pathogenic_count = numeric(AML_genes),
                                       # donor_pathogenic_count = numeric(AML_genes),
                                       # Condition = gene_list_point$Condition.s.,
                                       # Clinical.Significance.Last.Reviewed = gene_list_point$Clinical.significance..Last.reviewed.,
                                       stringsAsFactors = F)
for(id in 1:num_AML_genes){

  chrom <- AML_genes_mutation_table$CHROM[id]

  recipient_chrom_id <- which(recipient_missense_stats$CHROM == chrom)
  donor_chrom_id <- which(donor_missense_stats$CHROM %in% chrom)

  if(length(recipient_chrom_id) > 0){
    
    if(AML_genes_mutation_table$txStart[id] >0){ # one single gene
      
      r_ind <- which(recipient_missense_stats$POS[recipient_chrom_id] >= AML_genes_mutation_table$txStart[id] &
                       recipient_missense_stats$POS[recipient_chrom_id] <= AML_genes_mutation_table$txEnd[id])
      if(length(r_ind) >0){
        
        AML_genes_mutation_table$AMLPatient_mutation_POS_count[id] <- length(r_ind)
        AML_genes_mutation_table$AMLpatient_count[id] <- sum(recipient_missense_stats$NumDiff[recipient_chrom_id[r_ind]])
        
      }
      
    }else{ # multiple genes
      
      r_ind <- which(recipient_missense_stats$POS[recipient_chrom_id] >= AML_genes_mutation_table$startIndex[id] &
                       recipient_missense_stats$POS[recipient_chrom_id] <= AML_genes_mutation_table$endIndex[id])
      if(length(r_ind) > 0){
        
        AML_genes_mutation_table$AMLPatient_mutation_POS_count[id] <- length(r_ind)
        AML_genes_mutation_table$AMLpatient_count[id] <- sum(recipient_missense_stats$NumDiff[recipient_chrom_id[r_ind]])
        
      }
      
    }
    
  }
  
  if(length(donor_chrom_id) > 0){
    
    if(AML_genes_mutation_table$txStart[id] >0){ # one single gene
      
      d_ind <- which(donor_missense_stats$POS[donor_chrom_id] >= AML_genes_mutation_table$txStart[id] &
                       donor_missense_stats$POS[donor_chrom_id] <= AML_genes_mutation_table$txEnd[id])
      if(length(d_ind) >0){
        
        AML_genes_mutation_table$donor_mutation_POS_count[id] <- length(d_ind)
        AML_genes_mutation_table$donor_count[id] <- sum(donor_missense_stats$NumDiff[donor_chrom_id[d_ind]])
        
      }
      
    }else{ # multiple genes
      
      d_ind <- which(donor_missense_stats$POS[donor_chrom_id] >= AML_genes_mutation_table$startIndex[id] &
                       donor_missense_stats$POS[donor_chrom_id] <= AML_genes_mutation_table$endIndex[id])
      if(length(d_ind) > 0){
        
        AML_genes_mutation_table$donor_mutation_POS_count[id] <- length(d_ind)
        AML_genes_mutation_table$donor_count[id] <- sum(donor_missense_stats$NumDiff[donor_chrom_id[d_ind]])
        
      }
      
    }
    
  }

}

write.csv(AML_genes_mutation_table, file = "../ClinVar/AML_gene_all_region_mutation_counts.csv", row.names = F)

###### AML patient counts at each gene - cloud
source('util.R', echo = FALSE)

# load("../Data/GRCh38_gene_list.RData")

# library(vcfR)

VCF_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/wwang/hli_annotated/"
output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants_wwVersion/"

all_vcf_files <- list.files(VCF_file_dir, pattern = "\\.vcf$")

file_names <- all_vcf_files[grepl("_R_", all_vcf_files)] # 216 donors; 240 recipients

# AML patients only
load("../Data/HLI_available_pairs_dis_table.RData")

AML_ids <- which(Available_paired_table$Disease %in% "AML") # 196
RID <- sapply(1:length(file_names), function(x) as.numeric(unlist(strsplit(file_names[x], "_"))[4]))

AML_RID_ids <- which(RID %in% Available_paired_table$RID[AML_ids])
file_names <- file_names[AML_RID_ids]
#####
num_files <- length(file_names)  
AML_gene_stats <- data.frame(GroupID = numeric(),
                             CHROM=character(0), POS=numeric(0), REF = character(0),
                             ALT.A=numeric(0), ALT.T=numeric(0), ALT.G = numeric(0), ALT.C = numeric(0),
                             ALT.insertion = numeric(0), ALT.deletion = numeric(0), 
                             Single_double = numeric(0),
                             NumDiff = numeric(0), stringsAsFactors = F)


for(id in 1:num_files){
  
  vcf_file <- paste0(VCF_file_dir, file_names[id])
  # vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  # temp_info <- vcf_info@fix
  # vcf_missense <- temp_info[temp_info[,"FILTER"] == "PASS", c("CHROM", "POS", "REF","ALT")]
  # 
  vcf_missense <- read.table(vcf_file)
  colnames(vcf_missense) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "ANN", "INFO", "GT")
  
  vcf_missense$CHROM <- as.character(vcf_missense$CHROM)
  vcf_missense$POS <- as.numeric(as.character(vcf_missense$POS))
  vcf_missense$REF <- as.character(vcf_missense$REF)
  vcf_missense$ALT <- as.character(vcf_missense$ALT)
  vcf_missense$QUAL <- as.numeric(as.character(vcf_missense$QUAL))
  vcf_missense$ANN <- as.character(vcf_missense$ANN)
  vcf_missense$INFO <- as.character(vcf_missense$INFO)
  vcf_missense$GT <- as.character(vcf_missense$GT)
  
  for(jd in 1:num_AML_genes){
    
    chrom <- AML_genes_mutation_table$CHROM[id]
    
    chrom_id <- which(vcf_missense$CHROM == chrom)
    # donor_chrom_id <- which(donor_missense_stats$CHROM %in% chrom)
    
    if(length(chrom_id) > 0){
      
      if(AML_genes_mutation_table$txStart[id] >0){ # one single gene
        
        ind <- which(vcf_missense$POS[chrom_id] >= AML_genes_mutation_table$txStart[id] &
                         vcf_missense$POS[chrom_id] <= AML_genes_mutation_table$txEnd[id])
        if(length(ind) >0){
          
          AML_genes_mutation_table$AMLPatient_mutation_POS_count[id] <- length(ind)
          AML_genes_mutation_table$AMLpatient_count[id] <- sum(recipient_missense_stats$NumDiff[recipient_chrom_id[r_ind]])
          
        }
        
      }else{ # multiple genes
        
        r_ind <- which(recipient_missense_stats$POS[recipient_chrom_id] >= AML_genes_mutation_table$startIndex[id] &
                         recipient_missense_stats$POS[recipient_chrom_id] <= AML_genes_mutation_table$endIndex[id])
        if(length(r_ind) > 0){
          
          AML_genes_mutation_table$AMLPatient_mutation_POS_count[id] <- length(r_ind)
          AML_genes_mutation_table$AMLpatient_count[id] <- sum(recipient_missense_stats$NumDiff[recipient_chrom_id[r_ind]])
          
        }
        
      }
      
    }
    
  }
  
  
  
  num_variants <- dim(vcf_missense)[1]
  
  # print(head(vcf_missense))
  if(dim(missense_stats)[1] > 0){ # if not the first sample
    vcf_chrom_pos <- sapply(1:dim(vcf_missense)[1], function(x) paste0(vcf_missense[x, c("CHROM","POS")], collapse = ""))
    all_chrom_pos <- sapply(1:dim(missense_stats)[1], function(x) paste0(missense_stats[x, c("CHROM","POS")], collapse = ""))
    
    IDs_in_vcf <- which(vcf_chrom_pos %in% all_chrom_pos)
    IDs_in_all <- which(all_chrom_pos %in% vcf_chrom_pos) 
    
    missense_mutations <- count_missense_mutation_type(vcf_missense[IDs_in_vcf,])
    
    ### existing records increase the numbers
    ordered_IDs_in_all <- sapply(1:length(IDs_in_vcf), function(x) which(all_chrom_pos[IDs_in_all] == vcf_chrom_pos[IDs_in_vcf[x]]))
    
    missense_stats$NumDiff[IDs_in_all[ordered_IDs_in_all]] <- missense_stats$NumDiff[IDs_in_all[ordered_IDs_in_all]] + 1
    
    missense_stats[IDs_in_all[ordered_IDs_in_all], c("ALT.A", "ALT.T", "ALT.G", "ALT.C", "ALT.insertion", "ALT.deletion", "Single_double")] <- 
      missense_stats[IDs_in_all[ordered_IDs_in_all], c("ALT.A", "ALT.T", "ALT.G", "ALT.C", "ALT.insertion", "ALT.deletion", "Single_double")] +
      missense_mutations[, c("ALT.A", "ALT.T", "ALT.G", "ALT.C", "ALT.insertion", "ALT.deletion", "Single_double")]
    
    ### new records added to the bottom
    if(length(IDs_in_vcf) < length(vcf_chrom_pos)){ # new POS
      
      num_new_pos <- length(vcf_chrom_pos) -  length(IDs_in_vcf)
      missense_mutations <- count_missense_mutation_type(vcf_missense[-IDs_in_vcf,])
      
      temp_stats <- data.frame(vcf_missense[-IDs_in_vcf, c("CHROM", "POS")], missense_mutations, NumDiff = numeric(num_new_pos)+1, 
                               stringsAsFactors = F)
      
      temp_stats$CHROM <- as.character(temp_stats$CHROM)
      temp_stats$POS <- as.numeric(as.character(temp_stats$POS))
      
      missense_stats <- rbind(missense_stats, temp_stats)
      
    }
    
  } else { # the first sample
    
    
    missense_mutations <- count_missense_mutation_type(vcf_missense)
    # het_index <- as.numeric(which(sapply(1:num_variants, function(x) nchar(vcf_missense[x, "REF"]) != nchar(vcf_missense[x, "ALT"]))))
    # het_num <- length(het_index)
    # hom_index <- (1:num_variants)[-het_index]
    # hom_num <- length(hom_index)
    # missense_mutatations <- data.frame(REF = character(num_variants), ALT.A = numeric(num_variants), 
    #                                   ALT.T=numeric(num_variants), ALT.G=numeric(num_variants), 
    #                                   ALT.C=numeric(num_variants), stringsAsFactors = F)
    # # A - 1; T - 2; G - 3; C - 4;
    # 
    # # homozygous variants
    # missense_mutatations$REF[hom_index] <- vcf_missense[hom_index, "REF"]
    # for(hom_id in hom_index){
    #   missense_mutatations[hom_id, which(c("A", "T", "G", "C") %in% vcf_missense[hom_id, "ALT"])+1] <- 
    #     missense_mutatations[hom_id, which(c("A", "T", "G", "C") %in% vcf_missense[hom_id, "ALT"])+1] + 1
    # }
    # 
    # # heterozygous variants
    # missense_mutatations$REF[het_index] <- vcf_missense[het_index, "REF"]
    # for(het_id in het_index){
    #   missense_mutatations[het_id, which(c("A", "T", "G", "C") %in% unlist(strsplit(vcf_missense[het_id, "ALT"], ",")))+1] <- 
    #     missense_mutatations[het_id, which(c("A", "T", "G", "C") %in% unlist(strsplit(vcf_missense[het_id, "ALT"], ",")))+1] + 1
    # }
    
    
    missense_stats <- rbind(missense_stats, cbind(vcf_missense[, c("CHROM", "POS")], missense_mutations, NumDiff = numeric(num_variants)+1))
    
    missense_stats$CHROM <- as.character(missense_stats$CHROM)
    missense_stats$POS <- as.numeric(as.character(missense_stats$POS))
    # missense_stats[, c("ALT.A", "ALT.T", "ALT.G", "ALT.C", "NumDiff")] <- as.numeric(as.character(missense_stats[, c("ALT.A", "ALT.T", "ALT.G", "ALT.C", "NumDiff")]))
    
  }
  
}

recipient_missense_stats <- missense_stats
save(recipient_missense_stats, file = paste0(output_dir, "AMLPatients_missesense_stats_wwVersion.RData"))