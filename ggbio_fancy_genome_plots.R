library("ggbio")
# load(system.file("data", "hg19IdeogramCyto.rda", package="biovizBase", mustWork=TRUE))
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(biovizBase)

# source("https://bioconductor.org/biocLite.R")
# biocLite("Homo.sapiens")
# library(Homo.sapiens)
# class(Homo.sapiens)


## Plotting a single series in bar format
# if (!require("IdeoViz", quietly=TRUE, warn.conflicts = FALSE)) {
#   source("https://bioconductor.org/biocLite.R")
#   biocLite("IdeoViz")
#   library("IdeoViz", verbose=F, warn.conflicts =F)
# }

# data(binned_singleSeries) # GRanges class
# data(hg18_ideo)
# data(UCSC.HG38.Human.CytoBandIdeogram) # from RCircos
# data(hg38IdeogramCyto, package = "biovizBase")
# hg38_cytoband <- UCSC.HG38.Human.CytoBandIdeogram # cytoBandIdeo table downloaded previously and stored as a data.frame.
# names(hg38_cytoband) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")
# hg38_cytoband <- makeGRangesFromDataFrame(hg38_cytoband, keep.extra.columns = T)  # Ideogram format of HG38
# biovizBase::isIdeogram(hg38_cytoband)
# save(hg38_cytoband, UCSC.HG38.Human.CytoBandIdeogram, file = "../Data/hg38_cytobandIdeogram.RData")
#### hg38 cytoband data saved.
# Ideogram(hg38_cytoband, xlabel = T, subchr = "chr1")

hg38_cytoband <- getIdeogram("hg38")
Ideogram(hg38_cytoband, xlabel = T, subchr = "chr1")
# data(genesymbol, package = "biovizBase")
#### Gene list
genome <- BSgenome.Hsapiens.UCSC.hg38

hg38_known_genes <- TxDb.Hsapiens.UCSC.hg38.knownGene
knownGene_tx <- extractUpstreamSeqs(genome, hg38_known_genes)

refGene_txdb <- makeTxDbFromUCSC("hg38", "refGene")
refGene_up1000seqs <- extractUpstreamSeqs(genome, refGene_txdb)
# load("../Data/GRCh38_gene_list.RData")
GRCh38_known_gene_list <- read.delim("../UCSCgenome/database/knownGene.txt", header  = F)
kg_table <- read.delim("../UCSCgenome/database/hgTables.txt")
# GeneSymbols <- sapply(GRCh38_known_gene_list$name, function(x) kg_table$geneSymbol[which(as.character(x) == as.character(kg_table$X.kgID))])

num_genes <- dim(GRCh38_known_gene_list)[1]
GeneSymbols <- data.frame(genesymbol = character(num_genes),
                          stringsAsFactors = F)
for(id in 1:num_genes){
  index <- which(as.character(GRCh38_known_gene_list$name[id]) == as.character(kg_table$X.kgID))
  if(length(index) == 1){
  GeneSymbols$genesymbol[id] <- as.character(kg_table$geneSymbol[index])
  } else GeneSymbols$genesymbol[id] <- NA
}
GRCh38_known_gene_list$genesymbol <- GeneSymbols$genesymbol
colnames(GRCh38_known_gene_list) <- c("name","Chromosome", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                                      "exonCount", "exonStarts", "exonEnds", "proteinID", "alignID", "genesymbol")
# rownames(GRCh38_known_gene_list) <- GeneSymbols$genesymbol
GRCh38_known_genesymbol <- makeGRangesFromDataFrame(GRCh38_known_gene_list, seqnames.field = "Chromosome",
                                                    start.field = "txStart", end.field = "txEnd", 
                                                    strand.field = "strand",
                                                    keep.extra.columns = T)
save(GRCh38_known_genesymbol, GRCh38_known_gene_list, file = "../Data/hg38_known_gene_symbols.RData")
#### HG38 Known gene list saved based on hgTable from USCS



### hg38 knwon gene names from HUGO Gene Nomenclature Committee
GRCh38_known_gene_list <- read.delim("../UCSCgenome/database/knownGene.txt", header  = F)
HGNC_table <- read.delim("../UCSCgenome/HUGO_Gene_Nomenclature_Committee.txt")
# GeneSymbols <- sapply(GRCh38_known_gene_list$name, function(x) kg_table$geneSymbol[which(as.character(x) == as.character(kg_table$X.kgID))])

num_genes <- dim(GRCh38_known_gene_list)[1]
GeneSymbols <- data.frame(genesymbol = character(num_genes),
                          stringsAsFactors = F)
for(id in 1:num_genes){
  index <- which(as.character(GRCh38_known_gene_list$name[id]) == as.character(HGNC_table$UCSC.ID.supplied.by.UCSC.))
  if(length(index) == 1){
    GeneSymbols$genesymbol[id] <- as.character(HGNC_table$Approved.Symbol[index])
  } else GeneSymbols$genesymbol[id] <- NA
}

GRCh38_known_gene_list$genesymbol <- GeneSymbols$genesymbol
colnames(GRCh38_known_gene_list) <- c("name","Chromosome", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd",
                                      "exonCount", "exonStarts", "exonEnds", "proteinID", "alignID", "genesymbol")
# rownames(GRCh38_known_gene_list) <- GeneSymbols$genesymbol
GRCh38_known_genesymbol <- makeGRangesFromDataFrame(GRCh38_known_gene_list, seqnames.field = "Chromosome",
                                                    start.field = "txStart", end.field = "txEnd", 
                                                    strand.field = "strand",
                                                    keep.extra.columns = T)
save(GRCh38_known_genesymbol, GRCh38_known_gene_list, file = "../Data/hg38_known_gene_symbols_HGNC.RData")

#####
MHC_regions <- GRCh38_known_gene_list[grepl("HLA-A|HLA-B|HLA-C|HLA-DRB|HLA-DQB", GRCh38_known_gene_list$genesymbol),]
MHC_regions <- MHC_regions[grep("chr6$", MHC_regions$Chromosome),]
MHC_regions$Chromosome <- as.character(MHC_regions$Chromosome)
# GRCh38_known_gene_list$genesymbol[grepl("HLA-A|HLA-B|HLA-C|HLA-DRB|HLA-DQB", GRCh38_known_gene_list$genesymbol)]
GR_MHC_regions <- makeGRangesFromDataFrame(MHC_regions, seqnames.field = "Chromosome",
                                           start.field = "txStart", end.field = "txEnd", 
                                           strand.field = "strand",
                                           keep.extra.columns = T)
seqlengths(GR_MHC_regions) <- chr.len[names(seqlengths(GR_MHC_regions))]
save(GR_MHC_regions, MHC_regions, file = "../Data/hg38_MHC_region_gene_symbols_HGNC.RData")
#######



#### all missense 
# library(ggplot2)
known_MiHA_table <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T)

load("../WW_MiHA/summary/Data/missense_summary_DtoR.RData")

aGVHD_mc <- do.call("rbind", aGVHD_SNP_list)
aGVHD_mc$POS2 <- aGVHD_mc$POS+1
aGVHD_mc$POS3 <- aGVHD_mc$POS
aGVHD_mc$group <- "aGVHD"
aGVHD_mc$groupColor <- "#D55E00"
nGVHD_mc <- do.call("rbind", nGVHD_SNP_list)
nGVHD_mc$POS2 <- nGVHD_mc$POS+1
nGVHD_mc$POS3 <- nGVHD_mc$POS
nGVHD_mc$NumDiff <- -nGVHD_mc$NumDiff 
nGVHD_mc$group <- "non-aGVHD"
nGVHD_mc$groupColor <- "#0072B2"
all_mc <- rbind(aGVHD_mc, nGVHD_mc)
# groupColor <- c(rep("#D55E00", length(aGVHD_mc$POS2)), rep("#0072B2", length(nGVHD_mc$POS2)))
# all_mc <- all_mc[-which(all_mc$CHROM == "chrX"), ] # sex chromosomes
# all_mc <- all_mc[-which(all_mc$CHROM == "chrY"), ]

GR_all_mc <- makeGRangesFromDataFrame(all_mc, keep.extra.columns = T, 
                                      ignore.strand = T, seqnames.field = "CHROM",
                                      start.field = "POS", end.field = "POS"
                                      )

chr.len = seqlengths(BSgenome.Hsapiens.UCSC.hg38)  # get chromosome lengths
# remove X,Y,M and random chromosomes
chr.len = chr.len[grep("_|M", names(chr.len), invert = T)]
seqlengths(GR_all_mc) <- chr.len
###### Ideogram + missense mismatches
# library(ggbio)
# library(scales)
load("../Data/hg38_cytobandIdeogram.RData")
# load("../Data/hg38_known_gene_symbols_HGNC.RData")
load("../Data/hg38_MHC_region_gene_symbols_HGNC.RData")
# hg38<- keepSeqlevels(hg38_cytoband, paste0("chr", c(1:22)))  ## extract only chr1~22
# head(hg38)
# p <- ggplot(hg38) + layout_karyogram(cytoband = F) + ## base karyogram for 22 chromosomes
#   layout_karyogram(GR_all_mc, geom = "rect", color = "red")
# 
# autoplot(hg38, layout = "karyogram", cytoband = TRUE)

# seqlevels(GR_all_mc, force=TRUE) = paste0("chr", 1:22) ## extract only autosomes

known_MiHA_SNPs <- makeGRangesFromDataFrame(known_MiHA_table, keep.extra.columns = T,
                                            ignore.strand = T, seqnames = "Chr",
                                            start.field = "Pos", end.field = "Pos")
seqlengths(known_MiHA_SNPs) <- chr.len[names(seqlengths(known_MiHA_SNPs))]
# pp <- autoplot(known_MiHA_SNPs, layout = "karyogram", geom = "rect") + scale_color_discrete("brown")## SNPs in the chromosome space
hg38<- keepSeqlevels(hg38_cytoband, paste0("chr", c(1:22, "X", "Y")))  ## extract only chr1~22
pp <- autoplot(hg38, layout = "karyogram", cytoband = T)#+ scale_fill_giemsa()

# pp <- autoplot(seqinfo(GR_all_mc))
# colorValues <- c(alpha(c("#D55E00", "#0072B2"), .1), 
#                  getOption("biovizBase")$cytobandColor[c("acen", "gneg", "gpos100", "gpos25", "gpos50", "gpos75", "gvar","stalk")] 
#                  )
# names(colorValues) <- NULL
pp + layout_karyogram(known_MiHA_SNPs, layout = "karyogram", size = 1, color = "red") + #scale_color_discrete("red") +
  layout_karyogram(GR_MHC_regions, layout = "karyogram", size = 1, color = "#009E73") +
  layout_karyogram(GR_all_mc, aes(x = POS3, y = NumDiff, color = factor(group)),
                   ylim = c(10,40), rect.height = 5, geom="point", size = 0.5, alpha = 0.1) + 
  # scale_colour_manual(values = alpha(c("#D55E00", "#0072B2"), .1)) +
  # scale_fill_manual(values = alpha(c("#D55E00", "#0072B2"), .1)) +
  theme(#axis.title.x=element_blank(),
    legend.position = "none")

  # scale_fill_continuous(guide = guide_legend(title = "Groups"))

##################################
# log ratio
##################################
# library(ggplot2)
source("util.R")
known_MiHA_table <- read.table("../WW_MiHA/Known_MiHA_coordinates.txt", header = T)

load("../WW_MiHA/summary/missense_summary_DtoR.RData")

all_MiHA_SNP <- data.frame()
log_ratio_all <- data.frame()
for(id in 1:22){
  chr <- paste0("chr", id)
  
  aGVHD_Mutation <- aGVHD_SNP_list[[chr]]
  nGVHD_Mutation <- nGVHD_SNP_list[[chr]]
  
  known_MiHA_chr <- known_MiHA_table[which(known_MiHA_table$Chr %in% chr), ]
  known_MiHA_chr$NumDiff <- rep(0, length(known_MiHA_chr$Pos))
  names(known_MiHA_chr)[5] <- "POS"
  # known_MiHA_chr$POS <- factor(known_MiHA_chr$POS)
  # MiHA_SNP <- known_MiHA_chr$POS
  # known_MiHA_chr$group <- "aGVHD"
  
  # mc <- data.frame(position = c(aGVHD_Mutation$POS, nGVHD_Mutation$POS),
  #                  frequency = c(aGVHD_Mutation$NumDiff, -nGVHD_Mutation$NumDiff),
  #                  group = c(rep("aGVHD", length(aGVHD_Mutation$POS)), rep("nGVHD", length(nGVHD_Mutation$POS))), 
  #                  stringsAsFactors = F)
  inter_pos <- sort(intersect(aGVHD_Mutation$POS, nGVHD_Mutation$POS), decreasing = F)
  
  log_ratio <- data.frame(CHROM = rep(chr, length(aGVHD_Mutation$POS+length(nGVHD_Mutation$POS)), 1), 
                          POS = c(aGVHD_Mutation$POS, nGVHD_Mutation$POS), 
                          log_ratio = numeric(length(aGVHD_Mutation$POS) + length(nGVHD_Mutation$POS)),
                          group = character(length(aGVHD_Mutation$POS) + length(nGVHD_Mutation$POS)),
                          stringsAsFactors = F)
  
  aGVHD_only_id <- which(!(aGVHD_Mutation$POS %in% inter_pos))
  nGVHD_only_id <- which(!(nGVHD_Mutation$POS %in% inter_pos))
  
  log_ratio$log_ratio[aGVHD_only_id] <- aGVHD_Mutation$NumDiff[aGVHD_only_id]
  log_ratio$group[aGVHD_only_id] <- "aGVHD_only"
  log_ratio$log_ratio[nGVHD_only_id + length(aGVHD_Mutation$POS)] <- -nGVHD_Mutation$NumDiff[nGVHD_only_id]
  log_ratio$group[nGVHD_only_id + length(aGVHD_Mutation$POS)] <- "nGVHD_only"
  # intersect(aGVHD_Mutation$POS[aGVHD_only_id],  nGVHD_Mutation$POS[nGVHD_only_id])
  
  for(inter_id in 1:length(inter_pos)){
    
    aGVHD_inter_id <- which(aGVHD_Mutation$POS %in% inter_pos[inter_id])
    nGVHD_inter_id <- which(nGVHD_Mutation$POS %in% inter_pos[inter_id])
    
    log_ratio$log_ratio[aGVHD_inter_id] <- log10(aGVHD_Mutation$NumDiff[aGVHD_inter_id] / nGVHD_Mutation$NumDiff[nGVHD_inter_id])
    log_ratio$group[aGVHD_inter_id] <- "log_odds"
  }
  log_ratio_all <- rbind(log_ratio_all, log_ratio)
}

dim(log_ratio_all)
log_ratio_all$POS2 <- log_ratio_all$POS
log_ratio_all <- log_ratio_all[-which(log_ratio_all$group == ""), ]
# groupOnly_count <- log_ratio_all[-which(log_ratio_all$group=="log_odds"),]
# groupOnly_count$group <- factor(groupOnly_count$group , levels=factor(groupOnly_count$group))
GR_log_ratio_groupOnly <- makeGRangesFromDataFrame(log_ratio_all[-which(log_ratio_all$group=="log_odds"),], keep.extra.columns = T, 
                                      ignore.strand = T, seqnames.field = "CHROM",
                                      start.field = "POS", end.field = "POS")
GR_log_ratio_group <- makeGRangesFromDataFrame(log_ratio_all[which(log_ratio_all$group=="log_odds"),], keep.extra.columns = T, 
                                                   ignore.strand = T, seqnames.field = "CHROM",
                                                   start.field = "POS", end.field = "POS")

chr.len = seqlengths(BSgenome.Hsapiens.UCSC.hg38)  # get chromosome lengths
# remove X,Y,M and random chromosomes
chr.len = chr.len[grep("_|M|X|Y", names(chr.len), invert = T)]
seqlengths(GR_log_ratio_group) <- chr.len
seqlengths(GR_log_ratio_groupOnly) <- chr.len

known_MiHA_SNPs <- makeGRangesFromDataFrame(known_MiHA_table, keep.extra.columns = T,
                                            ignore.strand = T, seqnames = "Chr",
                                            start.field = "Pos", end.field = "Pos")
seqlengths(known_MiHA_SNPs) <- chr.len[names(seqlengths(known_MiHA_SNPs))]
# pp <- autoplot(known_MiHA_SNPs, layout = "karyogram", geom = "rect") + scale_color_discrete("brown")## SNPs in the chromosome space
# pp <- autoplot(seqinfo(GR_log_ratio_groupOnly))
hg38_auto<- keepSeqlevels(hg38_cytoband, paste0("chr", c(1:22)))  ## extract only chr1~22
pp <- autoplot(hg38_auto, layout = "karyogram", cytoband = T)

pp + layout_karyogram(known_MiHA_SNPs, layout = "karyogram", size = 1, color = "brown") +
  layout_karyogram(GR_MHC_regions, layout = "karyogram", size = 1, color = "#009E73") +
  layout_karyogram(GR_log_ratio_groupOnly, aes(x = POS2, y = log_ratio, color=factor(group)),
                   ylim = c(10,40), rect.height = 5, geom="point", size = 1, alpha = 0.1) + 
  # scale_colour_manual(value = alpha(c("#D55E00", "#0072B2"), .2)) +
  # scale_fill_manual(value = alpha(c("#D55E00", "#0072B2"), .2)) +
  theme(#axis.title.x=element_blank(),
        legend.position = "none")


### Log-ratio
# pp <- autoplot(seqinfo(GR_log_ratio_group))
pp + layout_karyogram(known_MiHA_SNPs, layout = "karyogram", size = 1, color = "brown") +
  layout_karyogram(GR_MHC_regions, layout = "karyogram", size = 1, color = "#009E73") +
  layout_karyogram(GR_log_ratio_group, aes(x = POS2, y = log_ratio, color=group),
                   ylim = c(10,40), rect.height = 5, geom="point", color = "#D55E00", size = 1, alpha = 0.1) + 
  # scale_colour_manual(values = alpha(c("#D55E00"), .2)) +
  # scale_fill_manual(values = alpha(c("#D55E00"), .2)) +
  theme(axis.title.x=element_blank(),
        legend.position = "none")

#### circular plot - overall MiHAs

all_MiHA_SNP <- data.frame()
for(id in 1:22){
  chr <- paste0("chr", id)
  
  aGVHD_Mutation <- aGVHD_SNP_list[[chr]]
  nGVHD_Mutation <- nGVHD_SNP_list[[chr]]
  
  known_MiHA_chr <- known_MiHA_table[which(known_MiHA_table$Chr %in% chr), ]
  known_MiHA_chr$NumDiff <- rep(0, length(known_MiHA_chr$Pos))
  names(known_MiHA_chr)[5] <- "POS"
  # known_MiHA_chr$POS <- factor(known_MiHA_chr$POS)
  # MiHA_SNP <- known_MiHA_chr$POS
  # known_MiHA_chr$group <- "aGVHD"
  
  mc <- data.frame(position = c(aGVHD_Mutation$POS, nGVHD_Mutation$POS),
                   frequency = c(aGVHD_Mutation$NumDiff, -nGVHD_Mutation$NumDiff),
                   group = c(rep("aGVHD", length(aGVHD_Mutation$POS)), rep("nGVHD", length(nGVHD_Mutation$POS))), 
                   stringsAsFactors = F)
  # mc$frequency[mc$position == aGVHD_Mutation$POS] <- aGVHD_Mutation$NumDiff
  # mc$frequency[mc$position == nGVHD_Mutation$POS] <- -nGVHD_Mutation$NumDiff
  # ggplot(aGVHD_Mutation, aes(x = POS, ymax = NumDiff, ymin = 0)) + 
  #   geom_linerange()
  if(min(mc$position) > 1){
    mc <- rbind(mc, data.frame(position = 1,
                               frequency = 0,
                               group = "nGVHD"))
  }
  
  if(max(mc$position) < chr.len[chr]){
    mc <- rbind(mc, data.frame(position = chr.len[chr],
                               frequency = 0,
                               group = "nGVHD"))
  }
  mc$CHROM <- chr
  # mc$group[which(mc$position %in% known_MiHA_chr$POS)] <- "MiHA"
  
  ### MiHA SNPs by chromosome
  if(length(known_MiHA_chr$POS) >= 1){
    index <- which(mc$position %in% known_MiHA_chr$POS)
    if(length(index) > 0){
      MiHA_SNP <- lapply(1:length(index), function(x) cbind(mc[index[x], ], known_MiHA_chr[which(known_MiHA_chr$POS == mc$position[index[x]]), ]))
      MiHA_SNP <- do.call(rbind.data.frame, MiHA_SNP)

      names(MiHA_SNP)[3] <- "Group"

      all_MiHA_SNP <- rbind(all_MiHA_SNP, MiHA_SNP)
      # print(
      #   ggplot(MiHA_SNP, aes(x = MiHAs, ymax = frequency, ymin = 0, color = Group)) +
      #     geom_linerange(size = 3) +
      #     geom_hline(yintercept = 0) +
      #     ggtitle(paste0(chr, " - MiHAs")) +
      #     theme_bw() +
      #     xlab("") +
      # scale_fill_discrete(name = "Group")
      # )
    }
  }

}

ggplot(all_MiHA_SNP, aes(x = MiHAs, ymax = frequency, ymin = 0, color = Group)) +
  geom_linerange(size = 3) +
  geom_hline(yintercept = 0) +
  ggtitle("Known MiHAs") +
  theme_bw() +
  xlab("") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 + 
                                     theme(axis.text.x = element_text(angle = 90, hjust = -1))))
all_MiHA_SNP$POS2 <- all_MiHA_SNP$POS
GR_all_MiHAs <- makeGRangesFromDataFrame(all_MiHA_SNP, keep.extra.columns = T, 
                                                   ignore.strand = T, seqnames.field = "CHROM",
                                                   start.field = "POS", end.field = "POS")

seqlengths(GR_all_MiHAs) <- chr.len[names(seqlengths(GR_all_MiHAs))]


# data("CRC", package = "biovizBase")
hg38sub_dataframe <- data.frame(CHROM = names(chr.len),
                                POS = chr.len, 
                                stringsAsFactors = F) 
hg38sub_dataframe$POS1 <- 1
hg38sub <- makeGRangesFromDataFrame(hg38sub_dataframe, keep.extra.columns = F,
                                    ignore.strand = T, seqnames.field = "CHROM",
                                    start.field = "POS1", end.field = "POS") 
seqlengths(hg38sub) <- chr.len[names(seqlengths(hg38sub))]
seqinfo(hg38sub)

ggplot() + layout_circle(GR_all_MiHAs, geom = "ideo", fill = "gray70", radius = 7, trackWidth = 3) + 
  layout_circle(geom = "point", color = "red", radius = 14,
                trackWidth = 3, grid = TRUE, aes(y = frequency))

ggplot() + layout_circle(hg38sub, geom = "ideo", fill = "gray70")  +#Ideogram
  layout_circle(GR_all_MiHAs, aes(x=POS2, y = frequency), 
                geom = "point", trackWidth = 3, radius = 5, grid = TRUE) +#somatic mutation
  # scale_colour_manual(values = alpha(c("#D55E00", "#0072B2"), .2)) +
  # scale_fill_manual(values = alpha(c("#D55E00", "#0072B2"), .2)) + 
  layout_circle(hg38sub, geom = "scale", radius = 15, size = 0.5, scale.type = "sci") +#Scale
  layout_circle(hg38sub, geom = "text", aes(label = seqnames), radius = 20, vjust = 4, size = 3)#label

################################
# quantsmooth
##############
# source("https://bioconductor.org/biocLite.R")
# # biocLite("BiocUpgrade")
# biocLite("quantsmooth")


################################
# chromPlot
#############
# source("https://bioconductor.org/biocLite.R")
# biocLite("chromPlot")
library("chromPlot")
data(hg_gap)
head(hg_gap)
chromPlot(gaps=hg_gap)
library(GenomicFeatures)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txgr <- transcripts(txdb)
chromPlot(gaps=hg_gap, annot1=txgr)

data_file7 <- system.file("extdata", "monocitosDEChr19-21.txt", package = "chromPlot")
monocytes  <- read.table(data_file7, sep="\t", header=TRUE, stringsAsFactors=FALSE)

data_file1 <- system.file("extdata", "hg19_refGeneChr19-21.txt", package = "chromPlot")
refGeneHg  <- read.table(data_file1, sep="\t", header=TRUE, stringsAsFactors=FALSE)
refGeneHg$Colors <- "red"
chromPlot(gaps=hg_gap, bands=hg_cytoBandIdeo, annot1=refGeneHg, 
          segment=monocytes, chrSide=c(-1,1,1,1,1,1,1,1), figCols=3, chr=c(19:21))