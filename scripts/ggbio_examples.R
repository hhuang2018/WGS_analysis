library("ggbio")
# load(system.file("data", "hg19IdeogramCyto.rda", package="biovizBase", mustWork=TRUE))
# library(BSgenome.Hsapiens.UCSC.hg38)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(GenomicRanges)
library(biovizBase)

source("https://bioconductor.org/biocLite.R")
biocLite("Homo.sapiens") ## hg19
library(Homo.sapiens)
class(Homo.sapiens)

## gene plot
data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand = TRUE)

p.txdb <- autoplot(Homo.sapiens, which  = wh)
p.txdb

autoplot(Homo.sapiens, which  = wh, label.color = "black", color = "brown",
         fill = "brown")
autoplot(Homo.sapiens, which = wh, gap.geom = "chevron")

autoplot(Homo.sapiens, which = wh, stat = "reduce")

columns(Homo.sapiens)

autoplot(Homo.sapiens, which  = wh, columns = c("TXNAME", "GO"), names.expr = "TXNAME::GO")

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr.txdb <- crunch(txdb, which = wh)
colnames(values(gr.txdb))[4] <- "model"
grl <- split(gr.txdb, gr.txdb$tx_id)
names(grl) <- sample(LETTERS, size = length(grl), replace = TRUE)
grl

autoplot(grl, aes(type = model))
ggplot() + geom_alignment(grl, type = "model")

## circular plot
data("CRC", package = "biovizBase")
head(hg19sub)
autoplot(hg19sub, layout = "circle", fill  = "gray70")

p <- ggbio() + circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "scale", size = 2) +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
p

p <- ggbio(trackWidth = 10, buffer = 0, radius = 10) + circle(hg19sub, geom = "ideo", fill = "gray70") +
  circle(hg19sub, geom = "scale", size = 2) +
  circle(hg19sub, geom = "text", aes(label = seqnames), vjust = 0, size = 3)
p

grl <- split(crc.gr, values(crc.gr)$individual) ## need "unit", load grid
library(grid)
crc.lst <- lapply(grl, function(gr.cur){
  print(unique(as.character(values(gr.cur)$individual)))
  cols <- RColorBrewer::brewer.pal(3, "Set2")[2:1]
  names(cols) <- c("interchromosomal", "intrachromosomal")
  p <- ggbio() + circle(gr.cur, geom = "link", linked.to = "to.gr",
                        aes(color = rearrangements)) +
    circle(hg19sub, geom = "ideo",
           color = "gray70", fill = "gray70") +
    scale_color_manual(values = cols)  +
    labs(title = (unique(values(gr.cur)$individual))) +
    theme(plot.margin = unit(rep(0, 4), "lines"))
}
)
arrangeGrobByParsingLegend(crc.lst, widths = c(4, 1), legend.idx = 1, ncol = 3)


# stacked karyogram
data(ideoCyto, package = "biovizBase")
autoplot(seqinfo(ideoCyto$hg19), layout = "karyogram")
biovizBase::isIdeogram(ideoCyto$hg19)
autoplot(ideoCyto$hg19, layout = "karyogram", cytoband = TRUE)

