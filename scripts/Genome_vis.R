lolliplot <- function (x, y = NULL, z = NULL, fillCol = NULL, labelCol = NULL, 
          txtAngle = 45, txtSize = 5, pntSize = 4, proteinColour = "#999999", 
          obsA.rep.fact = 5000, obsA.rep.dist.lmt = 500, obsA.attr.fact = 0.1, 
          obsA.adj.max = 0.1, obsA.adj.lmt = 0.5, obsA.iter.max = 50000, 
          obsB.rep.fact = 5000, obsB.rep.dist.lmt = 500, obsB.attr.fact = 0.1, 
          obsB.adj.max = 0.1, obsB.adj.lmt = 0.5, obsB.iter.max = 50000, 
          sideChain = FALSE, species = "hsapiens", maxLolliStack = NULL, 
          plotLayer = NULL, paletteA = NULL, paletteB = NULL, host = "www.ensembl.org", 
          out = "plot") 
{
  input <- lolliplot_qual(x, y, z)
  x <- input[[1]]
  y <- input[[2]]
  z <- input[[3]]
  transcriptID <- as.character(x$transcript_name[1])
  if (!is.null(y)) {
    y <- y[y$transcript_name == transcriptID, ]
  }
  gene <- as.character(x$gene[1])
  result <- lolliplot_transcriptID2codingSeq(transcriptID, 
                                             species = species, host = host)
  codingSeq <- result$coding
  cdsLen <- result$cds_length
  residueSeq <- lolliplot_DNAconv(codingSeq, to = "residue")
  if (sideChain == TRUE) {
    sidechain <- lolliplot_DNAconv(codingSeq, to = "sidechain")
    AAsequence <- cbind(sidechain, residueSeq)
    AAsequence <- as.data.frame(AAsequence)
    AAsequence$coord <- seq(from = 1, to = nrow(AAsequence))
  }
  else {
    AAsequence <- NULL
  }
  if (any(residueSeq %in% c("OPAL", "OCHRE", "AMBER"))) {
    stopRes <- c("OPAL", "OCHRE", "AMBER")
    residueSeq <- residueSeq[-which(residueSeq %in% stopRes)]
    if (!is.null(AAsequence)) {
      AAsequence <- AAsequence[-which(AAsequence$residueSeq %in% 
                                        stopRes), ]
    }
  }
  proteinLength <- length(residueSeq)
  if (!is.null(z)) {
    geneData <- lolliplot_constructGene(gene, z, proteinLength)
  }
  else {
    protein_domain <- lolliplot_fetchDomain(transcriptID, 
                                            species = species, host = host)
    geneData <- lolliplot_constructGene(gene, protein_domain, 
                                        proteinLength)
  }
  observed_mutation <- lolliplot_mutationObs(x, "top", fillCol, 
                                             labelCol, obsA.rep.fact, obsA.rep.dist.lmt, obsA.attr.fact, 
                                             obsA.adj.max, obsA.adj.lmt, obsA.iter.max)
  observed_mutation <- lolliplot_reduceLolli(observed_mutation, 
                                             maxLolliStack)
  if (!is.null(y)) {
    observed_mutation2 <- lolliplot_mutationObs(y, "bottom", 
                                                fillCol, labelCol, obsB.rep.fact, obsB.rep.dist.lmt, 
                                                obsB.attr.fact, obsB.adj.max, obsB.adj.lmt, obsB.iter.max)
    observed_mutation2 <- lolliplot_reduceLolli(observed_mutation2, 
                                                maxLolliStack)
  }
  else {
    observed_mutation2 <- NULL
  }
  plot <- lolliplot_buildMain(geneData, length, observed_mutation, 
                              observed_mutation2, fillCol, labelCol, txtAngle, txtSize, 
                              pntSize, proteinColour, AAsequence, plot_sidechain = sideChain, 
                              layers = plotLayer, paletteA = paletteA, paletteB = paletteB)
  dataOut <- list(gene = geneData, observed_mutation = observed_mutation, 
                  observed_mutation2 = observed_mutation2)
  output <- multi_selectOut(data = dataOut, plot = plot, out = out)
  return(output)
}


###########
## Plotting differentially methylated bases on an ideogram
###########

#' function for making ideogram for differential methylation values
#' requires methylKit, ggbio and GenomicRanges
#'
#' @example
#' library(BSgenome)
#' library("BSgenome.Hsapiens.UCSC.hg18")
#' chr.len = seqlengths(Hsapiens)  # get chromosome lengths
#' # remove X,Y,M and random chromosomes
#' chr.len = chr.len[grep("_|M|X|Y", names(chr.len), invert = T)] 
#' 
#' download.file("http://methylkit.googlecode.com/files/myDiff.rda", 
#'               destfile = "myDiff.rda")
#' load("myDiff.rda")
#' 
#' ideoDMC(myDiff, chrom.length = chr.len, difference = 25, qvalue = 0.01, 
#'        circos = TRUE, title = "test", hyper.col = "magenta", hypo.col = "green")

ideoDMC <- function(methylDiff.obj, chrom.length, difference = 25, 
                    qvalue = 0.01, circos = FALSE, title = "test", hyper.col = "magenta", 
                    hypo.col = "green") {
  if (!require("methylKit", quietly=TRUE, warn.conflicts = FALSE)) {
    library(devtools)
    install_github("al2na/methylKit", build_vignettes=FALSE, 
                   repos=BiocInstaller::biocinstallRepos(),
                   dependencies=TRUE)
    library("methylKit", verbose=F, warn.conflicts =F)
  }
  
  if (!require("optparse", quietly=TRUE, warn.conflicts = FALSE)) {
    install.packages("optparse", dependencies = TRUE)
    library("optparse", verbose=F, warn.conflicts =F)
  }
  
  if (!require("optparse", quietly=TRUE, warn.conflicts = FALSE)) {
    install.packages("optparse", dependencies = TRUE)
    library("optparse", verbose=F, warn.conflicts =F)
  }
  
  require(methylKit)
  require(GenomicRanges)
  require(ggbio)
  
  # chrom.length
  myIdeo <- GRanges(seqnames = names(chrom.length), ranges = IRanges(start = 1, 
                                                                     width = chrom.length))
  seqlevels(myIdeo) = names(chrom.length)
  seqlengths(myIdeo) = (chrom.length)
  
  
  hypo = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
                        type = "hypo")
  hyper = get.methylDiff(methylDiff.obj, difference = difference, qvalue = qvalue, 
                         type = "hyper")
  
  g.per = as(hyper, "GRanges")
  seqlevels(g.per, force=TRUE) = seqlevels(myIdeo)
  seqlengths(g.per)=(chrom.length)
  
  g.po = as(hypo, "GRanges")
  seqlevels(g.po, force=TRUE) = seqlevels(myIdeo)
  seqlengths(g.po)=(chrom.length)
  
  values(g.po)$id = "hypo"
  values(g.per)$id = "hyper"
  
  if (circos) {
    
    p <- ggplot() + layout_circle(myIdeo, geom = "ideo", fill = "gray70", 
                                  radius = 39, trackWidth = 2)
    
    
    p <- p + layout_circle(c(g.po, g.per), geom = "point", 
                           size = 1, aes(x = midpoint, 
                                         y = meth.diff, color = id), radius = 25, trackWidth = 30) +              
      scale_colour_manual(values = c(hyper.col, hypo.col))
    p + layout_circle(myIdeo, geom = "text", aes(label = seqnames), 
                      vjust = 0, radius = 55, trackWidth = 7) + labs(title = title)
    
  } else {
    
    p <- ggplot() + layout_karyogram(myIdeo, cytoband = FALSE)
    p + layout_karyogram(c(g.po, g.per), geom = "point", size = 1, 
                         aes(x = midpoint, 
                             y = meth.diff, color = id)) + scale_colour_manual(values = c(hyper.col, 
                                                                                          hypo.col)) + labs(title = title)
    # new alternative commented out
    #autoplot(c(g.po, g.per), layout = "karyogram", geom = "point", size = 0.65, 
    #aes(x = midpoint,y = meth.diff, color = id))  + scale_colour_manual(values = c(hyper.col, 
    #                                                                                        hypo.col)) + labs(title = title)
    
  }
}