source('stats.R', echo = FALSE)
source('util.R', echo = FALSE)
require(ggplot2)
stats_file_dir <- "../VariantStats/"
# Donor_file <- "n_197_D_62284237_annotated.vcf.gz"
# Recipient_file <- "n_197_R_1005928_annotated.vcf.gz"

all_RData_files <- list.files(stats_file_dir, pattern = "\\.RData$")

num_files <- length(all_RData_files)

ordered_chr_genes <- list()

for(id in 1:num_files){
  
  load(file = paste0(stats_file_dir, all_RData_files[id]))
  
  pdf_file_name <- gsub("_diff_variants.RData", "_chr_boxplot.pdf", all_RData_files[id])
  pdf(file = paste0(stats_file_dir, pdf_file_name))
  
  for(chr in 1:22){
    eval(parse(text = paste0("chr_stats <- ChromosomeStats[[\"chr", chr, "\"]]")))
    
#     plot.param <- reshape(chr_stats[, c(3, 5, 7, 9, 11, 13, 15)], dir = "long",
#                           idvar = "total_diff", 
#                           varying=c("intron_diff","exon_diff","utr_diff","promoter_diff", "insertion_diff", "deletion_diff"),
#                           v.names="y")
#     p.box <- ggplot(data = chr_stats[, c(3, 5, 7, 9, 11, 13, 15)])
#     p.box + geom_boxplot()
    # xlab_axis <- as.factor(colnames(chr_stats[, c(3, 5, 7, 9, 11, 13, 15)]))
    bp <- boxplot(chr_stats[, c(3, 5, 7, 9, 11, 13, 15)],
            main=paste0("Chromosome ", chr), 
            xlab="Variants type",
            ylab="Raw log-score", 
            axis = FALSE,
            axisnames = FALSE)
    
#     text(bp, labels=xlab_axis, par("usr")[3], 
#         srt=45, adj=1, xpd=TRUE)
#     axis(2, at=xlab_axis, labels=FALSE)
    # axis(2)
#     for(feat_region in colnames(chr_stats)[-1]){
#       ordered_gene_features <- sort_chr_genes(chr_stats, feat_region)
#     }
  } 
  dev.off()
  
}