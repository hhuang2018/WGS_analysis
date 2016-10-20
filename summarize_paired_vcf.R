source('util.R', echo = FALSE)

load("../Data/GRCh38_gene_list.RData")

if (!require("vcfR", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("vcfR", dependencies = TRUE)
  library("vcfR", verbose=F, warn.conflicts =F)
}

vcf_fp <- "../"
vcf_file <- "test_stats_DtoR.vcf"

vcf_info <- read.vcfR(paste0(vcf_fp, vcf_file), verbose = FALSE)

if(is.null(colnames(vcf_info@gt))){
  if(dim(vcf_info@gt)[2] == 2){
    colnames(vcf_info@gt) <- c("FORMAT", "SAMPLE")
  }else{
    colnames(vcf_info@gt) <- c("FORMAT", sapply(1:(dim(vcf_info@gt)[2]-1), function(x) paste0("SAMPLE", x)))
  }
}

vcf_gt <- extract.gt(vcf_info, element = "GT")
CHROMs <- unique(vcf_info@fix[,1])
num_CHROMs <- length(CHROMs)

summary_table <- data.frame(CHROM = character(num_CHROMs),
                             NumVariants = numeric(num_CHROMs),
                             NumGenes = character(num_CHROMs),
#                              ALT = character(num_rows),
#                              GT = numeric(num_rows),
#                              groupID = character(num_rows),
                             stringsAsFactors = FALSE)
Gene_list <- list()

for(chr in 1:num_CHROMs){
  ptm <- proc.time() 
  
  chrom <- CHROMs[chr]
  
  chr_index <- which(vcf_info@fix[, 1] == chrom)
  vcf_gt <- extract.gt(vcf_info, element = "GT")
  vcf_gt <- vcf_gt[chr_index, ]
  
  num_rows <- length(chr_index)

  summary_table$CHROM[chr] <- chrom
  summary_table$NumVariants[chr] <- num_rows
  
  #######
  # Check all the genes in the list
  #######
  GRCh38_gene_list_chr <- GRCh38_gene_list[which(GRCh38_gene_list$chrom == chrom), ]
  gene_names <- unique(GRCh38_gene_list_chr$GeneName)
  num_genes <- length(gene_names)
  
  Gene_summary <- data.frame(GeneName = vector(mode = "character", length = num_genes), 
                             total_num = vector(mode = "integer", length = num_genes),
                             intron_variants_num =  vector(mode = "integer", length = num_genes),
                             exon_variants_num =  vector(mode = "integer", length = num_genes),
                             utr_variants_num = vector(mode = "integer", length = num_genes),
                             promoter_variants_num = vector(mode = "integer", length = num_genes),
                             stringsAsFactors = F)
  
  geneVariants <- as.data.frame(vcf_info@fix[vcf_info@fix[, 1] == chrom, ], stringsAsFactors = FALSE)
  names(geneVariants) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")
  geneVariants$POS <- as.integer(geneVariants$POS)
  
  counter <- 0
  for (id in 1:num_genes){
    
    GeneName <- as.character(gene_names[id])
    geneInfo <- GRCh38_gene_list_chr[which(GRCh38_gene_list_chr$GeneName == GeneName)[1], ]
    
    gene_region_index <- which((geneVariants$POS >= (geneInfo$txStart-999)) & (geneVariants$POS <= (geneInfo$txEnd+999)))
    if(length(gene_region_index) > 0){
      counter <- counter + 1
      gene_chr <- geneVariants[gene_region_index, ]
      gene_stats <- gene_variant_stats(geneInfo, gene_chr)
      
      Gene_summary$GeneName[counter] <- GeneName
      Gene_summary[counter, 2:6] <- gene_stats$variants_stats[1:5]
    }
    
  }
  
  summary_table$NumGenes[chr] <- counter
  
  Gene_summary <- Gene_summary[-which(Gene_summary$GeneName == ""), ]
  
  Gene_list <- append(Gene_list, list(Gene_summary))
  
  proc.time() - ptm
}
save(summary_table, Gene_list, file = "../Output/paired_VCF_summary/test_stats_DtoR.RData")

load("../Output/paired_VCF_summary/test_stats_DtoR.RData")
