source('util.R', echo = FALSE)

load("../Data/GRCh38_gene_list.RData")

library(vcfR)

VCF_file_dir <- "/mnt/scratch/hhuang/hli_vcf_annotated/"
Donor_file <- "a_208_D_528085277_annotated.vcf.gz"
Recipient_file <- "a_208_R_1018111_annotated.vcf.gz"

output_filename <- paste0(unlist(strsplit(Donor_file, "_"))[1], "_",
                          unlist(strsplit(Donor_file, "_"))[2], "_",
                          "R_D_diff_variants")
pdf(paste0("../Output/",output_filename,".pdf"))

for(chr in 1:22){
  cat("Chromosome ", chr)
  Chrom <- paste0("chr", chr)
  
  # donor's VCF - Chromosome 
  vcf_file_D <- paste0(VCF_file_dir, Donor_file)
  vcf_D <- read.vcfR(vcf_file_D, verbose = FALSE)
  
  # donor_table <- as.data.frame(vcf_D@fix, stringsAsFactors = FALSE)
  
  D_chr <- as.data.frame(vcf_D@fix[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  D_chr$POS <- as.integer(D_chr$POS)
  
  rm(vcf_D)
  
  # recipient's VCF - Chromosome 
  vcf_file_R <- paste0(VCF_file_dir, Recipient_file)
  vcf_R <- read.vcfR(vcf_file_R, verbose = FALSE)
  
  R_chr <- as.data.frame(vcf_R@fix[vcf_R@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
  R_chr$POS <- as.integer(R_chr$POS)
  
  rm(vcf_R)
  
  #######
  # Check all the genes in the list
  #######
  GRCh38_gene_list_chr <- GRCh38_gene_list[which(GRCh38_gene_list$chrom == Chrom), ]
  gene_names <- unique(GRCh38_gene_list_chr$GeneName)
  num_genes <- length(gene_names)
  
  diff_num <- data.frame(total = vector(mode = "integer", length = num_genes),
                         intron =  vector(mode = "integer", length = num_genes),
                         exon =  vector(mode = "integer", length = num_genes))
  
  for (id in 1:num_genes){
    
    GeneName <- gene_names[id]
    geneInfo <- GRCh38_gene_list_chr[which(GRCh38_gene_list_chr$GeneName == GeneName)[1], ]
    
    D_stats <- gene_variant_stats(geneInfo, D_chr)
    R_stats <- gene_variant_stats(geneInfo, R_chr)
    
    D_same_variants <- D_stats$all_variants[which(D_stats$all_variants$POS %in% R_stats$all_variants$POS), ]
    R_same_variants <- R_stats$all_variants[which(R_stats$all_variants$POS %in% D_stats$all_variants$POS), ]
    
    #   if(!identical(D_same_variants$ALT, R_same_variants$ALT)){
    #     
    #   }
    
    D_diff_variants <- D_stats$all_variants[-which(D_stats$all_variants$POS %in% R_stats$all_variants$POS), ]
    R_diff_variants <- R_stats$all_variants[-which(R_stats$all_variants$POS %in% D_stats$all_variants$POS), ]
    
    diff_num$total[id] <- dim(D_diff_variants)[1] + dim(R_diff_variants)[1]
    diff_num$intron[id] <- length(which(grepl("Intron", D_diff_variants$ID))) + length(which(grepl("Intron", R_diff_variants$ID)))
    diff_num$exon[id] <- length(which(grepl("Exon", D_diff_variants$ID))) + length(which(grepl("Exon", R_diff_variants$ID)))      
    
  }
  
  plot(diff_num$total, 
       main = paste0("Full Gene region (Chormosome", chr,")"),
       xlab = "Gene Index",
       ylab = "Number of different variants")
  plot(diff_num$intron, 
       main = paste0("Gene Intron regions (Chormosome", chr,")"),
       xlab = "Gene Index",
       ylab = "Number of different variants")
  plot(diff_num$exon, 
       main = paste0("Gene Exon regions (Chormosome", chr,")"),
       xlab = "Gene Index",
       ylab = "Number of different variants")
  
  save(diff_num, file = paste0("../Output/",output_filename, "_chr", chr,".RData"))
}

dev.off()
# save(diff_num, file = "../Data/a_197_R_D_diff_variants.RData")
