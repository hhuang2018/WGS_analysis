normalized_VariantStats_by_Chr <- function(VCF_file_dir, Donor_file, Recipient_file, out.dir){
  source('util.R', echo = FALSE)
  load("../Data/GRCh38_gene_list.RData")
  
  library(vcfR)
  require(ggplot2)
  
  # VCF_file_dir <- "/mnt/scratch/hhuang/hli_vcf_annotated/"
  # Donor_file <- "n_197_D_62284237_annotated.vcf.gz"
  # Recipient_file <- "n_197_R_1005928_annotated.vcf.gz"
  
  # donor's VCF - Chromosome 
  vcf_file_D <- paste0(VCF_file_dir, Donor_file)
  vcf_D <- read.vcfR(vcf_file_D, verbose = FALSE)
  
  # recipient's VCF - Chromosome 
  vcf_file_R <- paste0(VCF_file_dir, Recipient_file)
  vcf_R <- read.vcfR(vcf_file_R, verbose = FALSE)
  
  output_filename <- paste0(unlist(strsplit(Donor_file, "_"))[1], "_",
                            unlist(strsplit(Donor_file, "_"))[2], "_",
                            "R_D_diff_variants")
  #pdf(paste0("../Output/",output_filename,".pdf"), width = 15, height = 20)
  #pdf(paste0(output_filename,".pdf"), width = 15, height = 20)
  
  
  ChromosomeStats <- list()
  for(chr in 1:22){
    
    pdf(paste0(out.dir, output_filename,"chr_", chr, ".pdf"), width = 15, height = 20)
    
    cat("Chromosome ", chr, ": ")
    
    Chrom <- paste0("chr", chr)
    
    # donor's VCF - Chromosome 
    D_chr <- as.data.frame(vcf_D@fix[vcf_D@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
    D_chr$POS <- as.integer(D_chr$POS)
    
    # rm(vcf_D)
    
    # recipient's VCF - Chromosome 
    R_chr <- as.data.frame(vcf_R@fix[vcf_R@fix[,1] == Chrom, ], stringsAsFactors = FALSE)
    R_chr$POS <- as.integer(R_chr$POS)
    
    # rm(vcf_R)
    
    #######
    # Check all the genes in the list
    #######
    GRCh38_gene_list_chr <- GRCh38_gene_list[which(GRCh38_gene_list$chrom == Chrom), ]
    gene_names <- unique(GRCh38_gene_list_chr$GeneName)
    num_genes <- length(gene_names)
    
    diff_num <- data.frame(GeneName = gene_names,
                           total_all = vector(mode = "integer", length = num_genes),
                           total_diff = vector(mode = "integer", length = num_genes),
                           intron_all =  vector(mode = "integer", length = num_genes),
                           intron_diff =  vector(mode = "integer", length = num_genes),
                           exon_all =  vector(mode = "integer", length = num_genes),
                           exon_diff =  vector(mode = "integer", length = num_genes),
                           utr_all =  vector(mode = "integer", length = num_genes),
                           utr_diff =  vector(mode = "integer", length = num_genes),
                           promoter_all =  vector(mode = "integer", length = num_genes),
                           promoter_diff = vector(mode = "integer", length = num_genes), 
                           insertion_all =  vector(mode = "integer", length = num_genes),
                           insertion_diff = vector(mode = "integer", length = num_genes), 
                           deletion_all =  vector(mode = "integer", length = num_genes),
                           deletion_diff = vector(mode = "integer", length = num_genes))
    
    for (id in 1:num_genes){
      
      GeneName <- gene_names[id]
      geneInfo <- GRCh38_gene_list_chr[which(GRCh38_gene_list_chr$GeneName == GeneName)[1], ]
      
      D_stats <- gene_variant_stats(geneInfo, D_chr)
      R_stats <- gene_variant_stats(geneInfo, R_chr)
      
      full_gene_diff <- compare_variants(D_stats$all_variants, R_stats$all_variants, D_stats$gene_length)
      exon_diff <- compare_variants(D_stats$exon_variants, R_stats$exon_variants, D_stats$exon_length)
      intron_diff <- compare_variants(D_stats$intron_variants, R_stats$intron_variants, D_stats$intron_length)
      utr_diff <- compare_variants(D_stats$utr_variants, R_stats$utr_variants, D_stats$utr_length)
      promoter_diff <- compare_variants(D_stats$promoter_variants, R_stats$promoter_variants, D_stats$promoter_length)
      
      insertion_diff <- compare_variants(D_stats$insert_variants, R_stats$insert_variants, dim(D_stats$insert_variants)[1])
      deletion_diff <- compare_variants(D_stats$deletion_variants, R_stats$deletion_variants, dim(D_stats$deletion_variants)[1])
      
      if(!is.na(full_gene_diff$allPos_score)){ # Full gene region
        diff_num$total_all[id] <- full_gene_diff$log_allPos_score
        diff_num$total_diff[id] <- full_gene_diff$log_diffPos_score
      }else{
        diff_num$total_all[id] <- NA
        diff_num$total_diff[id] <- NA
      }
      
      if(!is.na(intron_diff$allPos_score)){ # intron regions
        diff_num$intron_all[id] <- intron_diff$log_allPos_score
        diff_num$intron_diff[id] <- intron_diff$log_diffPos_score
      } else{
        diff_num$intron_all[id] <- NA
        diff_num$intron_diff[id] <- NA
      }
      
      if(!is.na(exon_diff$allPos_score)){ # exon region
        diff_num$exon_all[id] <- exon_diff$log_allPos_score
        diff_num$exon_diff[id] <- exon_diff$log_diffPos_score
      } else{
        diff_num$exon_all[id] <- NA
        diff_num$exon_diff[id] <- NA
      }
      
      if(!is.na(utr_diff$allPos_score)){ # UTR regions
        diff_num$utr_all[id] <- utr_diff$log_allPos_score
        diff_num$utr_diff[id] <- utr_diff$log_diffPos_score
      } else{
        diff_num$utr_all[id] <- NA
        diff_num$utr_diff[id] <- NA
      }
      
      if(!is.na(promoter_diff$allPos_score)){ # Promoter regions
        diff_num$promoter_all[id] <- promoter_diff$log_allPos_score
        diff_num$promoter_diff[id] <- promoter_diff$log_diffPos_score
      } else{
        diff_num$promoter_all[id] <- NA
        diff_num$promoter_diff[id] <- NA
      }
      
      if(!is.na(insertion_diff$allPos_score)){ # insertion
        diff_num$insertion_all[id] <- insertion_diff$samePos_change 
        diff_num$insertion_diff[id] <- insertion_diff$diffPos_change
      } else{
        diff_num$insertion_all[id] <- NA
        diff_num$insertion_diff[id] <- NA
      }
      
      if(!is.na(deletion_diff$allPos_score)){ # deletion 
        diff_num$deletion_all[id] <- deletion_diff$samePos_change
        diff_num$deletion_diff[id] <- deletion_diff$diffPos_change
      } else{
        diff_num$deletion_all[id] <- NA
        diff_num$deletion_diff[id] <- NA
      }
    }
    
    chrVariants_stats <- plot.chrVariants(diff_num, "total", "diff") # full gene region
    multiplot(chrVariants_stats$p_RD_pair, chrVariants_stats$p_D_prime, chrVariants_stats$p_R_prime, cols = 1)
    
    chrVariants_stats <- plot.chrVariants(diff_num, "exon", "diff") # exon region
    multiplot(chrVariants_stats$p_RD_pair, chrVariants_stats$p_D_prime, chrVariants_stats$p_R_prime, cols = 1)
    
    chrVariants_stats <- plot.chrVariants(diff_num, "intron", "diff") # intron region
    multiplot(chrVariants_stats$p_RD_pair, chrVariants_stats$p_D_prime, chrVariants_stats$p_R_prime, cols = 1)
    
    chrVariants_stats <- plot.chrVariants(diff_num, "promoter", "diff") # promoter region
    multiplot(chrVariants_stats$p_RD_pair, chrVariants_stats$p_D_prime, chrVariants_stats$p_R_prime, cols = 1)
    
    chrVariants_stats <- plot.chrVariants(diff_num, "utr", "diff") # utr region
    multiplot(chrVariants_stats$p_RD_pair, chrVariants_stats$p_D_prime, chrVariants_stats$p_R_prime, cols = 1)
    
    chrVariants_stats <- plot.chrVariants(diff_num, "insertion", "diff") # insertion region
    multiplot(chrVariants_stats$p_RD_pair, chrVariants_stats$p_D_prime, chrVariants_stats$p_R_prime, cols = 1)
    
    chrVariants_stats <- plot.chrVariants(diff_num, "deletion", "diff") # deletion region
    multiplot(chrVariants_stats$p_RD_pair, chrVariants_stats$p_D_prime, chrVariants_stats$p_R_prime, cols = 1)
    
    #   ordered_total <- diff_num[order(diff_num$total_diff, decreasing = TRUE, na.last = TRUE), c(1,3)]
    #   ordered_total <- ordered_total[!is.na(ordered_total$total_diff), ]
    #   
    #   ordered_total$GeneName <- factor(ordered_total$GeneName, levels = factor(ordered_total$GeneName))
    #   p_total_diff <- ggplot(data = ordered_total, aes(x = GeneName, y = total_diff)) + 
    #     geom_point()  + 
    #     theme(text = element_text(size=3), axis.text.x = element_text(angle = 60, hjust = 1)) +
    #     # ggtitle("All genes with variants") +
    #     annotate("text", x = length(ordered_total$GeneName)/2, y = 11, size = 10, label = "All genes with variants") +
    #     ylab("Normalized log-scale score")
    #   
    #   D_prior_genes <- ordered_total[ordered_total$total_diff >0 ,]
    #   D_prior_genes$GeneName <- factor(D_prior_genes$GeneName, levels = factor(D_prior_genes$GeneName))
    #   pd_diff <- ggplot(data = D_prior_genes, aes(x = GeneName, y = total_diff)) + 
    #     geom_point()  + 
    #     theme(text = element_text(size=4.5), axis.text.x = element_text(angle = 60, hjust = 1)) +
    #     # ggtitle("Genes with donor-dominant variants") +
    #     annotate("text", x = length(D_prior_genes$GeneName)/2, y = 11, size = 10, label = "Genes with donor-dominant variants") +
    #     ylab("Log-scale score")
    #   
    #   R_prior_genes <- ordered_total[ordered_total$total_diff <0 ,]
    #   R_prior_genes$total_diff <- -R_prior_genes$total_diff
    #   R_prior_genes <- R_prior_genes[order(R_prior_genes$total_diff, decreasing = T, na.last = T), ]
    #   R_prior_genes$GeneName <- factor(R_prior_genes$GeneName, levels = factor(R_prior_genes$GeneName))
    #   pr_diff <- ggplot(data = R_prior_genes, aes(x = GeneName, y = total_diff)) + 
    #     geom_point()  + 
    #     theme(text = element_text(size=4.5), axis.text.x = element_text(angle = 60, hjust = 1)) +
    #     #ggtitle("Genes with Recipient-dominant variants") +
    #     annotate("text", x = length(R_prior_genes$GeneName)/2, y = 11, size = 10, label = "Genes with Recipient-dominant variants") +
    #     ylab("Log-scale score")
    
    #   plot(diff_num$total_diff, 
    #        main = paste0("Full Gene region (Chormosome", chr,")"),
    #        xlab = "Gene Index",
    #        ylab = "Normalized of different variants within the gene feature (log)")
    # 
    #   plot(diff_num$intron, 
    #        main = paste0("Gene Intron regions (Chormosome", chr,")"),
    #        xlab = "Gene Index",
    #        ylab = "Normalized of different variants within the gene feature")
    #   plot(diff_num$exon, 
    #        main = paste0("Gene Exon regions (Chormosome", chr,")"),
    #        xlab = "Gene Index",
    #        ylab = "Normalized of different variants within the gene feature")
    #   
    #   save(diff_num, file = paste0("../Output/",output_filename, "_chr", chr,".RData"))
    #   
    
    eval(parse(text=paste0("ChromosomeStats[[\"chr", chr, "\"]] <- diff_num")))
    dev.off()
    
    cat("Done!\n")
  }
  save(ChromosomeStats, file = paste0(out.dir, output_filename, ".RData"))
  # save(diff_num, file = "../Data/a_197_R_D_diff_variants.RData")
}