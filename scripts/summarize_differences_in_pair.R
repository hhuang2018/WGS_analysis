# library(vcfR)
library(ggplot2)
library(reshape2)
# source('util.R', echo = FALSE)

summary_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/Allele_summary_by_chrome/"
#summary_file_dir <- "../Output/Chrom_allele_summary/"

output_dir <- "../Output/difference_summary/"

all_files <- list.files(summary_file_dir)

aGVHD_summary_files <- all_files[grepl("a_", all_files)]
# all_vcf_files <- all_vcf_files[grepl("a_", all_vcf_files)]
nGVHD_summary_files <- all_files[grepl("n_", all_files)]

pdf_file_name <- "aGVHD_v_nGVHD_all_pairs.pdf"
pdf(file = paste0(output_dir, pdf_file_name))

for(chr in 1:22){
  cat("Chromosome ", chr, "\n")
  
  ptm <- proc.time()
  
  chrom <- paste0("chr", chr)
  
#   aGVHD_summary_files <- list.files(summary_file_dir, pattern = paste0("a_*", chrom, ".RData"))
#   # all_vcf_files <- all_vcf_files[grepl("a_", all_vcf_files)]
#   aGVHD_summary_files <- aGVHD_summary_files[grepl(chrom,  aGVHD_summary_files)]
#   
#   nGVHD_summary_files <- list.files(summary_file_dir, pattern = "n_*")
#   nGVHD_summary_files <- nGVHD_summary_files[grepl(chrom,  nGVHD_summary_files)]
  aGVHD_chr_summary_files <- aGVHD_summary_files[grepl(chrom,  aGVHD_summary_files)]
  
  nGVHD_chr_summary_files <- nGVHD_summary_files[grepl(chrom,  nGVHD_summary_files)]
  
  aGVHD_num <- length(aGVHD_chr_summary_files)
  nGVHD_num <- length(nGVHD_chr_summary_files)
  
  if(aGVHD_num == 0 || nGVHD_num == 0) break

  all_chr_summary <- matrix(data = 0, nrow = aGVHD_num+nGVHD_num, ncol = 5)
  colnames(all_chr_summary) <- c("0", "1", "2", "All", "Group")
  all_chr_summary <- as.data.frame(all_chr_summary)
  
  for(id in 1:aGVHD_num){
    
    load(paste0(summary_file_dir, aGVHD_chr_summary_files[id]))
    
    all_chr_summary[id, 1:3] <- table(genotype_table$GT)
    all_chr_summary$All[id] <- sum(all_chr_summary[id, 1:3])
    all_chr_summary$Group[id] <- "aGVHD"
  } 
  
  for(jd in (1+aGVHD_num):(nGVHD_num+aGVHD_num)){
    
    load(paste0(summary_file_dir, nGVHD_chr_summary_files[jd - aGVHD_num]))
    
    all_chr_summary[jd, 1:3] <- table(genotype_table$GT)
    all_chr_summary$All[jd] <- sum(all_chr_summary[jd, 1:3])
    all_chr_summary$Group[jd] <- "nGVHD"
  } 
  
  print(t.test(all_chr_summary$`0`[1:aGVHD_num], all_chr_summary$`0`[(1+aGVHD_num):(nGVHD_num+aGVHD_num)], alternative = "two.sided"))
  print(t.test(all_chr_summary$`1`[1:aGVHD_num], all_chr_summary$`1`[(1+aGVHD_num):(nGVHD_num+aGVHD_num)], alternative = "two.sided"))
  print(t.test(all_chr_summary$`2`[1:aGVHD_num], all_chr_summary$`2`[(1+aGVHD_num):(nGVHD_num+aGVHD_num)], alternative = "two.sided"))
  print(t.test(all_chr_summary$All[1:aGVHD_num], all_chr_summary$All[(1+aGVHD_num):(nGVHD_num+aGVHD_num)], alternative = "two.sided"))
  
  reshaped_summary <- melt(all_chr_summary, id.vars = "Group")
  colnames(reshaped_summary) <- c("Group", "AlleleNum", "Counts")
  p <- ggplot(reshaped_summary, aes(x = factor(AlleleNum), y = Counts))
  print(p + geom_boxplot(aes(fill = Group)) + 
    ggtitle(paste0("Number of different variants within each pair\n (Chromosome ", chr,")")) +
    labs(x="Number of different alleles at each site", y = "Counts of sites"))

  print(proc.time() - ptm)
  save(all_chr_summary, file = paste0(output_dir, "aGVHD_nGVHD_all_pairs_", chrom,".RData"))
}
dev.off()

