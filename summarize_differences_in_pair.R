# library(vcfR)
library(ggplot2)
library(reshape2)
# source('util.R', echo = FALSE)

#summary_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"
summary_file_dir <- "../Output/Chrom_allele_summary/"

output_dir <- "../Output/difference_summary/"

pdf_file_name <- "aGVHD_v_nGVHD.pdf"
pdf(file = paste0(output_dir, pdf_file_name))

for(chr in 1:22){
  cat("Chromosome ", chr, "\n")
  
  ptm <- proc.time()
  
  chrom <- paste0("chr", chr)
  
  aGVHD_summary_files <- list.files(summary_file_dir, pattern = "a_")
  # all_vcf_files <- all_vcf_files[grepl("a_", all_vcf_files)]
  aGVHD_summary_files <- aGVHD_summary_files[grepl(chrom,  aGVHD_summary_files)]
  
  nGVHD_summary_files <- list.files(summary_file_dir, pattern = "n_")
  nGVHD_summary_files <- nGVHD_summary_files[grepl(chrom,  nGVHD_summary_files)]
  
  aGVHD_num <- length(aGVHD_summary_files)
  nGVHD_num <- length(nGVHD_summary_files)
  
  if(aGVHD_num == 0 || nGVHD_num == 0) break

  all_summary <- matrix(data = 0, nrow = aGVHD_num+nGVHD_num, ncol = 5)
  colnames(all_summary) <- c("0", "1", "2", "All", "Group")
  all_summary <- as.data.frame(all_summary)
  
  for(id in 1:aGVHD_num){
    
    load(paste0(summary_file_dir, aGVHD_summary_files[id]))
    
    all_summary[id, 1:3] <- table(genotype_table$GT)
    all_summary$All[id] <- sum(all_summary[id, 1:3])
    all_summary$Group[id] <- "aGVHD"
  } 
  
  for(jd in (1+aGVHD_num):(nGVHD_num+aGVHD_num)){
    
    load(paste0(summary_file_dir, nGVHD_summary_files[jd - aGVHD_num]))
    
    all_summary[jd, 1:3] <- table(genotype_table$GT)
    all_summary$All[jd] <- sum(all_summary[jd, 1:3])
    all_summary$Group[jd] <- "nGVHD"
  } 
  
  t.test(all_summary$`0`[1:aGVHD_num], all_summary$`0`[(1+aGVHD_num):(nGVHD_num+aGVHD_num)], alternative = "two.sided")
  t.test(all_summary$`1`[1:aGVHD_num], all_summary$`1`[(1+aGVHD_num):(nGVHD_num+aGVHD_num)], alternative = "two.sided")
  t.test(all_summary$`2`[1:aGVHD_num], all_summary$`2`[(1+aGVHD_num):(nGVHD_num+aGVHD_num)], alternative = "two.sided")
  t.test(all_summary$All[1:aGVHD_num], all_summary$All[(1+aGVHD_num):(nGVHD_num+aGVHD_num)], alternative = "two.sided")
  

  reshaped_summary <- melt(all_summary, id.vars = "Group")
  colnames(reshaped_summary) <- c("Group", "AlleleNum", "Counts")
  p <- ggplot(reshaped_summary, aes(x = factor(AlleleNum), y = Counts))
  p + geom_boxplot(aes(fill = Group)) + 
    ggtitle(paste0("Number of different variants within each pair\n (Chromosome ", chr,")")) +
    labs(x="Number of different alleles at each site", y = "Counts of sites")

  print(proc.time() - ptm)
}
dev.off()