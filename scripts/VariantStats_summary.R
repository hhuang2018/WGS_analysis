### Varaint Stats summary

stats_dir <- "../Output/VariantStats/"

all_stats_files <- list.files(stats_dir, pattern = "\\.RData$")

num_files <- length(all_stats_files)
# a_exon_genes_donor <- character(length = 0)
# a_exon_genes_recipient <- character(length = 0)
# n_exon_genes_donor <- character(length = 0)
# n_exon_genes_recipient <- character(length = 0)

aGVHD_exon_genes_recipient <- matrix(data = "", nrow = 102, ncol = 100 *22)
aGVHD_exon_genes_donor <- matrix(data = "", nrow = 102, ncol = 100 *22)
nGVHD_exon_genes_recipient <- matrix(data = "", nrow = 103, ncol = 100 *22)
nGVHD_exon_genes_donor <- matrix(data = "", nrow = 103, ncol = 100 * 22)

aGVHD_intron_genes_recipient <- matrix(data = "", nrow = 102, ncol = 100 * 22 )
aGVHD_intron_genes_donor <- matrix(data = "", nrow = 102, ncol = 100 *22)
nGVHD_intron_genes_recipient <- matrix(data = "", nrow = 103, ncol = 100 * 22)
nGVHD_intron_genes_donor <- matrix(data = "", nrow = 103, ncol = 100 *22)

a_counter <- 0
n_counter <- 0
for(id in 1:num_files){
  cat("File #", id, "\n")
  eval(parse(text = paste0("load(\"", stats_dir, all_stats_files[id], "\")"))) 
  groupType <- unlist(strsplit(all_stats_files[id], "_"))[1]
  
  if(groupType == "a"){
    a_counter <- a_counter + 1
    for (chr in 1:22){
      cat("  --> Chromosome #", chr, "\n")
      eval(parse(text = paste0("chrStats <- ChromosomeStats[[\"chr", chr, "\"]]")))
      aGVHD_exon_genes_donor[a_counter, ((chr-1)*100+1):(chr*100)] <- as.character(chrStats$GeneName[order(chrStats$exon_diff, decreasing = T)[1:100]])
      aGVHD_exon_genes_recipient[a_counter, ((chr-1)*100+1):(chr*100)] <- as.character(chrStats$GeneName[order(-chrStats$exon_diff, decreasing = T)[1:100]])
      #cat(length(aa$GeneName[which(chrStats$exon_diff < 0)]), "\t")
      
      aGVHD_intron_genes_donor[a_counter, ((chr-1)*100+1):(chr*100)] <- as.character(chrStats$GeneName[order(chrStats$intron_diff, decreasing = T)[1:100]])
      aGVHD_intron_genes_recipient[a_counter, ((chr-1)*100+1):(chr*100)] <- as.character(chrStats$GeneName[order(-chrStats$intron_diff, decreasing = T)[1:100]])
      
    }
  }else{
    n_counter <- n_counter + 1
    for (chr in 1:22){
      eval(parse(text = paste0("chrStats <- ChromosomeStats[[\"chr", chr, "\"]]")))
      nGVHD_exon_genes_donor[n_counter, ((chr-1)*100+1):(chr*100)] <- as.character(chrStats$GeneName[order(chrStats$exon_diff, decreasing = T)[1:100]])
      nGVHD_exon_genes_recipient[n_counter, ((chr-1)*100+1):(chr*100)] <- as.character(chrStats$GeneName[order(-chrStats$exon_diff, decreasing = T)[1:100]])
      #cat(length(aa$GeneName[which(chrStats$exon_diff < 0)]), "\t")
      
      nGVHD_intron_genes_donor[n_counter, ((chr-1)*100+1):(chr*100)] <- as.character(chrStats$GeneName[order(chrStats$intron_diff, decreasing = T)[1:100]])
      nGVHD_intron_genes_recipient[n_counter, ((chr-1)*100+1):(chr*100)] <- as.character(chrStats$GeneName[order(-chrStats$intron_diff, decreasing = T)[1:100]])
      # cat(length(aa$GeneName[which(chrStats$exon_diff < 0)]), "\t")
      
    }
  }
}

nGVHD_freq_exon_genes_donor <- matrix(data = "", nrow = 22, ncol = 100)
nGVHD_freq_exon_genes_recipient <- matrix(data = "", nrow = 22, ncol = 100)

nGVHD_freq_intron_genes_donor <- matrix(data = "", nrow = 22, ncol = 100)
nGVHD_freq_intron_genes_recipient <- matrix(data = "", nrow = 22, ncol = 100)

for(chr in 1:22){
  freq_nGVHD_exon_genes_donor <- table(nGVHD_exon_genes_donor[, ((chr-1)*100+1):(chr*100)])
  nGVHD_freq_exon_genes_donor[chr, ] <- names(freq_nGVHD_exon_genes_donor[order(freq_nGVHD_exon_genes_donor, decreasing = T)[1:100]])
  
  freq_nGVHD_exon_genes_recipient <- table(nGVHD_exon_genes_recipient[, ((chr-1)*100+1):(chr*100)])
  nGVHD_freq_exon_genes_recipient[chr, ] <- names(freq_nGVHD_exon_genes_recipient[order(freq_nGVHD_exon_genes_recipient, decreasing = T)[1:100]])
  
  freq_nGVHD_intron_genes_donor <- table(nGVHD_intron_genes_donor[, ((chr-1)*100+1):(chr*100)])
  nGVHD_freq_intron_genes_donor[chr, ] <- names(freq_nGVHD_intron_genes_donor[order(freq_nGVHD_intron_genes_donor, decreasing = T)[1:100]])
  
  freq_nGVHD_intron_genes_recipient <- table(nGVHD_intron_genes_recipient[, ((chr-1)*100+1):(chr*100)])
  nGVHD_freq_intron_genes_recipient[chr, ] <- names(freq_nGVHD_intron_genes_recipient[order(freq_nGVHD_intron_genes_recipient, decreasing = T)[1:100]])
  
}

write.csv(nGVHD_freq_exon_genes_donor, file = "../Output/nGVHD_exon_genes_donor.csv")
write.csv(nGVHD_freq_exon_genes_recipient, file = "../Output/nGVHD_exon_genes_recipient.csv")
write.csv(nGVHD_freq_intron_genes_donor, file = "../Output/nGVHD_intron_genes_donor.csv")
write.csv(nGVHD_freq_intron_genes_recipient, file = "../Output/nGVHD_intron_genes_recipient.csv")


aGVHD_freq_exon_genes_donor <- matrix(data = "", nrow = 22, ncol = 100)
aGVHD_freq_exon_genes_recipient <- matrix(data = "", nrow = 22, ncol = 100)

aGVHD_freq_intron_genes_donor <- matrix(data = "", nrow = 22, ncol = 100)
aGVHD_freq_intron_genes_recipient <- matrix(data = "", nrow = 22, ncol = 100)

for(chr in 1:22){
  freq_aGVHD_exon_genes_donor <- table(aGVHD_exon_genes_donor[1:a_counter, ((chr-1)*100+1):(chr*100)])
  aGVHD_freq_exon_genes_donor[chr, ] <- names(freq_aGVHD_exon_genes_donor[order(freq_aGVHD_exon_genes_donor, decreasing = T)[1:100]])
  
  freq_aGVHD_exon_genes_recipient <- table(aGVHD_exon_genes_recipient[1:a_counter, ((chr-1)*100+1):(chr*100)])
  aGVHD_freq_exon_genes_recipient[chr, ] <- names(freq_aGVHD_exon_genes_recipient[order(freq_aGVHD_exon_genes_recipient, decreasing = T)[1:100]])
  
  freq_aGVHD_intron_genes_donor <- table(aGVHD_intron_genes_donor[1:a_counter, ((chr-1)*100+1):(chr*100)])
  aGVHD_freq_intron_genes_donor[chr, ] <- names(freq_aGVHD_intron_genes_donor[order(freq_aGVHD_intron_genes_donor, decreasing = T)[1:100]])
  
  freq_aGVHD_intron_genes_recipient <- table(aGVHD_intron_genes_recipient[1:a_counter, ((chr-1)*100+1):(chr*100)])
  aGVHD_freq_intron_genes_recipient[chr, ] <- names(freq_aGVHD_intron_genes_recipient[order(freq_aGVHD_intron_genes_recipient, decreasing = T)[1:100]])
  
}

write.csv(aGVHD_freq_exon_genes_donor, file = "../Output/aGVHD_exon_genes_donor.csv")
write.csv(aGVHD_freq_exon_genes_recipient, file = "../Output/aGVHD_exon_genes_recipient.csv")
write.csv(aGVHD_freq_intron_genes_donor, file = "../Output/aGVHD_intron_genes_donor.csv")
write.csv(aGVHD_freq_intron_genes_recipient, file = "../Output/aGVHD_intron_genes_recipient.csv")