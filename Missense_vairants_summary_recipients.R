source('util.R', echo = FALSE)

# load("../Data/GRCh38_gene_list.RData")

library(vcfR)

VCF_file_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants/"
output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants/stats/"

all_vcf_files <- list.files(VCF_file_dir, pattern = "\\.vcf.gz$")

file_names <- all_vcf_files[grepl("_R_", all_vcf_files)] # 216 donors;  240 recipients

num_files <- length(file_names)  
missense_stats <- data.frame(CHROM=character(0), POS=numeric(0), NumDiff = numeric(0), stringsAsFactors = F)
for(id in 1:num_files){
  
  vcf_file <- paste0(VCF_file_dir, file_names[id])
  vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  temp_info <- vcf_info@fix
  vcf_missense <- temp_info[temp_info[,"FILTER"] == "PASS", c("CHROM", "POS", "REF","ALT")]
  
  if(dim(missense_stats)[1] > 0){
    vcf_chrom_pos <- sapply(1:dim(vcf_missense)[1], function(x) paste0(vcf_missense[x, c("CHROM","POS")], collapse = ""))
    all_chrom_pos <- sapply(1:dim(missense_stats)[1], function(x) paste0(missense_stats[x, c("CHROM","POS")], collapse = ""))
    
    inter_chrom_pos <- intersect(vcf_chrom_pos, all_chrom_pos)
    exist_chrom_pos_id <- which(all_chrom_pos %in% inter_chrom_pos)
    new_chrom_pos_id <- which(!vcf_chrom_pos %in% inter_chrom_pos)
    
    if(length(exist_chrom_pos_id) >0){
      missense_stats$NumDiff[exist_chrom_pos_id] <- missense_stats$NumDiff[exist_chrom_pos_id] + 1
      if(length(new_chrom_pos_id) > 0){
        temp_stats <- data.frame(vcf_missense[new_chrom_pos_id, c("CHROM", "POS")], NumDiff = numeric(length(new_chrom_pos_id))+1, 
                                 stringsAsFactors = F)
        missense_stats <- rbind(missense_stats, temp_stats)
      }
    }else{
      if(length(new_chrom_pos_id) > 0){
        temp_stats <- data.frame(vcf_missense[new_chrom_pos_id, c("CHROM", "POS")], NumDiff = numeric(length(new_chrom_pos_id))+1, 
                                 stringsAsFactors = F)
        missense_stats <- rbind(missense_stats, temp_stats)
      }
    }
    
  } else {
    missense_stats <- rbind(missense_stats, cbind(vcf_missense[, c("CHROM", "POS")], NumDiff = numeric(dim(vcf_missense)[1])+1))
    missense_stats$CHROM <- as.character(missense_stats$CHROM)
    missense_stats$POS <- as.numeric(as.character(missense_stats$POS))
    missense_stats$NumDiff <- as.numeric(as.character(missense_stats$NumDiff))
  }
}
recipient_missense_stats <- missense_stats
save(recipient_missense_stats, file = paste0(output_dir, "recipient_missesense_stats.RData"))