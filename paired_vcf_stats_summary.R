
library(vcfR)
source('util.R', echo = FALSE)

vcf_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"
#vcf_file_dir <- "../paired_vcf/"
all_vcf_files <- list.files(vcf_file_dir, pattern = "\\.vcf.gz$")
output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/paired_vcf_Stats_Summary/"

num_files <- length(all_vcf_files)

# SKAT genotype code: 
#    0 - both alelles are the same 
#    1 - only one allele is the same
#    2 - both alleles are different

for(id in 1:num_files){
  
  vcf_file <- paste0(vcf_file_dir, all_vcf_files[id])
  vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  
  vcf_gt <- extract.gt(vcf_info, element = "GT")
  
  groupID  <- gsub(".vcf.gz", "", all_vcf_files[id])
    
  num_rows <- dim(vcf_gt)[1]
  genotype_table <- data.frame(CHROM = character(num_rows),
                               POS = numeric(num_rows),
                               REF = character(num_rows),
                               ALT = character(num_rows),
                               GT = numeric(num_rows),
                               groupID = character(num_rows),
                               stringsAsFactors = FALSE)
  genotype_table$CHROM <- vcf_info@fix[, 1]
  genotype_table$POS <- as.numeric(vcf_info@fix[, 2])
  genotype_table$REF <- vcf_info@fix[, 4]
  genotype_table$ALT <- vcf_info@fix[, 5]
  genotype_table$groupID <- rep(groupID, times = num_rows)
  
  cat("File ", id, ": ", all_vcf_files[id], "\n")
  # cat("[")
  ptm <- proc.time() 
  is_same_gt <- sapply(1:num_rows, function(x) vcf_gt[x, 1] == vcf_gt[x, 2])
  
  diff_gt_index <- which(!is_same_gt)
  genotype_table$GT[diff_gt_index] <- sapply(1:length(diff_gt_index), function(x) is.same.gt(vcf_gt[diff_gt_index[x], ]))
  
  proc.time() - ptm
  #genotype_table$GT[is_same_gt] <- 0
  
#   for(rind in 1:num_rows){
#     
# #     progress <- round(rind / num_rows *100, digits = 2)
# #     if(progress > 0 && (progress %% 10 == 0)) {
# #       cat("=")
# #     }
#     
# #     genotype_table$CHROM[rind] <- vcf_info@fix[rind, 1]
# #     genotype_table$POS[rind] <- as.numeric(vcf_info@fix[rind, 2])
# #     genotype_table$REF[rind] <- vcf_info@fix[rind, 4]
# #     genotype_table$ALT[rind] <- vcf_info@fix[rind, 5]
# #     genotype_table$groupID[rind] <- groupID
#     
# #     gt_1 <- sort(as.numeric(unlist(strsplit(vcf_gt[rind, 1], "/"))))
# #     gt_2 <- sort(as.numeric(unlist(strsplit(vcf_gt[rind, 2], "/"))))
# #     
# #     uniq_gt_1 <- unique(gt_1)
# #     uniq_gt_2 <- unique(gt_2)
# #     
# #     if(length(uniq_gt_1) == 1){ # if the first one is homozygous
# #       
# #       if(length(uniq_gt_2) == 1){
# #         # both are homozygous
# #         if(uniq_gt_1 == uniq_gt_2){
# #           # both alleles are the same
# #           genotype_table$GT[rind] <- 0
# #         }else genotype_table$GT[rind] <- 2 # both alleles are different
# #         
# #       }else{
# #         # if the second one is heterozygous
# #         if(length(unique(c(uniq_gt_1, uniq_gt_2))) == 2 ){
# #           # one allele is the same
# #           genotype_table$GT[rind] <- 1
# #         }else genotype_table$GT[rind] <- 2 # both alleles are different
# #   
# #       }
# #     }else{ # if the first one is heterozygous #####
# #       
# #       if(length(uniq_gt_2) == 1){
# #         # if the second one is homozygous
# #         if(length(unique(c(uniq_gt_1, uniq_gt_2))) == 2 ){
# #           # one allele is the same
# #           genotype_table$GT[rind] <- 1
# #         }else genotype_table$GT[rind] <- 2 # both alleles are different
# #         
# #       }else{
# #         # if both are heterozygous
# #         # shared_allele_num <- length(intersect(uniq_gt_1, uniq_gt_2)) 
# #         switch(as.character(length(intersect(uniq_gt_1, uniq_gt_2)) ),
# #                "0" = genotype_table$GT[rind] <- 2, # both alleles are different
# #                "1" = genotype_table$GT[rind] <- 1, # one allele is different
# #                "2" = genotype_table$GT[rind] <- 0) # both alleles are the same 
# #       }
# #     }
#     # proc.time() - ptm
#     
#   }# inner for loop
  cat("ID Check Done!\n")
  save(genotype_table, file = paste0(output_dir, groupID, ".RData"))
  cat(groupID, ".RData is saved under ", output_dir, "\n")
}