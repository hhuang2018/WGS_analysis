library(vcfR)

vcf_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"
#vcf_file_dir <- "../paired_vcf/"
all_vcf_files <- list.files(vcf_file_dir, pattern = "\\.vcf.gz$")
output_dir <- "../Output/SKAT_format"

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
  for(rind in 1:num_rows){
    
    genotype_table$CHROM[rind] <- vcf_info@fix[rind, 1]
    genotype_table$POS[rind] <- as.numeric(vcf_info@fix[rind, 2])
    genotype_table$REF[rind] <- vcf_info@fix[rind, 4]
    genotype_table$ALT[rind] <- vcf_info@fix[rind, 5]
    genotype_table$groupID[rind] <- groupID
    
    gt_1 <- sort(as.numeric(unlist(strsplit(vcf_gt[rind, 1], "/"))))
    gt_2 <- sort(as.numeric(unlist(strsplit(vcf_gt[rind, 2], "/"))))
    
    uniq_gt_1 <- unique(gt_1)
    uniq_gt_2 <- unique(gt_2)
    
    if(length(uniq_gt_1) == 1){ # if the first one is homozygous
      
      if(length(uniq_gt_2) == 1){
        # both are homozygous
        if(uniq_gt_1 == uniq_gt_2){
          # both alleles are the same
          genotype_table$GT[rind] <- 0
        }else genotype_table$GT[rind] <- 2 # both alleles are different
        
      }else{
        # if the second one is heterozygous
        if(length(unique(c(uniq_gt_1, uniq_gt_2))) == 2 ){
          # one allele is the same
          genotype_table$GT[rind] <- 1
        }else genotype_table$GT[rind] <- 2 # both alleles are different
  
      }
    }else{ # if the first one is heterozygous #####
      
      if(length(uniq_gt_2) == 1){
        # if the second one is homozygous
        if(length(unique(c(uniq_gt_1, uniq_gt_2))) == 2 ){
          # one allele is the same
          genotype_table$GT[rind] <- 1
        }else genotype_table$GT[rind] <- 2 # both alleles are different
        
      }else{
        # if both are heterozygous
        # shared_allele_num <- length(intersect(uniq_gt_1, uniq_gt_2)) 
        switch(as.character(length(intersect(uniq_gt_1, uniq_gt_2)) ),
               "0" = genotype_table$GT[rind] <- 2, # both alleles are different
               "1" = genotype_table$GT[rind] <- 1, # one allele is different
               "2" = genotype_table$GT[rind] <- 0) # both alleles are the same 
      }
    }
  }# inner for loop
  
  save(genotype_table, file = paste0(output_dir, groupID, ".RData"))
}