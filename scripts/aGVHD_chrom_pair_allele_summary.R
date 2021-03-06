library(vcfR)
source('util.R', echo = FALSE)

vcf_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/paired_vcf/"
#vcf_file_dir <- "../paired_vcf/"
all_vcf_files <- list.files(vcf_file_dir, pattern = "\\.vcf.gz$")
all_vcf_files <- all_vcf_files[grepl("a_", all_vcf_files)]

output_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/Allele_summary_by_chrome/"

existed_file <- list.files(output_dir, pattern = "\\_chr1.RData$")
existed_file <- gsub("_chr1.RData", ".vcf.gz", existed_file)
existed_file <- existed_file[grepl("a_", existed_file)]

all_vcf_files <- all_vcf_files[!(all_vcf_files %in% existed_file)]

num_files <- length(all_vcf_files)

# SKAT genotype code: 
#    0 - both alelles are the same 
#    1 - only one allele is the same
#    2 - both alleles are different


for(id in 1:num_files){
  
  vcf_file <- paste0(vcf_file_dir, all_vcf_files[id])
  vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  
  for(chr in 1:22){
    
    ptm <- proc.time() 
    
    chrom <- paste0("chr", chr)
    
    chr_index <- which(vcf_info@fix[, 1] == chrom)
    vcf_gt <- extract.gt(vcf_info, element = "GT")
    vcf_gt <- vcf_gt[chr_index, ]
    
    groupID  <- gsub(".vcf.gz", "", all_vcf_files[id])
    
    num_rows <- dim(vcf_gt)[1]
    genotype_table <- data.frame(CHROM = character(num_rows),
                                 POS = numeric(num_rows),
                                 REF = character(num_rows),
                                 ALT = character(num_rows),
                                 GT = numeric(num_rows),
                                 groupID = character(num_rows),
                                 stringsAsFactors = FALSE)
    genotype_table$CHROM <- rep(chrom, times = num_rows)
    genotype_table$POS <- as.numeric(vcf_info@fix[chr_index, 2])
    genotype_table$REF <- vcf_info@fix[chr_index, 4]
    genotype_table$ALT <- vcf_info@fix[chr_index, 5]
    genotype_table$groupID <- rep(groupID, times = num_rows)
    
    #rm(vcf_info)
    
    cat("File ", id, ": ", all_vcf_files[id], " Chromosome ", chr, "\n")
    # cat("[")
    
    is_same_gt <- sapply(1:num_rows, function(x) vcf_gt[x, 1] == vcf_gt[x, 2])
    
    diff_gt_index <- which(!is_same_gt)
    genotype_table$GT[diff_gt_index] <- sapply(1:length(diff_gt_index), function(x) is.same.gt(vcf_gt[diff_gt_index[x], ]))
    
    # aa <- proc.time() - ptm
    print(proc.time() - ptm)
    
    cat("ID Check Done!\n")
    save(genotype_table, file = paste0(output_dir, groupID, "_", chrom,".RData"))
    cat(groupID, "_", chrom, ".RData is saved under ", output_dir, "\n", sep = "")
    rm(vcf_gt)
  }
  rm(vcf_info)
}