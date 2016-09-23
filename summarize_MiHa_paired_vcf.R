source('util.R', echo = FALSE)

load("../Data/ID_table.RData")

if (!require("vcfR", quietly=TRUE, warn.conflicts = FALSE)){
  install.packages("vcfR", dependencies = TRUE)
  library("vcfR", verbose=F, warn.conflicts =F)
}

paired_DtoR_fp <- "/mnt/cloudbiodata_nfs_1/hli_scratch/wwang/MiHAIP/DtoR/"
paired_files <- list.files(path = paired_DtoR_fp, pattern = "\\.txt$")
output_fp <- "/mnt/cloudbiodata_nfs_2/users/hhuang/MiHA_missense_Summary/"

num_files <- length(paired_files)

aGVHD_SNP_list <- list()
nGVHD_SNP_list <- list()

for(id in 1:num_files){
  
  paired_info <- read.table(paste0(paired_DtoR_fp, paired_files[id]))
  
  groupID <- as.integer(gsub(".txt", "", paired_files[id]))

  colnames(paired_info) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
  paired_info$POS <- as.integer(paired_info$POS)
  
  CHROMs <- unique(paired_info$CHROM)
  num_CHROMs <- length(CHROMs)
  
  summary_table <- data.frame(CHROM = character(num_CHROMs),
                              NumVariants = numeric(num_CHROMs),
                              # NumGenes = character(num_CHROMs),
                              stringsAsFactors = FALSE)

  
  for(chr in 1:num_CHROMs){
    # ptm <- proc.time() 
    
    chrom <- CHROMs[chr]
    # chrom_num <- as.numeric(gsub("ch))
    chr_index <- which(paired_info$CHROM == chrom)
    
    num_rows <- length(chr_index)
    
    summary_table$CHROM[chr] <- chrom
    summary_table$NumVariants[chr] <- num_rows
    
    #######
    # Check all the missense variants in the list
    #######
    SNP_summary <- data.frame(CHROM = vector(mode = "character", length = num_rows),
                              POS = vector(mode = "integer", length = num_rows), 
                              NumDiff = vector(mode = "integer", length = num_rows),
                              # total_num = vector(mode = "integer", length = num_rows),
                              stringsAsFactors = F)
    counter <- 0
    
    SNP_summary$CHROM <- rep(chrom, length = num_rows)
    SNP_summary$POS <- paired_info$POS[chr_index]
    SNP_summary$NumDiff <- SNP_summary$NumDiff + 1
    
    # aGvHD group 
    if(ID_table$Group[which(ID_table$GroupID %in% groupID)[1]] == "a"){
      
      # SNP_summary$POS 
      # if it's the first chromosome in the list; or the chromosome table appears for the first time
      if(length(aGVHD_SNP_list) == 0 | !is.element(chrom, names(aGVHD_SNP_list))){
        
        eval(parse(text=paste0("aGVHD_SNP_list <- append(aGVHD_SNP_list, list(", chrom, " = SNP_summary))")))
        
      }else {
        
        # else make union SNP list, and the intersected SNPs +1 for aGVHD, -1 for non-GVHD
        
        original_list <- aGVHD_SNP_list[[chrom]]
        
        samePOS_index_original <- which(original_list$POS %in% intersect(original_list$POS, SNP_summary$POS[1:100]))
        samePOS_index_new <- which(SNP_summary$POS %in% intersect(original_list$POS, SNP_summary$POS[1:100]))
        
        original_list$NumDiff[samePOS_index_original] <- original_list$NumDiff[samePOS_index_original] + 1
        
        aGVHD_SNP_list[[chrom]] <- rbind(original_list, SNP_summary[-samePOS_index_new, ])

      }
        
    }else { # non-GvHD group -1
      
      # SNP_summary$NumDiff <- -SNP_summary$NumDiff
      # SNP_summary$POS 
      # if it's the first chromosome in the list; or the chromosome table appears for the first time
      if(length(nGVHD_SNP_list) == 0 | !is.element(chrom, names(nGVHD_SNP_list))){ # aGVHD +1
        
        eval(parse(text=paste0("nGVHD_SNP_list <- append(nGVHD_SNP_list, list(", chrom, " = SNP_summary))")))
        
      }else {
        
        # else make union SNP list, and the intersected SNPs +1 for aGVHD, -1 for non-GVHD
        original_list <- nGVHD_SNP_list[[chrom]]
        
        samePOS_index_original <- which(original_list$POS %in% intersect(original_list$POS, SNP_summary$POS))
        samePOS_index_new <- which(SNP_summary$POS %in% intersect(original_list$POS, SNP_summary$POS))
        
        original_list$NumDiff[samePOS_index_original] <- original_list$NumDiff[samePOS_index_original] + 1
        
        nGVHD_SNP_list[[chrom]] <- rbind(original_list, SNP_summary[-samePOS_index_new, ])

      } # if-else-newList
    } # if-else-GvHD
  } # for-chrom
  save(summary_table, file = paste0(output_fp, "Summary_", groupID, ".RData"))
}

save(aGVHD_SNP_list, nGVHD_SNP_list, file = paste0(output_fp,"missense_summary_DtoR.RData"))
# load("../Output/paired_VCF_summary/test_stats_DtoR.RData")
