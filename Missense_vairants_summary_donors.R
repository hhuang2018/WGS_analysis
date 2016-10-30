source('util.R', echo = FALSE)

# load("../Data/GRCh38_gene_list.RData")

library(vcfR)

VCF_file_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants/"
output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants/stats/"

all_vcf_files <- list.files(VCF_file_dir, pattern = "\\.vcf.gz$")

file_names <- all_vcf_files[grepl("_D_", all_vcf_files)] # 216 donors; 240 recipients

num_files <- length(file_names)  
missense_stats <- data.frame(CHROM=character(0), POS=numeric(0), REF = character(0),
                             ALT.A=numeric(0), ALT.T=numeric(0), ALT.G = numeric(0), ALT.C = numeric(0),
                             ALT.insertion = numeric(0), ALT.deletion = numeric(0), 
                             NumDiff = numeric(0), stringsAsFactors = F)
for(id in 1:num_files){
  
  vcf_file <- paste0(VCF_file_dir, file_names[id])
  vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  temp_info <- vcf_info@fix
  vcf_missense <- temp_info[temp_info[,"FILTER"] == "PASS", c("CHROM", "POS", "REF","ALT")]
  
  num_variants <- dim(vcf_missense)[1]
  
  if(dim(missense_stats)[1] > 0){ # if not the first sample
    vcf_chrom_pos <- sapply(1:dim(vcf_missense)[1], function(x) paste0(vcf_missense[x, c("CHROM","POS")], collapse = ""))
    all_chrom_pos <- sapply(1:dim(missense_stats)[1], function(x) paste0(missense_stats[x, c("CHROM","POS")], collapse = ""))
    
    inter_chrom_pos <- intersect(vcf_chrom_pos, all_chrom_pos) # intersection of vcf and all_stats
    
    stats_exist_chrom_pos_id <- which(all_chrom_pos %in% inter_chrom_pos)
    new_chrom_pos_id <- which(!vcf_chrom_pos %in% inter_chrom_pos)
    vcf_exist_chrom_pos_id <- which(vcf_chrom_pos %in% inter_chrom_pos)
    
    if(length(stats_exist_chrom_pos_id) >0){ # existing POS
      missense_stats$NumDiff[stats_exist_chrom_pos_id] <- missense_stats$NumDiff[stats_exist_chrom_pos_id] + 1
      
      missense_mutations <- count_missense_mutation_type(vcf_missense[vcf_exist_chrom_pos_id,])
      
      same_order_id <- which(vcf_chrom_pos[vcf_exist_chrom_pos_id] == all_chrom_pos[stats_exist_chrom_pos_id]) ### 
      if(length(same_order_id) >0){ # same ordered ID (POS order)
        
        missense_stats[stats_exist_chrom_pos_id[same_order_id], c("ALT.A", "ALT.T", "ALT.G", "ALT.C")] <- 
          missense_stats[stats_exist_chrom_pos_id[same_order_id], c("ALT.A", "ALT.T", "ALT.G", "ALT.C")] + 
          missense_mutations[same_order_id, c("ALT.A", "ALT.T", "ALT.G", "ALT.C")]
        
      }
      if(length(same_order_id) < length(vcf_exist_chrom_pos_id)){ # if there's different ordered ID
        
        ################# debug this part ###################
        if(length(same_order_id) >0){ # if there's at least one same order ID
          
          vcf_diff_order <- vcf_exist_chrom_pos_id[-same_order_id]
          all_diff_order <- sapply(vcf_diff_order, function(x) which(all_chrom_pos[stats_exist_chrom_pos_id] == vcf_chrom_pos[x]))
          
          missense_stats[stats_exist_chrom_pos_id[all_diff_order], c("ALT.A", "ALT.T", "ALT.G", "ALT.C")] <- 
            missense_stats[stats_exist_chrom_pos_id[all_diff_order], c("ALT.A", "ALT.T", "ALT.G", "ALT.C")] + 
            missense_mutations[vcf_diff_order, c("ALT.A", "ALT.T", "ALT.G", "ALT.C")]
          
        }else { # If all orders are not the same
          
          vcf_diff_order <- vcf_exist_chrom_pos_id
          all_diff_order <- sapply(vcf_diff_order, function(x) which(all_chrom_pos[stats_exist_chrom_pos_id] == vcf_chrom_pos[x]))
          
          missense_stats[stats_exist_chrom_pos_id[all_diff_order], c("ALT.A", "ALT.T", "ALT.G", "ALT.C")] <- 
            missense_stats[stats_exist_chrom_pos_id[all_diff_order], c("ALT.A", "ALT.T", "ALT.G", "ALT.C")] + 
            missense_mutations[, c("ALT.A", "ALT.T", "ALT.G", "ALT.C")]
          
        }
        
      } # if there's different ordered ID
      
    } 
    
    if(length(new_chrom_pos_id) > 0){ # new POS
      
      missense_mutations <- count_missense_mutation_type(vcf_missense[new_chrom_pos_id,])
      
      temp_stats <- data.frame(vcf_missense[new_chrom_pos_id, c("CHROM", "POS")], missense_mutations, NumDiff = numeric(length(new_chrom_pos_id))+1, 
                               stringsAsFactors = F)
      
      temp_stats$CHROM <- as.character(temp_stats$CHROM)
      temp_stats$POS <- as.numeric(as.character(temp_stats$POS))
      
      missense_stats <- rbind(missense_stats, temp_stats)
      
    }
    
  } else { # the first sample
    
    
    missense_mutations <- count_missense_mutation_type(vcf_missense)
    # het_index <- as.numeric(which(sapply(1:num_variants, function(x) nchar(vcf_missense[x, "REF"]) != nchar(vcf_missense[x, "ALT"]))))
    # het_num <- length(het_index)
    # hom_index <- (1:num_variants)[-het_index]
    # hom_num <- length(hom_index)
    # missense_mutatations <- data.frame(REF = character(num_variants), ALT.A = numeric(num_variants), 
    #                                   ALT.T=numeric(num_variants), ALT.G=numeric(num_variants), 
    #                                   ALT.C=numeric(num_variants), stringsAsFactors = F)
    # # A - 1; T - 2; G - 3; C - 4;
    # 
    # # homozygous variants
    # missense_mutatations$REF[hom_index] <- vcf_missense[hom_index, "REF"]
    # for(hom_id in hom_index){
    #   missense_mutatations[hom_id, which(c("A", "T", "G", "C") %in% vcf_missense[hom_id, "ALT"])+1] <- 
    #     missense_mutatations[hom_id, which(c("A", "T", "G", "C") %in% vcf_missense[hom_id, "ALT"])+1] + 1
    # }
    # 
    # # heterozygous variants
    # missense_mutatations$REF[het_index] <- vcf_missense[het_index, "REF"]
    # for(het_id in het_index){
    #   missense_mutatations[het_id, which(c("A", "T", "G", "C") %in% unlist(strsplit(vcf_missense[het_id, "ALT"], ",")))+1] <- 
    #     missense_mutatations[het_id, which(c("A", "T", "G", "C") %in% unlist(strsplit(vcf_missense[het_id, "ALT"], ",")))+1] + 1
    # }
    
    
    missense_stats <- rbind(missense_stats, cbind(vcf_missense[, c("CHROM", "POS")], missense_mutations, NumDiff = numeric(num_variants)+1))
    
    missense_stats$CHROM <- as.character(missense_stats$CHROM)
    missense_stats$POS <- as.numeric(as.character(missense_stats$POS))
    # missense_stats[, c("ALT.A", "ALT.T", "ALT.G", "ALT.C", "NumDiff")] <- as.numeric(as.character(missense_stats[, c("ALT.A", "ALT.T", "ALT.G", "ALT.C", "NumDiff")]))

  }
}

donor_missense_stats <- missense_stats
save(donor_missense_stats, file = paste0(output_dir, "donor_missesense_stats_updated.RData"))

for(id in 1:num_files){
  
  vcf_file <- paste0(VCF_file_dir, file_names[id])
  vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  temp_info <- vcf_info@fix
  vcf_missense <- temp_info[temp_info[,"FILTER"] == "PASS", c("CHROM", "POS", "REF","ALT")]
  
  index <- which(missense_stats[, "CHROM"] == "chr4" & missense_stats[, "POS"] == 68647129)
  index2 <- which(vcf_missense[, "CHROM"] == "chr4" & vcf_missense[, "POS"] == 68647129)
  # index <- which(missense_mutations[, "CHROM"] == "chr4" & missense_stats[, "POS"] == 68647129)
  if(length(index2) > 0){
    cat(file_names[id], "\n")
    print(missense_stats[index, ])
    print(vcf_missense[index2, ])
    cat("\n")
  }
  
}