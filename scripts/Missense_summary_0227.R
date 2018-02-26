source('util.R', echo = FALSE)

library(vcfR)
library(dplyr)

### Varaint Stats summary

stats_dir <- "../Output/VariantStats/"

all_stats_files <- list.files(stats_dir, pattern = "\\.RData$")

num_files <- length(all_stats_files)
# a_exon_genes_donor <- character(length = 0)
# a_exon_genes_recipient <- character(length = 0)
# n_exon_genes_donor <- character(length = 0)
# n_exon_genes_recipient <- character(length = 0)

disease_type_file <- "../ClinVar/HLI_DiseaseTypes_2.csv"

disease_type <- read.csv(disease_type_file)


VCF_file_dir <- "../Missense_VCF/"
# VCF_file_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants/"
output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants/stats/"

all_vcf_files <- list.files(VCF_file_dir, pattern = "\\.vcf.gz$")

file_names <- all_vcf_files[grepl("_D_", all_vcf_files)] # 216 donors; 240 recipients

extract_INFO <- function(INFO){
  
  Gene_info <- strsplit(strsplit(INFO, ";")[[1]][3], "\\|")[[1]][c(4,5)] # Gene Name; GeneID;
  
  return(Gene_info)
}
  

num_files <- length(file_names)  
missense_gene_stats <- data.frame(GeneName=character(0), 
                                  GeneID=character(0), 
                                  Count = numeric(0),
                                  stringsAsFactors = F)
AML_missense_gene_stats <- data.frame(GeneName=character(0), 
                                      GeneID=character(0), 
                                      Count = numeric(0),
                                      stringsAsFactors = F)
for(id in 1:num_files){
  
  vcf_file <- paste0(VCF_file_dir, file_names[id])
  vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  temp_info <- vcf_info@fix
  vcf_missense <- temp_info[temp_info[,"FILTER"] == "PASS", c("CHROM", "POS", "REF","ALT", "INFO")]
  
  gene_vcf_INFO <- lapply(vcf_missense[, "INFO"], extract_INFO)
  
  gene_vcf_INFO2 <- do.call(rbind, gene_vcf_INFO)
  gene_vcf_INFO2 <- as.data.frame(gene_vcf_INFO2, stringsAsFactors = F)
  colnames(gene_vcf_INFO2) <- c("GeneName", "GeneID")
  
  missense_gene_count1 <- count_(gene_vcf_INFO2, vars = c("GeneName", "GeneID"))
  
  missense_gene_count <- as.data.frame(missense_gene_count1, stringAsFactors = F)
  
  out.file.name <- strsplit(file_names[id], "\\.")[[1]][1]
  save(missense_gene_count, file = paste0(output_dir, out.file.name, "gene_missense_mutation_stats.RData"))
  
  if(id == 1){# the first file
    
    missense_gene_stats <- missense_gene_count
    
  }else{
    
    temp_stats1 <- setdiff(missense_gene_stats[, c(1,2)], missense_gene_count[, c(1,2)])
    temp_stats2 <- setdiff(missense_gene_count[, c(1,2)], missense_gene_stats[, c(1,2)])
    
    # todo: merge to the 
    
  }
  
  #### AML cohort
  Group_ID <- as.numeric(strsplit(file_names[id], "_")[[1]][2])
  index <- which(Group_ID %in% disease_type$GroupID)
  # dis_type <- ifelse(length(index)!=0, disease_type$disease[index], "Other")
  if(ifelse(length(index)!=0, disease_type$disease[index], "Other") == "AML"){ # AML cohort
    
    print("AML")
    
    if(dim(AML_missense_gene_stats)[1] == 0){ # first record 
      
      AML_missense_gene_stats <- missense_gene_count
      
    } else{
      
      
      
    }
     
  }#else{print("Non_AML")}
}

donor_missense_stats <- missense_stats
save(donor_missense_stats, file = paste0(output_dir, "donor_missesense_stats_updated_1030.RData"))