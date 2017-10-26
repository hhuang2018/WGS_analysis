source('util.R', echo = FALSE)

# load("../Data/GRCh38_gene_list.RData")

# library(vcfR)

VCF_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/wwang/MiHAIP/DtoR/"
output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants_wwVersion/"

all_vcf_files <- list.files(VCF_file_dir, pattern = "\\.txt$")

file_names <- all_vcf_files[grepl("_R_", all_vcf_files)] # 216 donors; 240 recipients

# AML patients only
load("../Data/HLI_available_pairs_dis_table.RData")

AML_ids <- 1:length(Available_paired_table$Disease)#which(Available_paired_table$Disease %in% "AML") # 196
RID <- sapply(1:length(file_names), function(x) as.numeric(unlist(strsplit(file_names[x], "_"))[4]))

AML_RID_ids <- which(RID %in% Available_paired_table$RID[AML_ids])
file_names <- file_names[AML_RID_ids]
#####
num_files <- length(file_names)  
missense_mm_stats <- vector(mode = "integer", length = num_files)
for(id in 1:num_files){
  
  vcf_file <- paste0(VCF_file_dir, file_names[id])
  # vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  # temp_info <- vcf_info@fix
  # vcf_missense <- temp_info[temp_info[,"FILTER"] == "PASS", c("CHROM", "POS", "REF","ALT")]
  # 
  vcf_missense <- read.table(vcf_file)
  
  missense_mm_stats[id] <- dim(vcf_missense)[1]
  
   
}

save(missense_mm_stats, file = paste0(output_dir, "All_missense_mismatch_stats_wwVersion.RData"))
