source('util.R', echo = FALSE)

# load("../Data/GRCh38_gene_list.RData")

# library(vcfR)

VCF_file_dir <- "/mnt/cloudbiodata_nfs_1/hli_scratch/wwang/MiHAIP/DtoR/"
output_dir <- "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants_wwVersion/"

all_vcf_files <- list.files(VCF_file_dir, pattern = "\\.txt$")

#file_names <- all_vcf_files[grepl("_R_", all_vcf_files)] # 216 donors; 240 recipients

for(id in 1:num_files){
  
  vcf_file <- paste0(VCF_file_dir, all_vcf_files[id])
  # vcf_info <- read.vcfR(vcf_file, verbose = FALSE)
  # temp_info <- vcf_info@fix
  # vcf_missense <- temp_info[temp_info[,"FILTER"] == "PASS", c("CHROM", "POS", "REF","ALT")]
  # 
  vcf_missense <- read.table(vcf_file)
  
  missense_mm_stats[id] <- dim(vcf_missense)[1]
  
   
}

save(missense_mm_stats, file = paste0(output_dir, "All_missense_mismatch_stats_wwVersion.RData"))
