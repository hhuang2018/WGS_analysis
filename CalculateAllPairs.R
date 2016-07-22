source('stats.R', echo = FALSE)

VCF_file_dir <- "/mnt/scratch/hhuang/hli_vcf_renamed/"
# Donor_file <- "n_197_D_62284237_annotated.vcf.gz"
# Recipient_file <- "n_197_R_1005928_annotated.vcf.gz"

all_vcf_files <- list.files(VCF_file_dir, pattern = "\\.vcf.gz$")

# num_files <- length(all_vcf_files)
out.dir <- "/mnt/scratch/hhuang/VariantStats/"

fileName_list <- strsplit(all_vcf_files, "_")
groupIDs <- sapply(fileName_list, function(x) paste0(x[1], "_", x[2]))

duplicated_IDs <- which(duplicated(groupIDs))

num_pairs <- length(duplicated_IDs)

for(id in 1:num_pairs){
  pair_ids <- which(groupIDs %in% groupIDs[duplicated_IDs[id]])
  Donor_file <- ifelse(fileName_list[[pair_ids[1]]][3] == "D", all_vcf_files[pair_ids[1]], all_vcf_files[pair_ids[2]])
  Recipient_file <- ifelse(fileName_list[[pair_ids[1]]][3] == "R", all_vcf_files[pair_ids[1]], all_vcf_files[pair_ids[2]])
  
  cat(id, "-th Group #: ", groupIDs[duplicated_IDs[id]], "\n")
  cat("Donor File: ", Donor_file, "\n")
  cat("Recipient File: ", Recipient_file, "\n")
  normalized_VariantStats_by_Chr(VCF_file_dir, Donor_file, Recipient_file, out.dir)
  
}