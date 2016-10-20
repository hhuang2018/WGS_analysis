# list .vcf.gz files
original_file_dir = "/mnt/cloudbiodata_nfs_1/hli_scratch/hhuang/hli_vcf_annotated/"
destination_dir = "/mnt/cloudbiodata_nfs_2/users/hhuang/vcf_missense_variants/"
snpEff_dir = "/home/hhuang/tools/snpEff/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
filenames <- gsub(".vcf.gz", "", filenames)

annotated_filenames <- list.files(destination_dir, pattern = "\\.vcf.gz$")
annotated_filenames <- annotated_filenames[grepl("a_", annotated_filenames)]
annotated_filenames <- gsub("_annotated.vcf.gz", "", annotated_filenames)

Not_annotated_index <- which(!(filenames %in% annotated_filenames))

num_files <- length(Not_annotated_index)

##### need csv format of summary
for(id in 1:num_files){
  system(paste0("java -Xmx16g -jar ", snpEff_dir, "SnpSift.jar filter \"ANN[*].EFFECT has 'missense_variant'\" ",
                original_file_dir, filenames[Not_annotated_index[id]], ".vcf.gz | bgzip > ", 
                destination_dir, filenames[Not_annotated_index[id]],"_missense_any.vcf.gz"), 
         intern = TRUE)
}

# java -jar /home/hhuang/tools/snpEff/SnpSift.jar filter "ANN[*].EFFECT has 'missense_variant'" hli_vcf_annotated/a_102_D_41790023_annotated.vcf.gz | bgzip >  vcf_missense_variables/a_102_D_41790023_missense_any.vcf.gz