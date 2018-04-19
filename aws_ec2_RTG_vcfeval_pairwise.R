
require(foreach)

# list .vcf.gz files
original_file_dir = "~/data/HLI/"
destination_dir = "~/data/HLI_output/"

filenames <- list.files(original_file_dir, pattern = "\\.vcf.gz$")
groupnames <- sapply(1:length(filenames), function(x) paste0(unlist(strsplit(filenames[x], "_"))[1:2], collapse ="_"))

unique_groups <- unique(groupnames)

num_files <- length(unique_groups)

##### need csv format of summary
foreach(id = 1:num_files) %dopar% {
  
  index <- which(grepl(unique_groups[id], filenames))
  
  if(length(index) == 2){
    
    donor_file <- filenames[index[grepl('_D_', filenames[index])]]
    recipient_file <- filenames[index[grepl('_R_', filenames[index])]]
  }
  
  system(paste0("rtg vcfeval -b ", original_file_dir, "snpEff.jar ann -csvStats ",
                destination_dir, filenames[Not_annotated_index[id]], "_snpEff_summary.csv -v ",
                "GRCh38.82 ", original_file_dir, filenames[Not_annotated_index[id]], ".vcf.gz > ", 
                destination_dir, filenames[Not_annotated_index[id]],"_annotated.vcf"), 
         intern = TRUE)
  system(paste0("rm ", destination_dir, filenames[Not_annotated_index[id]],"_annotated.vcf"))
  system(paste0("rm ", destination_dir, filenames[Not_annotated_index[id]], "_snpEff_summary.genes.txt"))
  # system(paste0("mv snpEff_summary.csv ",  destination_dir, filenames[Not_annotated_index[id]], "_snpEff_summary.csv"))
  #save(t1, file = paste0(destination_dir, filenames[Not_annotated_index[id]],"_annotated.out.RData"))
}