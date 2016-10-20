#!/usr/bin/env Rscript

####################### Parse INPUT options #######################
if (!require("optparse", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("optparse", dependencies = TRUE)
  library("optparse", verbose=F, warn.conflicts =F)
}

if (!require("vcfR", quietly=TRUE, warn.conflicts = FALSE)) {
  install.packages("vcfR", dependencies = TRUE)
  library("vcfR", verbose=F, warn.conflicts =F)
}


# Source files
# file.sources = list.files("lib", pattern="*.R$",
#                           #file.sources = list.files(paste(Sys.getenv('MWAS_DIR'),'/lib',sep=''), pattern="*.R$",
#                           full.names=TRUE, ignore.case=TRUE)
# invisible(sapply(file.sources, source, .GlobalEnv))
source('util.R', echo = FALSE)


option_list = list(
  make_option(c("-i", "--input_fp"), type="character", default=NULL, 
              help="input VCF file path [required]", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, 
              help="output directory [required]", metavar="character"),
  make_option(c("-s", "--summarize", type="logical", action="store_true", default=FALSE,
                help="Summarize the input VCF file", metavar="logical"))
)

opts <- parse_args(OptionParser(option_list=option_list),
                   args=commandArgs(trailing=TRUE))

# create output directory if needed
if(opts$outdir != ".") dir.create(opts$outdir, showWarnings=FALSE, recursive=TRUE)

if(opts$summarize){ # if needs to summarize

    vcf_info <- read.vcfR(opts$input_fp, verbose = FALSE)
    
    if(is.null(colnames(vcf_info@gt))){
      if(dim(vcf_info@gt)[2] == 2){
        colnames(vcf_info@gt) <- c("FORMAT", "SAMPLE")
      }else{
        colnames(vcf_info@gt) <- c("FORMAT", sapply(1:(dim(vcf_info@gt)[2]-1), function(x) paste0("SAMPLE", x)))
      }
    }
    
    vcf_gt <- extract.gt(vcf_info, element = "GT")

    
  
}