# >> brew install hdf5 --enable-cxx
#library(devtools)
#install_github("mannau/h5")
# Install GridLMM
# install_github('deruncie/GridLMM')
#source("http://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
########

library(h5)

library(rhdf5)

library(data.table)
library(GridLMM)

#vignette(topic = 'Running_GridLMM_GWAS',package='GridLMM')


# read HD5F files
# The function is from https://www.kaggle.com/jeffmoser/read-hdf-function-for-r  
# by Jeff Moserread_hdf function for R

# Adapted from http://pandas.pydata.org/pandas-docs/stable/io.html#io-external-compatibility 
read_hdf <- function(h5path, dataframe_name=NULL) {
  h5File <- H5Fopen(h5path, flags="H5F_ACC_RDONLY")
  listing <- h5ls(h5File)
  
  if (is.null(dataframe_name)) {
    dataframe_name <- listing$name[1]
  }
  
  group_name <- paste0("/", dataframe_name)
  
  # Filter to just the requested dataframe:
  listing <- listing[listing$group == group_name,]
  
  # Find all data nodes, values are stored in *_values and corresponding column
  # titles in *_items
  data_nodes <- grep("_values$", listing$name)
  name_nodes <- grep("_items$", listing$name)
  data_paths = paste(listing$group[data_nodes], listing$name[data_nodes], sep = "/")
  name_paths = paste(listing$group[name_nodes], listing$name[name_nodes], sep = "/")
  columns = list()
  for (idx in seq(data_paths)) {
    # NOTE: matrices returned by h5read have to be transposed to to obtain
    # required Fortran order!
    data <- data.frame(t(h5read(h5File, data_paths[idx])))
    names <- t(h5read(h5File, name_paths[idx]))
    entry <- data.frame(data)
    colnames(entry) <- names
    columns <- append(columns, entry)
  }

  data <- data.frame(columns)
  
  # If "axis0" is specified, we can return the dataframe columns in the original order:
  if ("axis0" %in% listing$name) {
    orig_col_order <- h5read(h5File, paste0(group_name, "/axis0"))
    data <- data[orig_col_order]
  }
  
  H5Fclose(h5File)
  return(data)
}


read_hdf5 <- function(h5File) {
    listing <- rhdf5::h5ls(h5File)
    # Find all data nodes, values are stored in *_values and corresponding
    # column titles in *_items
    data_nodes <- grep("_values", listing$name)
    name_nodes <- grep("_items", listing$name)

    data_paths <- paste(listing$group[data_nodes], listing$name[data_nodes],
                       sep = "/")
    name_paths <- paste(listing$group[name_nodes], listing$name[name_nodes],
                       sep = "/")

    columns <- list()
    for (idx in seq(data_paths)) {
        data <- data.frame(t(rhdf5::h5read(h5File, data_paths[idx])))
        names <- t(rhdf5::h5read(h5File, name_paths[idx]))
        entry <- data.frame(data)
        colnames(entry) <- names
        columns <- append(columns, entry)
    }

    data <- data.frame(columns)

    return(data)
}

train <- read_hdf("../input/train.h5")

head(train)

###
chr_fp = '../../2018WGS/Data/MismatchEncoded/presabs_based/EncodedMatrix/'

chr_Encoded_file = 'chr1_EncodedMatrix_95filtered_presabs_based.h5'

chr_Encoded_mat <- h5file(paste0(chr_fp, chr_Encoded_file), 'r')
list.datasets(chr_Encoded_mat, recursive = T)
vals <- readDataSet(chr_Encoded_mat)

#testvec <- chr_Encoded_mat[1:2, ]
#h5attr(chr_Encoded_mat)

# Load data.
filePath = paste0(chr_fp, chr_Encoded_file)
h5ls(filePath)

dataset = h5read(file=filePath, name="chr_1")

dataset1 <- read_hdf5(filePath)

block0_values = t(dataset$block0_values)
block1_values = t(dataset$block1_values)
colnames(block0_values) = dataset$block0_items
colnames(block1_values) = dataset$block1_items
dataset = cbind(block0_values, block1_values)
