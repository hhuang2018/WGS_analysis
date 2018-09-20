
count_based_fp = 'Documents/NGSProject/2018WGS/Data/GWASH/unimputed/LD_pruned/Encoding/Count_based/mRMR/count_based_agvhi24_SNPs.csv'

feat_table <- read.csv(count_based_fp, header = T, stringsAsFactors = F)

dist_matrix <- dist(feat_table[, -1], method="manhattan") # 
## PCA
library(FactoMineR)

# apply PCA
pca3 = PCA(USArrests, graph = FALSE)

# matrix with eigenvalues
pca3$eig

## MDS
fit <- cmdscale(d,eig=TRUE, k=2) 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  main="Metric	MDS",	type="n")
text(x, y, labels = row.names(mydata), cex=.7)
# Nonmetric MDS
# N rows (objects) x p columns (variables)
# each row identified by a unique row name

library(MASS)
d <- dist(mydata) # euclidean distances between the rows
fit <- isoMDS(d, k=2) # k is the number of dim
fit # view results

# plot solution 
x <- fit$points[,1]
y <- fit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
  main="Nonmetric MDS", type="n")
text(x, y, labels = row.names(mydata), cex=.7)