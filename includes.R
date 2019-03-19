library(ggbiplot)
library(arules)
library('Hmisc')
library("ggpubr")
library(sva)
library(limma)
library(plyr)
library("caret")
library('doParallel')
library('data.table')
cl <- makeCluster(4)
registerDoParallel(cl)
require("VennDiagram")
library(devtools)
require(ggplot2)
library(xlsx)

#For classifier
#install_github("GabrielHoffman/variancePartition")
library(variancePartition)
library('variancePartition')
library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)
#For plotting and for clustering
library(dynamicTreeCut)
library(gridExtra)
library(cowplot)
#in Gene Enrichment analysis
library(clusterProfiler)
library(DOSE)
library(anchors)
library(plyr)
library(stringr)
library(devtools)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)


#BiocManager::install("clusterProfiler", version = "3.8")
#install_github("mategarb/R.ROSETTA")

library(R.ROSETTA)

#For clustering and then indicating the different clusters on PCA used in Normalization script
library(RColorBrewer)
library(scales)
library(ggfortify)


