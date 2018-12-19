library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)


#Protein coding genes and those that are unannotated (To run feature selection only on protein coding genes)
filter=getAnnotatedGenes(genes,colnames(TARGET_GE_Classifier))
protein_coding=filter[[1]][which(filter[[1]]$gene_type =="protein_coding"),geneID]

#Try Feature selection with all RNAs or genes except psudo genes and IGV


#Have to run TARGETAMLPreprocess and TARGETAMLfunctions first 
TARGET_GE_Classifier_ProteinCoding=cbind.data.frame(TARGET_GE_Classifier[,append(protein_coding,filter[[2]])],decision)
#All Genes
TARGET_GE_Classifier=cbind.data.frame(TARGET_GE_Classifier,decision)
#MCFS and Boruta are run on ULAM
write.table(TARGET_GE_Classifier,"MCFS_AML_GE_Classifier",row.names=TRUE)
#After filtering out according to CPM normalization
#matched 3  folder is output of MCFS for all genes after filtering according edgeR
#matched 4 folder is output of MCFS for protein coding and unannotated genes after filtering according edgeR
MCFSFeatures=FilterFeatures("MCFSoutput/matched 4/matched__RI.csv",200)
MCFSFeatures=FilterFeatures("MCFSoutput/matched 3/matched__RI.csv",200)


#---Before applying the CPM filter using edge R method
boruta.train=readRDS("TARGET_AML_GE_BORUTA.rds")
#After filtering out according to CPM using edge R method in normalization
boruta.train=readRDS("TARGET_AML_GE_BORUTA1.rds")
BorutaFeatures=extractFeaturesBoruta(boruta.train)
  
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
ensmblid <- getBM(attributes=c('hgnc_symbol','chromosome_name','start_position','end_position'),
  filters=('ensembl_gene_id'), values =MCFSFeatures, mart = ensembl)

ensmblid <- getBM(attributes=c('hgnc_symbol','chromosome_name','start_position','end_position'),
                  filters=('ensembl_gene_id'), values = MCFSFeatures[which(!(MCFSFeatures %in% genes$geneID  ))], mart = ensembl)
filteredAll=getAnnotatedGenes(genes,MCFSFeatures)

genes <- readGeneDB("gencode.v20.annotation.gtf")
genes$geneID=gsub("(*)\\.([0-9]+)","\\1",genes$gene_id)
#20 Features when using all genes 
#30 Features when running protein coding genes on MCFS
filteredMCFS=getAnnotatedGenes(genes,MCFSFeatures[1:20])
filteredMCFS=getAnnotatedGenes(genes,MCFSFeatures[1:30])
#MCFS genes after removing all pseudo genes and comparing accuracies than getting the optimal number of features
filteredMCFSNoPseudo=getAnnotatedGenes(genes,classifierMCFSgenesNoPseudo[1:17])
write.table(filteredMCFS[[1]][,-c("attributes")],"filteredGenesMCFS",row.names=TRUE)
filteredBoruta=getAnnotatedGenes(genes,BorutaFeatures)
write.table(filteredBoruta[[1]][,-c("attributes")],"filteredGenesBoruta",row.names=TRUE)

plotNames=gsub("(TARGET-20-((.)+-(R|T)))","\\2",rownames(TARGET_GE_Classifier))

plotDistributionCategory(list(filteredMCFS[[1]]$gene_type,filteredBoruta[[1]]$gene_type),c("MCFS","Boruta"),"Distribution_Classes_FS",getwd(),"")

#Only MCFS after removing the pseudo genes from the first 20 genes
plotDistributionCategory(list(filteredMCFS[[1]]$gene_type),c("MCFS"),"Distribution_Classes_FS_MCFS_WithoutPseudoGenes","ResultsRules/AllGenes","")

#Plot distribution of expirement 2 with removing psuedo genes from 200 MCFS list then getting 17 as the optimal number of features
plotDistributionCategory(list(filteredMCFSNoPseudo[[1]]$gene_type),c("MCFS"),"Distribution_Classes_FS","ResultsRules/AllGenes","")


#Without PseudoGenes MCFS
classifierMCFSgenes=append(filteredMCFS[[1]][which(!(grepl("*\\_pseudogene",filteredMCFS[[1]]$gene_type))),]$geneID,filteredMCFS[[2]])
#Without PseudoGenes MCFS
classifierMCFSgenesNoPseudo=append(filteredAll[[1]][which(!(grepl("*\\_pseudogene",filteredAll[[1]]$gene_type))),]$geneID,filteredAll[[2]])

#All Genes MCFS
classifierMCFSgenes=append(filteredMCFS[[1]]$geneID,filteredMCFS[[2]])
#Ordering according to MCFS Features
classifierMCFSgenes=classifierMCFSgenes[!is.na(match(MCFSFeatures, classifierMCFSgenes))]

annotation=annotateInOrder(filteredMCFS,classifierMCFSgenes)

#classifierMCFSgenes=MCFSFeatures[1:30]
#Without PseudoGenes
classifierBorutagenes=append(filteredBoruta[[1]][which(!(grepl("*\\_pseudogene",filteredBoruta[[1]]$gene_type))),]$geneID,filteredBoruta[[2]])
#All Genes
classifierBorutagenes=append(filteredBoruta[[1]]$geneID,filteredBoruta[[2]])

my.plots <- vector(2, mode='list')

plotPCA(getwd(),TARGET_GE_Classifier[,classifierMCFSgenes],decision,"PCA_MCFS")
my.plots[[1]]=recordPlot()
plotPCA(getwd(),TARGET_GE_Classifier[,classifierBorutagenes],decision,"PCA_Boruta")
my.plots[[2]]=recordPlot()
#savePDF(my.plots,"PCA_AllGenes",getwd())
savePDF(my.plots,"PCA_WithoutPseudogenes",getwd())

TARGET_GE_Classifier_MCFS=TARGET_GE_Classifier[,append(classifierMCFSgenes,"decision")]
#colnames(TARGET_GE_Classifier_MCFS)=append(filteredMCFS[[1]]$gene_name[which(filteredMCFS[[1]]$geneID %in% classifierMCFSgenes)],append(filteredMCFS[[2]],"decision"))
colnames(TARGET_GE_Classifier_MCFS)=append(annotation,"decision") 

TARGET_GE_Classifier_Boruta=TARGET_GE_Classifier[,append(classifierBorutagenes,"decision")]
colnames(TARGET_GE_Classifier_Boruta)=append(filteredBoruta[[1]]$gene_name[which(filteredBoruta[[1]]$geneID %in% classifierBorutagenes)],append(filteredBoruta[[2]],"decision"))
#Compare accuracies on all MCFS features even with the pseudoGenes included
Accuracies=compareAccuracies(TARGET_GE_Classifier,200,MCFSFeatures)
#Compare accuracies on all MCFS features but without pseudoGenes included
Accuracies=compareAccuracies(TARGET_GE_Classifier,length(classifierMCFSgenesNoPseudo),classifierMCFSgenesNoPseudo)
#Compare accuracies on all MCFS features when running MCFS on only protein coding genes 
Accuracies=compareAccuracies(TARGET_GE_Classifier,length(classifierMCFSgenesNoPseudo),classifierMCFSgenesNoPseudo)


resultRosettaMCFS=rosetta(TARGET_GE_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3)

recalculatedResultRosettaMCFS=recalculateRules(TARGET_GE_Classifier_MCFS,resultRosettaMCFS$main)

write.csv(recalculatedResultRosettaMCFS,paste("ResultsRules/AllGenes/RulesAllGenes-",substr(date(),1,10),".csv",sep=""))
saveLineByLine(recalculatedResultRosettaMCFS,  paste("ResultsRules/AllGenes/NetworksAllGenes-",substr(date(),1,10),".txt",sep=""))

write.csv(recalculatedResultRosettaMCFS,paste("ResultsRules/ProteinCoding/RulesAllGenes-",substr(date(),1,10),".csv",sep=""))
saveLineByLine(recalculatedResultRosettaMCFSGenetic, paste("ResultsRules/AllGenes/NetworksAllGenes-",substr(date(),1,10),".txt",sep=""))

resultRosettaMCFSGenetic=rosetta(TARGET_GE_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3,reducer="Genetic",ruleFiltration=TRUE, ruleFiltrAccuracy=c(0,0.8),ruleFiltrSupport=c(1,3))
#resultRosettaMCFSGenetic=rosetta(TARGET_GE_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3,reducer="Genetic")
recalculatedResultRosettaMCFSGenetic=recalculateRules(TARGET_GE_Classifier_MCFS,resultRosettaMCFSGenetic$main)
write.csv(recalculatedResultRosettaMCFS,"ResultsRules/AllGenes/RulesGenetic.csv")
saveLineByLine(recalculatedResultRosettaMCFSGenetic, "ResultsRules/AllGenes/NetworksAllGenesGenetic.txt")

resultRosettaBoruta=rosetta(TARGET_GE_Classifier_Boruta,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency", discreteParam=3)
clusteredRulesMCFS=clusterRules(recalculatedResultRosettaMCFS,rownames(TARGET_GE_Classifier_MCFS))
colnames(clusteredRulesMCFS)<-plotNames
my.plots=vector(1, mode='list');
svg(paste("ResultsRules/AllGenes/HeatMapAllGenes-",substr(date(),1,10),".svg",sep=""))
clusters=heatmap.F(t(clusteredRulesMCFS), colors=c('white','white','white','white','blue','blue'),distmethod='pearson')
dev.off()
write.csv(clusters,paste("ResultsRules/AllGenes/Clusters-",substr(date(),1,10),".csv",sep = ""))
my.plots[[1]]=recordPlot()
savePDF(my.plots,paste("HeatMapAllGenes-",substr(date(),1,10),sep=""),"ResultsRules/AllGenes/")

