#Adult Cohort processing
#You have to run Normalization adults first

library("rmcfs")
library("Boruta")
library("R.ROSETTA")
library(biomaRt)



#Linda_GE_Classifier2 is ready now to apply feature selection and everything

write.table(Linda_GE_Classifier2,"MCFS_AML_Adults_GE_Classifier",row.names=TRUE)


Linda_Adults_MCFSFeatures=FilterFeatures("Linda_Adults_Cohort_Results/MCFSoutput/LindaunmatchedAdults1__RI.csv",1000)
names(Linda_GE_Classifier2)[names(Linda_GE_Classifier2)=="decision_Linda2"]="decision"
Linda_Adult_Accuracies=compareAccuracies(Linda_GE_Classifier2,200,as.character(Linda_Adult_MCFSFeatures))

genes=readGeneDB("gencode.v20.annotation.gtf")


Linda_filtered_Adults_MCFS=getAnnotatedGenes(genes,Linda_Adult_MCFSFeatures[1:50])
#All Genes MCFS
Linda_Adults_classifierMCFSgenes=append(Linda_filtered_Adults_MCFS[[1]]$geneID,Linda_filtered_Adults_MCFS[[2]])


#Ordering according to MCFS Features
Linda_Adults_classifierMCFSgenes=Linda_Adults_classifierMCFSgenes[!is.na(match(Linda_Adults_MCFSFeatures, Linda_Adults_classifierMCFSgenes))]


Linda_Adults_annotation=annotateInOrder(Linda_filtered_Adults_MCFS,Linda_Adults_classifierMCFSgenes)

Linda_GE_Adults_Classifier_MCFS=Linda_GE_Classifier2[,append(Linda_Adults_classifierMCFSgenes,"decision")]

colnames(Linda_GE_Adults_Classifier_MCFS)<-append(Linda_Adults_annotation,"decision")

Linda_Adults_resultRosettaMCFSGenetic=rosetta(Linda_GE_Adults_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=10,reducer="Genetic",ruleFiltration=TRUE, discreteParam=3)

Linda_Adults_recalculatedResultRosettaMCFSGenetic=recalculateRules(Linda_GE_Adults_Classifier_MCFS,Linda_Adults_resultRosettaMCFSGenetic$main)


Linda_Adults_resultRosettaMCFS=rosetta(Linda_GE_Adults_Classifier_MCFS,classifier="StandardVoter",discrete=FALSE, discreteMethod="EqualFrequency",cvNum=10, discreteParam=3)

Linda_Adults_recalculatedResultRosettaMCFS=recalculateRules(Linda_GE_Adults_Classifier_MCFS,Linda_Adults_resultRosettaMCFS$main)


Linda_Adults_clusteredRulesMCFS=clusterRules(Linda_Adults_recalculatedResultRosettaMCFS,rownames(Linda_GE_Adults_Classifier_MCFS))

Linda_Adults_clusteredRulesMCFSGenetic=clusterRules(Linda_Adults_recalculatedResultRosettaMCFSGenetic,rownames(Linda_GE_Adults_Classifier_MCFS))

writeOutput("Linda_Adults_Cohort_Results/",Sys.Date(),"AllGenes","Genetic",Linda_Adults_clusteredRulesMCFSGenetic,Linda_Adults_recalculatedResultRosettaMCFSGenetic,Linda_Adults_resultRosettaMCFSGenetic$main,enrichmentLinda,"enrichment",FALSE)

writeOutput("Linda_Adults_Cohort_Results/",Sys.Date(),"AllGenes","Johnson",Linda_clusteredRulesMCFS,Linda_recalculatedResultRosettaMCFS,Linda_resultRosettaMCFS$main,enrichmentLinda,"enrichment",FALSE)



writeOutput("Integrated_Cohort_Results/",Sys.Date(),"AllGenes","Genetic",Integrated_Matched_clusteredRulesMCFSGenetic,Integrated_Matched_recalculatedResultRosettaMCFSGenetic,Integrated_Matched_resultRosettaMCFSGenetic$main,
            enrichmentLinda_Adults,"enrichment",FALSE)

writeOutput("Linda_Adults_Cohort_Results/",Sys.Date(),"AllGenes","Johnson",Linda_Adults_clusteredRulesMCFS,Linda_Adults_recalculatedResultRosettaMCFS,Linda_Adults_resultRosettaMCFS$main,
            enrichmentLinda_Adults,"enrichment",FALSE)


plotDistributionCategory(list(Linda_filtered_Adults_MCFS[[1]]$gene_type),list(c("MCFS")),"Distribution_Classes_FS",paste(getwd(),"/Linda_Adults_Cohort_Results",sep=""),"")


#Enrichment Analysis and finding which features from Rosetta is enriched from the significant biological processes

enrichmentLinda_Adults=GOenrichment(Linda_Adult_MCFSFeatures,colnames(Linda_GE_Classifier2)[1:length(Linda_GE_Classifier2)-1],"ENSEMBL")

Linda_Adults_RosettaEnrichment=findReleventBP(enrichmentLinda_Adults,append(Linda_Adult_MCFSFeatures[1:50],"decision"),Linda_GE_Classifier2,1.0)

BP=unique(sapply(as.data.frame(enrichmentLinda)$Description, function(x) x))

plotEnrichment(enrichmentLinda,"Enrichment Pediatric cohort All Genes","/Linda_Cohort_Results/2019-04-01/AllGenes/.csv")
plotGeneRulesEnrichment(Linda_RosettaEnrichment,"Linda_Cohort_Results/","2019-04-02","AllGenes")

#---------------------------Trying to predict the 36 samples with the TARGET_Classifier--------------------------------------

Linda_GE_Classifier_temp=Linda_GE_Classifier[,c(MCFSFeatures[1:20],"decision_Linda")]
filteredMCFS_Linda=getAnnotatedGenes(genes,colnames(Linda_GE_Classifier_temp)[1:20])
annotation_Linda=annotateInOrder(filteredMCFS_Linda,classifierMCFSgenes)
colnames(Linda_GE_Classifier_temp)=append(annotation_Linda,"decision") 
predictClass(Linda_GE_Classifier_temp, resultRosettaMCFS$main)
#------------------------------#What are the common genes between TARGET and Linda_Classifier---------------------------------------------------------------------------
geneLS=list(MCFS=MCFSFeatures,Boruta=Linda_MCFSFeatures)
intersection=plotVenn(geneLS,list(c("darkmagenta", "darkblue"),c(0.5,0.5),2,2, c("TARGET-MCFS","LC-MCFS")))
intersection_annotated=getAnnotatedGenes(genes,intersection$`MCFS:Boruta`)

#-------------------------------Predict Normals------------------------------------------
#You have to run normalizationPediatrics first
Normals=Linda_GE_Classifier[which(grepl("BM.*",rownames (Linda_GE_Classifier))),]
decision_Normals=rep("Diagnosis",dim(Normals)[1])
Normals=cbind(Normals,decision_Normals)
Linda_GE_Classifier_temp=Normals[,c(Linda_MCFSFeatures[1:10],"decision_Normals")]
filteredMCFS_Linda=getAnnotatedGenes(genes,colnames(Linda_GE_Classifier_temp)[1:10])
annotation_Linda=annotateInOrder(filteredMCFS_Linda,Linda_classifierMCFSgenes)
colnames(Linda_GE_Classifier_temp)=append(annotation_Linda,"decision") 
predictClass(Linda_GE_Classifier_temp, Linda_resultRosettaMCFSGenetic$main)

#Checking runs of Lindaunmatched6 and Lindamatched3
x=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindamatched3/Lindamatched3__RI.csv",50)
y=FilterFeatures("Linda_Cohort_Results/MCFSoutput/Lindaunmatched6/Lindaunmatched6__RI.csv",50)
