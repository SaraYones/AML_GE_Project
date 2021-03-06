savePDF=function(myplots,variable,filepath,discretizemethod="")
{
  
  
  
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}

#Write output for AML project (Gene expression)

writeOutput=function(path,date,folder1,folder2,clusteredRulesMCFS,recalculatedResultRosettalocal,rules,enrichment,enrichmentTitle,flag)
{
  #Create folder in the path
  ifelse(!dir.exists(paste(path,date,sep="")), dir.create(paste(path,date,sep="")), FALSE)
 # dir.create(paste(path,date,sep=""))
  ifelse(!dir.exists(paste(path,date,"/",folder1,sep="")), dir.create(paste(path,date,"/",folder1,sep="")), FALSE)
  ifelse(!dir.exists(paste(path,date,"/",folder1,"/",folder2,sep="")), dir.create(paste(path,date,"/",folder1,"/",folder2,sep="")), FALSE)  
  #dir.create(paste(path,date,folder,sep=""))
  temp=paste(path,date,"/",folder1,"/",folder2,sep="")
  write.csv(recalculatedResultRosettalocal,paste(temp,"/RulesAllGenes-",Sys.Date(),".csv",sep=""))
  print("Hello")
  saveLineByLine(rules,  paste(temp,"/NetworksAllGenes-",Sys.Date(),".txt",sep=""))
  print("Hi")
  graphics.off()
  my.plots=vector(1, mode='list');
  #svg(paste(temp,"/HeatMapAllGenes-",Sys.Date(),".svg",sep=""))
  pdf(paste(temp,"/HeatMapAllGenes-",Sys.Date(),".pdf",sep=""), onefile=TRUE)
  View(clusteredRulesMCFS)
   clusters=heatmap.F(t(clusteredRulesMCFS), colors=c('white','blue'),distmethod='pearson')
   print("Iam here")
  write.csv(clusters,paste(temp,"/Clusters-",Sys.Date(),".csv",sep = ""))
 # my.plots[[1]]=recordPlot()
  #savePDF(my.plots,paste("HeatMapAllGenes",Sys.Date()),paste(temp,"/",sep=""))
  
  
    dev.off()
  
plotEnrichment(enrichment,paste(temp,"/GOenrichment-",Sys.Date(),".pdf",sep=""),enrichmentTitle)
    plot("Afer Enrichment")
  #Merge output clusters with metadata exploratory
  #paste "TARGET-20-  to the names of the clusters then match it with the name of the clusters, get the values and add it to a new coloumn in metadata_exploratory
  if(flag==TRUE)
  {
  metadata_for_exploratory$clusters=clusters[match(rownames(metadata_for_exploratory),paste(rep("TARGET-20-",length(clusters)),names(clusters),sep=""))]
  #Order based on the order of the clusters 
  i1=match(paste(rep("TARGET-20-",length(clusters)),names(clusters),sep=""),rownames(metadata_for_exploratory))
  ordered_metadata_for_exploratory=metadata_for_exploratory[i1,]
  write.csv(ordered_metadata_for_exploratory,paste(temp,"/MergedMetaData-",Sys.Date(),".csv",sep=""))
  
  return(ordered_metadata_for_exploratory)
  
  
  }
}
