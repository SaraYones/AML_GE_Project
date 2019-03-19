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

writeOutput=function(path,date,folder,clusteredRulesMCFS,recalculatedResultRosettalocal,rules,enrichment,enrichmentTitle,flag)
{
  #Create folder in the path
  dir.create(paste(path,date,sep=""))
  ifelse(!dir.exists(paste(path,date,"/",folder,sep="")), dir.create(paste(path,date,"/",folder,sep="")), FALSE)
  
  #dir.create(paste(path,date,folder,sep=""))
  temp=paste(path,date,"/",folder,sep="")
  write.csv(recalculatedResultRosettalocal,paste(temp,"/RulesAllGenes-",Sys.Date(),".csv",sep=""))
  saveLineByLine(rules,  paste(temp,"/NetworksAllGenes-",Sys.Date(),".txt",sep=""))
  
  graphics.off()
  my.plots=vector(1, mode='list');
  #svg(paste(temp,"/HeatMapAllGenes-",Sys.Date(),".svg",sep=""))
  pdf(paste(temp,"/HeatMapAllGenes-",Sys.Date(),".pdf",sep=""), onefile=TRUE)
  View(clusteredRulesMCFS)
   clusters=heatmap.F(t(clusteredRulesMCFS), colors=c('white','white','blue','blue'),distmethod='pearson')
  write.csv(clusters,paste(temp,"/Clusters-",Sys.Date(),".csv",sep = ""))
 # my.plots[[1]]=recordPlot()
  #savePDF(my.plots,paste("HeatMapAllGenes",Sys.Date()),paste(temp,"/",sep=""))
  
  
    dev.off()
  
plotEnrichment(enrichment,paste(temp,"/GOenrichment-",Sys.Date(),".pdf",sep=""),enrichmentTitle)
    
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
