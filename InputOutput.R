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

writeOutput=function(path,date,folder,clusteredRulesMCFS)
{
  #Create folder in the path
  dir.create(paste(path,date,sep=""))
  ifelse(!dir.exists(paste(path,date,"/",folder,sep="")), dir.create(paste(path,date,"/",folder,sep="")), FALSE)
  
  #dir.create(paste(path,date,folder,sep=""))
  temp=paste(path,date,"/",folder,sep="")
  write.csv(recalculatedResultRosettaMCFS,paste(temp,"/RulesAllGenes-",Sys.Date(),".csv",sep=""))
  saveLineByLine(recalculatedResultRosettaMCFS,  paste(temp,"/NetworksAllGenes-",Sys.Date(),".txt",sep=""))
  
  my.plots=vector(1, mode='list');
  svg(paste(temp,"/HeatMapAllGenes-",Sys.Date(),".svg",sep=""))
  clusters=heatmap.F(t(clusteredRulesMCFS), colors=c('white','white','white','white','blue','blue'),distmethod='pearson')
  write.csv(clusters,paste(temp,"/Clusters-",Sys.Date(),".csv",sep = ""))
  my.plots[[1]]=recordPlot()
  savePDF(my.plots,paste("HeatMapAllGenes",Sys.Date()),paste(temp,"/",sep=""))
  #  dev.off()
  
  #Merge output clusters with metadata exploratory
  #paste "TARGET-20-  to the names of the clusters then match it with the name of the clusters, get the values and add it to a new coloumn in metadata_exploratory
  
  metadata_for_exploratory$clusters=clusters[match(rownames(metadata_for_exploratory),paste(rep("TARGET-20-",length(clusters)),names(clusters),sep=""))]
  #Order based on the order of the clusters 
  i1=match(paste(rep("TARGET-20-",length(clusters)),names(clusters),sep=""),rownames(metadata_for_exploratory))
  ordered_metadata_for_exploratory=metadata_for_exploratory[i1,]
  write.csv(ordered_metadata_for_exploratory,paste(temp,"/MergedMetaData-",Sys.Date(),".csv",sep=""))
  
  return(ordered_metadata_for_exploratory)
  
  
  
}
