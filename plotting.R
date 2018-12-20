plotPCAmeta=function(GE,metadata,variable,filepath,discretizemethod)
{
  
  
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    
    
    temp=metadata[[i]]
    temp=temp[,variable]
    
    if(variable=="age:"|variable=="Age.at.Diagnosis.in.Days")
    {
      temp=discretizeAge(temp,discretizemethod)
      
      
    }
    
    
    
    #print(GE[i])
    #print(temp[,variable])
    my.plots[[i]]=plotPCA(filepath,GE[[i]],temp,variable)
    
  }
  
  graphics.off()
  
  pdf(paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}


plotPCA= function(filepath,GeneExpressionMatrixlocal,Groups,variable){
  
  myplots=NULL
  #pdf(filepath)
  # log transform 
  # apply PCA - scale. = TRUE is highly 
  # advisable, but default is FALSE. 
  
  #ir.pca <- prcomp(GeneExpressionMatrixlocal,
  #                center = TRUE,
  #               scale. = TRUE) 
  ir.pca <- prcomp(GeneExpressionMatrixlocal,
                   center = TRUE) 
  
  
  g <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
                groups =Groups, ellipse = FALSE, 
                circle = FALSE,var.axes = FALSE)
  g <- g + scale_color_discrete(name = variable)
  g <- g + theme(legend.direction = 'horizontal', 
                 legend.position = 'top')
  #print(g)
  #dev.off()
  plot(g)
  myplots=recordPlot()
  return(myplots)
}

plotVenn=function(geneLS,parameters)
{
  
  VENN.LIST <- geneLS
  venn.plot <- venn.diagram(VENN.LIST , NULL, fill=parameters[[1]], alpha=parameters[[2]], cex = as.numeric(parameters[[3]]), cat.fontface=4, category.names=parameters[[5]], main="Gene Lists")
  
  # To plot the venn diagram we will use the grid.draw() function to plot the venn diagram
  grid.draw(venn.plot)
  
  # To get the list of gene present in each Venn compartment we can use the gplots package
  require("gplots")
  
  a <- venn(VENN.LIST, show.plot=FALSE)
  
  # You can inspect the contents of this object with the str() function
  inters <- attr(a,"intersections")
  return(inters)
  
  
}


plotDistributionCategory=function(GE,titles,variable,filepath,discretizemethod)
{  
  
  print(length(GE))
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    op <- par(mar = c(16,7,7,4) + 0.1)
    end_point = 0.5 + length(GE[[i]]) #+ length(GE[[i]])-1 #this is the line which does the trick (together with barplot "space = 1" parameter)
    
    plt=barplot(table(GE[[i]]),xlab = "",xaxt = "n",space=1,cex.axis = 0.5 ,cex.names = 0.5,cex.main=0.8,main=titles[i])
    par(op)
    text( plt, par("usr")[3]-0.25, 
          srt = 60, adj= 1, xpd = TRUE,
          labels = names(table(GE[[i]])), cex=0.65)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
    
  }
  
  graphics.off()
  # update_geom_default("point", list(size=1))
  # theme_set(theme_grey(base_size=6))
  sc <- c(0.5,0.75,1)
  
  pdf(width=7*sc[1],height=7*sc[1],pointsize=1,paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}

#--------Plotting Functions--------------------
plotDistributionCategory=function(GE,titles,variable,filepath,discretizemethod)
{  
  
  print(length(GE))
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    op <- par(mar = c(16,7,7,4) + 0.1)
    end_point = 0.5 + length(GE[[i]]) #+ length(GE[[i]])-1 #this is the line which does the trick (together with barplot "space = 1" parameter)
    
    plt=barplot(table(GE[[i]]),xlab = "",xaxt = "n",space=1,cex.axis = 0.5 ,cex.names = 0.5,cex.main=0.8,main=titles[i])
    par(op)
    text( plt, par("usr")[3]-0.25, 
          srt = 60, adj= 1, xpd = TRUE,
          labels = names(table(GE[[i]])), cex=0.65)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
    
  }
  
  graphics.off()
  # update_geom_default("point", list(size=1))
  # theme_set(theme_grey(base_size=6))
  sc <- c(0.5,0.75,1)
  
  pdf(width=7*sc[1],height=7*sc[1],pointsize=1,paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}


#--------Plotting Functions--------------------
plotDistributionCategory=function(GE,titles,variable,filepath,discretizemethod)
{  
  
  print(length(GE))
  
  my.plots <- vector(length(GE), mode='list')
  
  
  for(i in 1:length(GE)){
    op <- par(mar = c(16,7,7,4) + 0.1)
    end_point = 0.5 + length(GE[[i]]) #+ length(GE[[i]])-1 #this is the line which does the trick (together with barplot "space = 1" parameter)
    
    plt=barplot(table(GE[[i]]),xlab = "",xaxt = "n",space=1,cex.axis = 0.5 ,cex.names = 0.5,cex.main=0.8,main=titles[i])
    par(op)
    text( plt, par("usr")[3]-0.25, 
          srt = 60, adj= 1, xpd = TRUE,
          labels = names(table(GE[[i]])), cex=0.65)
    
    myplots=recordPlot()
    my.plots[[i]]=myplots
    
  }
  
  graphics.off()
  # update_geom_default("point", list(size=1))
  # theme_set(theme_grey(base_size=6))
  sc <- c(0.5,0.75,1)
  
  pdf(width=7*sc[1],height=7*sc[1],pointsize=1,paste(filepath,variable,"-",discretizemethod,".pdf",sep=""), onefile=TRUE)
  for (my.plot in my.plots) {
    replayPlot(my.plot)
  }
  graphics.off()
  
}



