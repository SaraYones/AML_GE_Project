if (!requireNamespace("BiocManager", quietly=TRUE))
install.packages("BiocManager")
BiocManager::install("topGO")


BiocManager::install("clusterProfiler", version = "3.8")
de <- classifierMCFSgenes
background<-colnames(TARGET_GE_Classifier)[-dim(TARGET_GE_Classifier)[2]]

gene.df <- bitr(de, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

x<-org.Hs.egENSEMBL
yy <- enrichGO(de,'org.Hs.eg.db',keyType= 'ENSEMBL', ont="BP",universe=background, pvalueCutoff=0.1)
filteredBoruta=getAnnotatedGenes(genes,BorutaFeatures)



library(org.Hs.eg.db)
data(geneList)
gene <- names(geneList)[abs(geneList) >= 1]
gene.df <- bitr(gene, fromType ="ENTREZID",
                toType = c("ENSEMBL"),
                OrgDb = org.Hs.eg.db)


GOenrichment=function(features,background,flag)
{
  
  if (flag=="ENSEMBL")
  {
    
  ego <- enrichGO(gene          = features,
                  universe      = background,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = 'ENSEMBL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.01)
  }
  return(ego)
  
  
  
}



de=Linda_MCFSFeatures[1:700]
background=colnames(Linda_GE_Classifier)[1:length(Linda_GE_Classifier)-1]
ego <- enrichGO(gene          = de,
                universe      = background,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)


de=Integrated_Matched_MCFSFeatures[1:20]
background=colnames(Integrated_Matched_Classifier)[1:length(Integrated_Matched_Classifier)-1]
ego <- enrichGO(gene          = de,
                universe      = background,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)

de=MCFSFeatures[1:500]
background=colnames(TARGET_GE_Classifier)[1:length(TARGET_GE_Classifier)-1]
ego <- enrichGO(gene          = de,
                universe      = background,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.01)



data(geneList, package="DOSE")
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

geneList.df <- bitr(names(geneList), fromType = "ENTREZID",
                toType = c("ENSEMBL", "SYMBOL"),
                OrgDb = org.Hs.eg.db)


head(gene.df)
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

ego2 <- enrichGO(gene         = gene.df$ENSEMBL,
                 universe      = geneList.df$ENSEMBL,
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "CC",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

yy <- enrichGO(gene         = de,
              universe      = background,
               OrgDb         = org.Hs.eg.db,
               keyType       = 'ENSEMBL',
               ont           = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)

ggo <- groupGO(gene     = de,
              # universe      = background,
               OrgDb    = org.Hs.eg.db,
               keyType       = 'ENSEMBL',
               ont      = "BP",
               level    = 3)

gene.df <- bitr(de, fromType = "ENSEMBL",
                toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)

#Try with EntrezID (was done by bioMart the conversion)

genes.with.id=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=de, mart= ensembl) 
trial=genes.with.id[which(!is.na(genes.with.id[which(genes.with.id$ensembl_gene_id %in% de),])),]
trial=trial[which(!is.na(trial$entrezgene)),]

genes.with.id.bg=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=background, mart= ensembl) 
trial.bg=genes.with.id.bg[which(!is.na(genes.with.id.bg[which(genes.with.id.bg$ensembl_gene_id %in% background),])),]
trial.bg=trial.bg[which(!is.na(trial.bg$entrezgene)),]


genes.with.id=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=de, mart= ensembl) 
trial=genes.with.id[which(!is.na(genes.with.id[which(genes.with.id$ensembl_gene_id %in% de),])),]
trial=trial[which(!is.na(trial$entrezgene)),]

genes.with.id.bg=getBM(attributes=c("ensembl_gene_id", "entrezgene"),filters ='ensembl_gene_id' ,values=background, mart= ensembl) 
trial.bg=genes.with.id.bg[which(!is.na(genes.with.id.bg[which(genes.with.id.bg$ensembl_gene_id %in% background),])),]
trial.bg=trial.bg[which(!is.na(trial.bg$entrezgene)),]


yy <- enrichGO(gene         = as.character(trial$entrezgene),
               universe      =as.character(trial.bg$entrezgene),
               OrgDb         = org.Hs.eg.db,
             #  keyType       = 'ENSEMBL',
               ont           = "CC",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.01,
               qvalueCutoff  = 0.05)

dotplot(yy, showCategory=30)
#Without using BG
x <- enrichDO(gene          = as.character(trial$entrezgene),
              ont           = "DO",
            #  keyType       = 'ENSEMBL',
              pvalueCutoff  = 0.01,
              pAdjustMethod = "BH",
              universe      = as.character(trial.bg$entrezgene),
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.01,
              readable      = FALSE)

#Try with gene symbols alone without any transformation


#Protein Coding Genes only

geneProtein=append(filteredMCFS[[1]]$geneID,filteredMCFS[[2]])
geneProteinbg=append(filter[[1]]$geneID,filter[[2]])

ego <- enrichGO(gene          = geneProtein,
             #   universe      = geneProteinbg,
                OrgDb         = org.Hs.eg.db,
                keyType       = 'ENSEMBL',
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)


#I didnt get any significant BP when using protein coding genes

