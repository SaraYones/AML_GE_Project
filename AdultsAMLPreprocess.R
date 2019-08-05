#Adults Preprocessing (Checking for confounders)


#Relapse unmatched

cases_Adults_Relapse1_unmatched<-read.xlsx("Linda_Cohort_Metadata/AdultsData/AdultsMetadata.xlsx", sheetName = "Relapse1")
#The coloumn names are marked with "." so we need to change to -
new_col_names=str_replace_all(rownames(Linda_GE_Classifier2), "\\.","-")
rownames(Linda_GE_Classifier2)<-new_col_names

files2_Linda<-read.xlsx("Linda_Cohort_Metadata/Svea/20190123_AML_Metadata.xls", sheetName = "Sheet1")

#Filtering Adults according to age If you want to check for the ages of the cases you include
ages=files2_Linda[which(!is.na(files2_Linda$Age.at.Sampling..Yrs.)),]
ages=ages[which(ages$Age.at.Sampling..Yrs.!="N/A"),]
ages$Age.at.Sampling..Yrs.=as.numeric(as.character(ages$Age.at.Sampling..Yrs.))

#remove all the control cases
Linda_GE_Classifier2=cbind.data.frame(Linda_GE_Classifier2,decision_Linda2)

Linda_GE_Classifier2=Linda_GE_Classifier2[-which(grepl("BM.*",rownames (Linda_GE_Classifier2))),]
#if i run the normalize with (run pipeline==FALSE) to remove all the nearzero cols

#because we normalized on the whole dataset so i decided not to remove until everything was normalized
remove_cols=nearZeroVar(Linda_GE_Classifier2)
if(length(remove_cols)!=0)
  Linda_GE_Classifier2=Linda_GE_Classifier2[,-remove_cols]


#check confounders
metadata_Linda=read.xlsx("Linda_Cohort_Metadata/Metadata_RNA_Sara.xls", sheetName = "Sheet1")

metadata_Adults_Linda<-read.xlsx("Linda_Cohort_Metadata/Svea/20190123_AML_Metadata.xls", sheetName = "Sheet1")


#metadata_Linda[,"Sample.ID"] <-str_replace_all(metadata_Linda[,"Sample.ID"] , "_|-", ".")
metadata_Adults_Linda=metadata_Adults_Linda[which(metadata_Adults_Linda[,"Sample.ID"] %in% rownames(Linda_GE_Classifier2) ),]

paste(getwd(),"/Linda_Adults_Cohort_Results/",sep="")

form_Linda_Adults <- ~ (1|Sex) + (1|Stage) + (1|Age.at.onset..Grouped.) +(1|Blast.count.....Grouped.)+(1|FAB.subtype)+(1|EFS..Grouped.)


checkVariableEffects(list(as.data.frame(Linda_GE_Classifier2[,1:dim(Linda_GE_Classifier2)[2]-1])),list(as.data.frame(metadata_Adults_Linda)),form_Linda,paste(getwd(),"/Linda_Adults_Cohort_Results/",sep=""),"variableEffectsAdults")






