library(minfi)
library(psych)
library(tidyverse)

files<-list.files(pattern = "Bval.rds",recursive = TRUE,full.names = T)
b1<-readRDS(files[[1]])
b1<-data.frame(b1)
b2<-readRDS(files[[2]])
b2<-data.frame(b2)
b3<-readRDS(files[[3]])
b3<-data.frame(b3)
b3$Row.names<-rownames(b3)
b3<-b3 %>% dplyr::select(Row.names,everything())
b3$Row.names<-sapply(strsplit(b3$Row.names,"_"),"[",1)
b3<-b3[!duplicated(b3$Row.names),]
allB<-merge(b1,b2,by=0,all=F)
allB<-merge(allB,b3,by=1,all=F)
write.csv(allB,"NCI_BetaValues.csv",row.names = F)
NCI_beta<-read.csv("NCI_BetaValues.csv", header = T,row.names = 1)
NCI_beta<-t(NCI_beta)
NCI_beta<-data.frame(NCI_beta)
antiProbe<-read.csv(file.choose(), header = T)
antiProbe<-antiProbe[,c("cpg","feature")]
GeneBval<-merge(antiProbe,NCI_beta,by.x=1,by.y=0,all=F)
GeneBval<-GeneBval[,-1]
GeneBval<-aggregate(. ~ feature, data =GeneBval, FUN = median)

write.csv(GeneBval, "NCI_GP34_SixMyelodCellstatePredicted_TSS1500_Median_betaValues.csv", row.names = F)
adNCI_beta<-readRDS("betas.rds")
antiProbe<-read.csv(file.choose(), header = T)
antiProbe<-antiProbe[,c("cpg","feature")]
GeneBval<-merge(antiProbe,adNCI_beta,by.x=1,by.y=0,all=F)
GeneBval<-GeneBval[,-1]
GeneBval<-aggregate(. ~ feature, data =GeneBval, FUN = median)

write.csv(GeneBval, "NCI_GP34_newSampleMalignantCellstate1000Gene_TSS1500_Median_betaValues.csv", row.names = F)

