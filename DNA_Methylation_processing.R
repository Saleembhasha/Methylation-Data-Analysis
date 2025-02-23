library(minfi)
library(GEOquery)
library(limma)
library(openxlsx)
library(tidyverse)
library(RColorBrewer)
library(missMethyl) # Can take a short time...
library(minfiData)
library(Gviz)
library(DMRcate)
library(DMRcatedata)
library(stringr)
#library(mCSEA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#library("optparse")



run_methyl<-function(targets, norm.method="preprocessRaw",prob.filter=FALSE,levels_beta=c("UCSC_RefGene_Name","TSS1500;TSS1500"),outdir="MethylResult",dm.Anls=FALSE){
  
  ## out dir
  outdir<-outdir
  if (!file.exists(outdir)) {
    dir.create(outdir)
  }
  
  ## read methylation idat files
  RGset <- read.metharray.exp(targets = targets, verbose = TRUE, force = TRUE)
  cat("reading idat file is completed!\n")
  
  ## pre processing
  if(norm.method=="preprocessRaw"){
    MSet <- preprocessRaw(RGset)
  }else if(norm.method=="preprocessIllumina"){
    MSet <- preprocessIllumina(RGset)
  }else if(norm.method=="preprocessSWAN"){
    MSet <- preprocessSWAN(RGset)
  }else if(norm.method=="preprocessFunnorm"){
    MSet <- preprocessFunnorm(RGset)
  }else if(norm.method=="preprocessQuantile"){
    MSet <- preprocessQuantile(RGset)
  }else if(norm.method=="preprocessQuantile"){
    MSet <- preprocessNoob(RGset)
  }
  cat("normalization process is completed!\n")
  ima <- annotation(MSet)[['array']]
  array_type =gsub("IlluminaHumanMethylation","",ima)
  array_type =gsub("v[0-9]","",array_type)
  ## get beta and M values
  RatioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
  saveRDS(RatioSet,file = file.path(outdir,"RatioSet.rds"))
  Gset <- mapToGenome(RatioSet)
  saveRDS(Gset,file = file.path(outdir,"Gset.rds"))
  Bval <- getBeta(Gset)
  saveRDS(Bval,file = file.path(outdir,"Bval.rds"))
  Mval <- getM(Gset)
  saveRDS(Mval,file = file.path(outdir,"Mval.rds"))
  CN <- getCN(Gset)
  cat("Beta and Mvalue converstin process is completed!\n")
  #save(RGset,MSet,RatioSet,Gset,Bval,Mval,CN, file = file.path(outdir,"DNAmethylationData.RData"))
  ##
  if(prob.filter==TRUE){
    ## filter# probe filtering
    message("probe filtering ...",Sys.time())
    amb.filter <- read.table(file.path("mnp_training-master","filter","amb_3965probes.vh20151030.txt"),header=F)
    epic.filter <- read.table(file.path("mnp_training-master","filter","epicV1B2_32260probes.vh20160325.txt"),header=F)
    snp.filter <- read.table(file.path("mnp_training-master","filter","snp_7998probes.vh20151030.txt"),header=F)
    xy.filter <- read.table(file.path("mnp_training-master","filter","xy_11551probes.vh20151030.txt"),header=F)
    rs.filter <- grep("rs",rownames(MSet))
    ch.filter <- grep("ch",rownames(MSet))
    
    # filter CpG probes
    remove <- unique(c(match(amb.filter[,1], rownames(MSet)),
                       match(epic.filter[,1], rownames(MSet)),
                       match(snp.filter[,1], rownames(MSet)),
                       match(xy.filter[,1], rownames(MSet)),
                       rs.filter,
                       ch.filter))
    
    MSet_filtered <- MSet[-remove,]
    RatioSet_filtered<- ratioConvert(MSet_filtered, what = "both", keepCN = TRUE)
    Gset_filtered <- mapToGenome(RatioSet_filtered)
    Bval_filtered <- getBeta(Gset_filtered)
    Mval_filtered <- getM(Gset_filtered)
    #save(RGset,MSet,RatioSet,Gset,Bval,Mval,CN, file = file.path(outdir,"DNAmethylationData_filtered.RData"))
  }
  
  ### annotation file
  if(array_type=="450K"){
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  }else if(array_type=="EPIC"){
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  }
  if(!is_empty(levels_beta)){
    group<-levels_beta
    ## get gene names form annotatoon data
    annoG<- data.frame(anno@listData)
    if(length(as.list(strsplit(group[2],";")[[1]]))>=2){
      grouplist<-c(group[2],sapply(strsplit(group[2],";"),"[",1))
      annoG<-data.frame(annoG[annoG$UCSC_RefGene_Group==grouplist[1]|annoG$UCSC_RefGene_Group==grouplist[2],])
    }else{
      annoG<-annoG[annoG$UCSC_RefGene_Group==group[2],]
    }
    annoG<-annoG[,c("Name", "UCSC_RefGene_Name")]
    annoG<-annoG[!annoG$UCSC_RefGene_Name=="",]
    annoG$UCSC_RefGene_Name<-sapply(strsplit(annoG$UCSC_RefGene_Name,";"),"[",1)
    rownames(annoG)<-annoG$Name
    annoG<-annoG[,"UCSC_RefGene_Name", drop=FALSE]
    
    ## calculate gene levels beta values
    beta_gene<-merge(annoG,beta,by=0,all=F)
    beta_gene<-beta_gene[,-1]
    beta_gene<-aggregate(. ~ UCSC_RefGene_Name, data =beta_gene, FUN = median)
    write.csv(beta_gene,"Gene_level_Beta_values.txt", row.names = F, quote = F,sep = "\t")
  }
  if(dm.Anls==FALSE){
    cat("not required differential methylation analysis")
  }else if(dm.Anls==prob){
    group <- factor(targets$Sample_Group)
    # use the above to create a design matrix
    design <- model.matrix(~0+group, data=targets)
    colnames(design) <- levels(group)[-1]
    # fit the actual linear model to the data
    fit <- lmFit(Mval_filtered, design)
    if(levels(group)>2){
      # fit the contrasts
      fit2 <- contrasts.fit(fit)
      # Rank genes
      fit2 <- eBayes(fit2)
      anno <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
      DMPs <- topTable(fit2, num=Inf, coef=1, genelist=anno)
      }else if(levels(group)>2){
        contMatrix <- contMatrix
        # fit the contrasts
        fit2 <- contrasts.fit(fit, contMatrix)
        # Rank genes
        fit2 <- eBayes(fit2)
        for (i in levels(group)) {
          anno <- ann450k[match(rownames(mVals),ann450k$Name),
                        c(1:4,12:19,24:ncol(ann450k))]
          DMPs <- topTable(fit2, num=Inf, coef=1, genelist=anno)
      }
      
    }
  }
}

## list options for CODEFACS
'option_list = list(
  make_option(c("-i", "--idatfiles"), type="character", default=NULL, 
              help="idat file dataframe", metavar="character"),
  make_option(c("-n", "--norm.method"), type="character", default="preprocessRaw", 
              help="Method for normalization", metavar="character"),
  make_option(c("-f", "--prob.filter"), type="character", default=FALSE, 
              help="probe filtering", metavar="character"),
  make_option(c("-n", "--levels_beta"), type="character", default=NULL, 
              help="gene levels of beta based on transcipt group ex.UCSC_RefGene_Name and TSS1500", metavar="character"),
  make_option(c("-e", "--extra_slurm_parameters"), type="character", default=NULL, 
              help="extra slurm parameters such like \"partition=XX, mem=50g, time=02:00:00\" [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="MethylResult", 
              help="output file name for beta and Mvalues", metavar="character"))'


## parse arguments
'opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)'

# read in the sample sheet for the experiment
'MB_meta<-read.csv(file.choose(),header=T)
MB_meta$sample_id<-paste0("CBTN-",MB_meta$sample_id)
MB_meta<-MB_meta[!MB_meta$subgroup=="MB, To be classified",]
SampleID<-MB_meta[,c("Kids_First_Biospecimen_ID", "Sample_Name","subgroup")]
targets <- read.metharray.sheet(getwd(), pattern="Sample_Sheet.csv")
targets<-targets[targets$Sample_Name %in% MB_meta$sample_id,]'


'targets <-list.files(pattern = "*idat.gz",full.names = T)
targets<-data.frame(targets)
targets$Basename<-sapply(strsplit(targets$targets,"_Grn|_Red"),"[",1)
targets$Sentrix_id<-gsub("./","",targets$Basename)
targets$SampleID<-sapply(strsplit(targets$Sentrix_id,"_"),"[",1)
targets<-targets[,-1]
targets<-targets[!duplicated(targets$Basename),]
methyl<-run_methyl(targets = targets,levels_beta=NULL)'

sample<-read.csv("MB_samples.csv",header = T)
for (plate in c("HumanMethylationEPICv2","HumanMethylationEPIC")) {
  targets<-sample[sample$Platform_methy==plate,]
  targets<-targets[,c("Sentrix_ID","Sentrix_position","idat_filename")]
  targets$Basename <- paste0("./",targets$idat_filename)
  methyl<-run_methyl(targets = targets,levels_beta=NULL,outdir = plate)
}
