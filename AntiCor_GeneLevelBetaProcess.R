
library(easypackages)



libraries("tidyverse","dplyr","GSVA","GSEABase","msigdbr","singscore","patchwork","cowplot","cluster","dendextend","factoextra",
          "ComplexHeatmap","presto","clustree","SummarizedExperiment","ggrepel","ggplot2","ggplotify","gridExtra","stringr","magrittr","forcats",
          "Rtsne","factoextra","NbClust","survival","survminer","scran","limma","RColorBrewer","viridis","colorRamp2","circlize","egg","ggmin") 

source("/data/asanigaris/Complete_Deconvolution/MB_Methylation/NCI_Methylation/MB_data/AntiCorrelatedGeneProbes.R")

message("Reading files...!")
data_dir="/data/asanigaris/Complete_Deconvolution/MB_Methylation"
cbtn_embl_beta<-readRDS(file.path(data_dir,"cbtn_embl_betaValues.rds"))
cbtn_beta<-read.csv(file.path(data_dir,"CBTN_ProbBetaVal_selected_Sample.csv"), header=TRUE, row.names = 1)
cell_data<-readRDS(file.path(data_dir,"allDattaSixCelltypeExp.rds"))
NCI_CBTN_beta<-readRDS("./Methylation_Data_processed/NCI_combined_plate_BetaValue.rds")
cell_deg_files<-readRDS(file.path(data_dir,"cellStateDEGlist.rds"))
message("Reading all files Done..!")

main_dir<-file.path(getwd(),"/Methylation_Data_processed/CBTN_EMBL_AntiCorProbeGeneResults")
if (!file.exists(main_dir)) {
  dir.create(main_dir, showWarnings = FALSE)
}

message("Anti coorelation process satarted....!")

# 1. Identify anti-cor probes
results_list <- identify_anticorr_probes(
  beta_file       = cbtn_embl_beta,
  cell_data_file  = cell_data,
  deg_files_list  = cell_deg_files,
  GeneNum         = c(500, 1000, 2000, 3000, 5000),
  Cells           = c("Malignant","Myeloid"),
  cor_cutoff      = -0.1,
  outdir          = file.path(main_dir,"AntiCorrelatedGeneProbes")
)

message("Gene level beta values process started...!")

# 2A. Option A: Use the in-memory 'results_list'
generate_gene_level_beta(
  beta_file      = NCI_CBTN_beta, 
  Cells          = c("Malignant","Myeloid"),
  GeneNum        = c(500,1000, 2000, 3000, 5000),
  cor_values     = c(-0.1, -0.2, -0.3, -0.4),
  outdir         = file.path(main_dir,"Genelevel_beta"),
  antiCorResults = results_list
)

message("Anticorrelation and Gene level beta value calculation is Done...!")





#dir<-"/data/asanigaris/Complete_Deconvolution/MB_Methylation/NCI_Methylation/MB_data/Methylation_Data_processed/AntiCorProbeGeneResults"
deg_dir<-"/data/asanigaris/Complete_Deconvolution/MB_Methylation"

Cells<-c("Malignant", "Myeloid")
for (Cell in Cells){
  tryCatch({
    indir<-file.path(main_dir, "Genelevel_beta")
    files<-list.files(path = file.path(indir,Cell),pattern = "*Gene_TSS1500_Median_betaValues_.*\\.csv", recursive = TRUE)
    if (Cell=="Malignant"){
      malig<-read.csv(file.path(deg_dir,"MalignantClusterNamedProteinCodingDEGs.csv"),header = T)
    }else if(Cell=="Myeloid"){
      Myel<-read.csv(file.path(deg_dir,"MyeloidClusterNamedProteinCodingDEGs.csv"), header = T)
    }
    for (file in files) {
      GeneNum<-sapply(strsplit(file,"Gene"), "[",1)
      GeneNum<-sapply(strsplit(GeneNum,"state"), "[",2)
      GeneNum<-as.numeric(GeneNum)
      CorNum<-sapply(strsplit(file,"_"), "[",5)
      CorNum<-sapply(strsplit( CorNum,"[.]"), "[",1)
      if (Cell=="Malignant"){
        maligGenes <- malig %>% filter(logFC > 0.01 & pval < 0.05) %>%
          group_by(group) %>% arrange(desc(logFC), .by_group=T) %>% top_n(n = GeneNum, wt = logFC)
        maligGenes<- maligGenes %>% split(x=.$feature, f=.$group)
        mGenelist<- list()
        for (i in names(maligGenes)){
          mGenelist[[i]]<- GeneSet(unique(maligGenes[[i]]), setName =i, organism="H")
        }
      }else if(Cell=="Myeloid"){
        MyelGenes <- Myel %>% filter(logFC > 0.01 & pval < 0.05) %>%
          group_by(group) %>% arrange(desc(logFC), .by_group=T) %>% top_n(n = GeneNum, wt = logFC)
        MyelGenes<- MyelGenes %>% split(x=.$feature, f=.$group)
        myGenelist<- list()
        for (i in names(MyelGenes)){
          myGenelist[[i]]<- GeneSet(unique(MyelGenes[[i]]), setName =i, organism="H")
        }
      }
      
      CellState_meth<-read.csv(file.path(indir,Cell, file),header = T,row.names = 1)
      outdir<-paste0(Cell,"_MethylClusters")
      outdir<- file.path(indir, Cell ,outdir)
      #dir.create(outdir)
      outdir<- file.path(outdir, GeneNum)
      if (!dir.exists(outdir)) {
        dir.create(outdir, recursive = TRUE)
      }
      
      meta_data<-read.csv("./Methylation_Data_processed/MB_sample_filter.csv",header = T,row.names = 1)
      #meta_data<-meta_data[meta_data$Subgroup=="Group3"| meta_data$Subgroup=="Group4",]
      #meta_data<-meta_data[!duplicated(meta_data$NIH_labels),]
      #meta_data<-meta_data[!grepl("CBTN",meta_data$NIH_labels),]
      #meta_data<-meta_data[!is.na(meta_data$Subgroup),]
      # GSVA Score
      if (Cell=="Malignant"){
        gcolc <- GeneSetCollection(mGenelist)
      } else if (Cell=="Myeloid"){
        gcolc <- GeneSetCollection(myGenelist)
      }
      
      ## GSVA soring method
      gsvaPar <- gsvaParam(data.matrix(CellState_meth), gcolc,
                           minSize=1, maxSize=5000, maxDiff = FALSE)
      SampScore <- gsva(gsvaPar)
      colnames(SampScore)<-gsub("X","", colnames(SampScore))
      cluster<-data.frame(t(SampScore))
      SampScore<-data.matrix(scale(SampScore))
      SampScore<- data.matrix(t(scale(t(SampScore))))
      
      if(Cell=="Malignant"){
        cluster$Malignant.Cluster <- apply(cluster, 1, function(x) {
          colnames(cluster)[which.min(x)]})
        cut_avg<-cluster[,"Malignant.Cluster", drop=F]
        meta_data<-merge(meta_data,cut_avg, by=0,all=F)
        rownames(meta_data)<-meta_data$Row.names
        meta_data<- meta_data[order(meta_data$Malignant.Cluster),]
        SampScore<-SampScore[,colnames(SampScore) %in% rownames(meta_data)]
        SampScore<-SampScore[,match(rownames(meta_data),colnames(SampScore))]
        save(SampScore,file =paste0(outdir,"/",Cell,"_",CorNum,"_betaScore.RData"))
        col = list(Malignant.Cluster=c("Neuronal"="#FF7F00","Neural.Crest"="#66A61E","Proliferative"="#C71587","Photoreceptor"="#1F78B4"),
                   #Malignant.Cluster=c("Nueronal1"="#FF7F00","Nueronal2"="#66A61E","Nueronal3"="#C71587","Photoreceptor"="#1F78B4"),
                   Subgroup=c("Group 3"="chartreuse3" ,"Group 4"= "deeppink2","Grp3/Grp4"="yellow3","SHH"="#FF7F00", "WNT"="#1F78B4"),
                   Subtype=c("MB_G34_I"="#004586","MB_G34_II"="#FF420E","MB_G34_III"="#FFD320","MB_G34_IV"="#579D1C","MB_G34_V"="#7E0021","MB_G34_VI"="#83CAFF","MB_G34_VII"="#314004","MB_G34_VIII"="#AECF00", 
                             "MB_SHH_1"="#4B1F6F","MB_SHH_2"="#FF950E","SHH_3"="#C5000B","SHH_4"="#DD4477","WNT"="#990099"),
                   Sex=c("F"="#17BECF","M"="#D62728")
                   #Age_group=c("Infants(0-3)"="#265DAB","Children(4-9)"="#059748","Toddlers(10-16)"="#00A2B3","Adults(17-40)"="#E5126F","Old(41-80)"="#7B3A96")
                   #Age_group=c("infant(0-3)"="#265DAB","child(4-16)"="#059748","adult(17-40)"="#E5126F","old(41-80)"="#7B3A96")
        )
        
        
        column_ha <- HeatmapAnnotation(Malignant.Cluster=as.factor(meta_data$Malignant.Cluster),
                                       Subgroup=as.factor(meta_data$Subgroup),
                                       Subtype=as.factor(meta_data$Subtype),
                                       Sex=as.factor(meta_data$Sex),
                                       #Age_group= as.factor(meta_data$Age_group),
                                       #show_annotation_name = c(Cluster = F,IDH.status=F ,Grade=F),
                                       annotation_legend_param = list(legend_direction = "horizontal", nrow = 4, by_row = F,
                                                                      title_position = "topcenter", title_gp = gpar(fontsize = 11)),
                                       col = col)
      }else if(Cell=="Myeloid"){
        cluster$Myeloid.Cluster <- apply(cluster, 1, function(x) {
          colnames(cluster)[which.min(x)]})
        cut_avg<-cluster[,"Myeloid.Cluster", drop=F]
        meta_data<-merge(meta_data,cut_avg, by=0,all=F)
        rownames(meta_data)<-meta_data$Row.names
        meta_data<- meta_data[order(meta_data$Myeloid.Cluster),]
        meta_data<- meta_data %>% mutate(Myeloid.Cluster=case_when(meta_data$Myeloid.Cluster=="Mac.M0"~"Macrophage.M0",
                                                                   meta_data$Myeloid.Cluster=="Mac.M1"~"Macrophage.M1",
                                                                   meta_data$Myeloid.Cluster=="Mac.M2"~"Macrophage.M2"))
        SampScore<-SampScore[,colnames(SampScore) %in% rownames(meta_data)]
        SampScore<-SampScore[,match(rownames(meta_data),colnames(SampScore))]
        save(SampScore,file =paste0(outdir,"/",Cell,"_",CorNum,"_betaScore.RData"))
        col = list(
          Myeloid.Cluster=c("Macrophage.M0" ="#FFA206","Macrophage.M1"="#2F4B7C","Macrophage.M2" ="#E31A1C"),
          #Myeloid.Cluster=c("Mac.M1"="#2F4B7C","Mac.M0"="#FFA206","Mac.M2"="#E31A1C"),
          Subgroup=c("Group 3"="chartreuse3" ,"Group 4"= "deeppink2","Grp3/Grp4"="yellow3","SHH"="#FF7F00", "WNT"="#1F78B4"),
          Subtype=c("MB_G34_I"="#004586","MB_G34_II"="#FF420E","MB_G34_III"="#FFD320","MB_G34_IV"="#579D1C","MB_G34_V"="#7E0021","MB_G34_VI"="#83CAFF","MB_G34_VII"="#314004","MB_G34_VIII"="#AECF00", 
                    "SHH_1"="#4B1F6F","SHH_2"="#FF950E","SHH_3"="#C5000B","SHH_4"="#DD4477","WNT"="#990099"),
          Sex=c("F"="#17BECF","M"="#D62728")
          #Age_group=c("Infants(0-3)"="#265DAB","Children(4-9)"="#059748","Toddlers(10-16)"="#00A2B3","Adults(17-40)"="#E5126F","Old(41-80)"="#7B3A96")
          #Age_group=c("infant(0-3)"="#265DAB","child(4-16)"="#059748","adult(17-40)"="#E5126F","old(41-80)"="#7B3A96")
        )
        
        
        column_ha <- HeatmapAnnotation(Myeloid.Cluster=as.factor(meta_data$Myeloid.Cluster),
                                       Subgroup=as.factor(meta_data$Subgroup),
                                       Subtype=as.factor(meta_data$Subtypes),
                                       Sex=as.factor(meta_data$Sex),
                                       #Age_group= as.factor(meta_data$Age_group),
                                       #show_annotation_name = c(Cluster = F,IDH.status=F ,Grade=F),
                                       annotation_legend_param = list(legend_direction = "horizontal", nrow = 4, by_row = F,
                                                                      title_position = "topcenter", title_gp = gpar(fontsize = 11)),
                                       col = col)
      }else if(Cell=="lymphoid"){
        cluster$Lymphoid.Cluster <- apply(cluster, 1, function(x) {
          colnames(cluster)[which.min(x)]})
        cut_avg<-cluster[,"Lymphoid.Cluster", drop=F]
        meta_data<-merge(meta_data,cut_avg, by=0,all=F)
        rownames(meta_data)<-meta_data$Row.names
        meta_data<- meta_data[order(meta_data$Lymphoid.Cluster),]
        SampScore<-SampScore[,colnames(SampScore) %in% rownames(meta_data)]
        SampScore<-SampScore[,match(rownames(meta_data),colnames(SampScore))]
        col = list(Lymphoid.Cluster=c("CD8Tcell"="#8965B3","CD4Tcell"="#2879CD","NKTcell"="#EE5761","CD8TRM"="#8DB032"),
                   Subgroup=c("Group3"="chartreuse3" ,"Group4"= "deeppink2","Grp3/Grp4"="yellow3","SHH"="#FF7F00", "WNT"="#1F78B4"),
                   Subtype=c("G34_I"="#004586","G34_II"="#FF420E","G34_III"="#FFD320","G34_IV"="#579D1C","G34_V"="#7E0021","G34_VI"="#83CAFF","G34_VII"="#314004","G34_VIII"="#AECF00", 
                             "SHH_1"="#4B1F6F","SHH_2"="#FF950E","SHH_3"="#C5000B","SHH_4"="#DD4477","WNT"="#990099"),
                   Sex=c("F"="#17BECF","M"="#D62728"),
                   #Age_group=c("Infants(0-3)"="#265DAB","Children(4-9)"="#059748","Toddlers(10-16)"="#00A2B3","Adults(17-40)"="#E5126F","Old(41-80)"="#7B3A96")
                   Age_group=c("infant(0-3)"="#265DAB","child(4-16)"="#059748","adult(17-40)"="#E5126F","old(41-80)"="#7B3A96"))
        
        
        column_ha <- HeatmapAnnotation(Lymphoid.Cluster=as.factor(meta_data$Lymphoid.Cluster),
                                       Subgroup=as.factor(meta_data$Subgroup),
                                       Subtype=as.factor(meta_data$Subtype),
                                       Sex=as.factor(meta_data$Sex),
                                       Age_group= as.factor(meta_data$Age_group),
                                       #show_annotation_name = c(Cluster = F,IDH.status=F ,Grade=F),
                                       annotation_legend_param = list(legend_direction = "horizontal", nrow = 4, by_row = F,
                                                                      title_position = "topcenter", title_gp = gpar(fontsize = 11)),
                                       col = col)
      }
      write.csv(meta_data, file.path(paste0(outdir,"/",Cell,"_",CorNum,"_metaCluster_Info.csv")), row.names = T)
      ## hetamap
      hp<- Heatmap(SampScore, name = "Values",
                   top_annotation = column_ha,
                   show_column_names = FALSE, 
                   cluster_columns=FALSE,
                   #cluster_rows = FALSE,
                   #cluster_columns = color_branches(col_dend, k =  5),
                   #clustering_distance_columns = "pearson",
                   #clustering_distance_columns = "spearman",
                   cluster_column_slices = FALSE,
                   show_column_dend = FALSE,
                   show_row_dend = FALSE,
                   column_split = meta_data$Malignant.Cluster,
                   #column_split =5,
                   #column_km = 5,
                   column_title = NULL,
                   col = colorRamp2(breaks = c(-2,-1, 0,1,2), colors = viridis(5)), 
                   heatmap_legend_param = list(legend_direction = "horizontal",color_bar = "continuous", #border = "black", 
                                               title = "Score", title_position = "topcenter", title_gp = gpar(fontsize = 14), legend_width = unit(1.5, "in")),
                   use_raster = T, raster_device = c("png")
      )
      if(Cell=="Malignant"){
        #CluColor=c("Nueronal1"="#FF7F00","Nueronal2"="#66A61E","Nueronal3"="#C71587","Photoreceptor"="#1F78B4")
        CluColor=c("Neuronal"="#FF7F00","Neural.Crest"="#66A61E","Proliferative"="#C71587","Photoreceptor"="#1F78B4")
      }else if(Cell=="Myeloid"){
        CluColor=c("Macrophage.M0" ="#FFA206","Macrophage.M1"="#2F4B7C","Macrophage.M2" ="#E31A1C")
        #CluColor=c("Mac.M1"="#2F4B7C","Mac.M0"="#FFA206","Mac.M2"="#E31A1C")
      }else if(Cell=="lymphoid"){
        CluColor=c("CD8Tcell"="#8965B3","CD4Tcell"="#2879CD","NKTcell"="#EE5761","CD8TRM"="#8DB032")
      }
      colnames(meta_data)[10]<-"Cluster"
      d2 <- meta_data %>% 
        group_by(Subgroup,Cluster) %>% 
        summarise(count = n()) %>% 
        mutate(percent=100*count/sum(count))
      d2<-d2[!is.na(d2$Subgroup),]
      mylColor=c("Macrophage.M0" ="#FFA206","Macrophage.M1"="#2F4B7C","Macrophage.M2" ="#E31A1C")
      #lymColor=c("CD8Tcell"="#8965B3","CD4Tcell"="#2879CD","NKTcell"="#EE5761","CD8TRM"="#8DB032")
      d2<-data.frame(d2)
      d2<-d2 %>% mutate(group=case_when(d2$Subgroup=="Group 3"~"GP3",
                                        d2$Subgroup=="Group 4"~"GP4"))
      pltNum=ggplot(d2, aes(x=group,y=percent , fill=forcats::fct_rev(Cluster)))+
        geom_col() + theme_bw()+
        xlab("Subgroup")+
        labs(y ="fraction of cell states (%)", tag = "b")+
        scale_y_continuous(labels=function(x) paste0(x,"%"))+
        scale_fill_manual(values =CluColor)+
        #theme_article()+
        coord_flip()+
        theme(legend.position = "none",legend.title = element_blank(),plot.tag = element_text(face = "bold", size = 20),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(),
              axis.line = element_line(colour = "black"),axis.text.x  = element_text(color="black",size=10),axis.title.x.bottom = element_blank(),plot.margin = margin(0,0,-0.5,0),
              axis.text=element_text(size=6.5, face="bold"),axis.text.y.left = element_text(face = "bold", size = 12) ,axis.title=element_text(size=12,face="bold"))
      label<- sort(unique(meta_data$Cluster))
      survPlt<-ggsurvplot(survfit(Surv(OS_months, OS_status)~Cluster, data=meta_data), 
                          xlab="Months",legend.title="Malignant clusters",
                          pval=TRUE, pval.size=8,legend.labs=label,
                          #pval.coord = c(125, 0.85),
                          palette = CluColor)
      survPlt$plot<-survPlt$plot+labs(tag = "c")+
        guides(colour = guide_legend( ncol = 2))+
        theme(legend.position = "none",plot.tag = element_text(face = "bold", size = 20), legend.title = element_blank(),
              legend.text = element_text(face = "bold", size = 12), axis.text.x = element_text(face = "bold", size = 12),
              axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(face = "bold", size = 14),
              axis.title.y.left = element_text(face = "bold", size = 14))
      
      pltfsuv<-wrap_plots(list(pltNum, survPlt$plot), ncol = 1,heights  = c(0.5,1))
      heatmap_grob <- grid.grabExpr(draw(hp,heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                                         adjust_annotation_extension = T, merge_legends = T))
      Hmap<-as.ggplot(heatmap_grob)
      Hmap<-Hmap+labs(tag = "a")+
        guides(colour = guide_legend(label.theme = element_text(size=10, face = "bold")))+
        theme(plot.tag = element_text(face = "bold", size = 20))
      HmPlt<-wrap_plots(list(Hmap, pltfsuv),ncol = 2, widths = c(2,1.2))&plot_annotation(tag_levels = "a")
      #saveRDS(maligPlt,file = file.path(getwd(),paste0(Cell,"ClusterPlots.rds")))
      'png(filename = file.path(paste0(outdir,"/",Cell,"_",CorNum,"_Clusterfigures.png")), width = 13.5, height = 6.5, units = "in", res = 1200,family = "ArialMT")
      print(HmPlt)
      dev.off()'
      
      subgroups<-c("Group 3","Group 4")
      SubSurv<-list()
      for (Sub in subgroups) {
        metaSub<-meta_data[meta_data$Subgroup==Sub,]
        metaSub<-metaSub[!is.na(metaSub$OS_status),]
        label<-sort(unique(metaSub$Cluster))
        surv<-ggsurvplot(survfit(Surv(OS_months, OS_status)~Cluster, data=metaSub), 
                         xlab="Months",ylab="Survival probability",legend.title="Malignant Clusters",title = Sub,
                         pval=TRUE, pval.size=8, legend.labs=label,palette =CluColor)
        surv$plot<-surv$plot+
          guides(color=guide_legend(ncol=2,title.position = "top"))+
          theme(legend.position = "bottom",legend.text = element_text(size = 14,face = "bold"),
                legend.title = element_blank(),text = element_text(size = 16,face = "bold"),plot.title = element_text(face = "bold", hjust = 0.5),
                axis.title.x.bottom = element_text(face = "bold"),axis.text.x = element_text(face = "bold"),axis.title.y.left = element_text(size = 16,face = "bold"),
                axis.text.y.left  = element_text(size = 16,face = "bold"),axis.line.x.bottom = element_line(linewidth =1),axis.line.y.left = element_line(linewidth = 1))
        SubSurv[[Sub]]<-surv$plot
      }
      SubCluSurv<-wrap_plots(SubSurv,ncol = length(SubSurv))&theme(legend.position = "none")
      figure4<-wrap_plots(list(HmPlt,SubCluSurv),nrow = 2, heights = c(1,0.5))
      png(filename = file.path(paste0(outdir,"/",Cell,"_",CorNum,"_Clusterfigures.png")), width = 12, height = 9.5, units = "in", res = 600,family = "ArialMT")
      print(figure4)
      dev.off()
    }
  } , error=function(e){message("Error in iteration for cell ", Cell, ": ", e$message)})
}


