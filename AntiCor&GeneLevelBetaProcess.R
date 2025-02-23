library(dplyr)
source("/data/asanigaris/Complete_Deconvolution/MB_Methylation/NCI_Methylation/MB_data/AntiCorrelatedGeneProbes.R")

message("Reading files...!")
data_dir="/data/asanigaris/Complete_Deconvolution/MB_Methylation"
cbtn_embl_beta<-readRDS(file.path(data_dir,"cbtn_embl_betaValues.rds"))
cell_data<-readRDS(file.path(data_dir,"allDattaSixCelltypeExp.rds"))
NCI_CBTN_beta<-readRDS("NCI_combined_plate_BetaValue.rds")
cell_deg_files<-readRDS(file.path(data_dir,"cellStateDEGlist.rds"))
message("Reading all files Done..!")

main_dir<-file.path(getwd(),"AntiCorProbeGeneResults")
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
