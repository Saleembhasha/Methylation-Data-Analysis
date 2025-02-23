library(minfi)
library(GEOquery)
library(limma)
library(openxlsx)
library(tidyverse)
library(RColorBrewer)
#library(missMethyl) # Can take a short time...
#library(minfiData)
library(Gviz)
library(DMRcate)
#library(DMRcatedata)
library(stringr)
#library(mCSEA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
#library("optparse")



run_methyl <- function(targets, norm.method = "preprocessRaw", prob.filter = FALSE, levels_beta = c("UCSC_RefGene_Name", "TSS1500;TSS1500"), outdir = "MethylResult", dm.Anls = FALSE) {
  
  # Create output directory
  if (!file.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }
  
  # Read methylation IDAT files
  RGset <- read.metharray.exp(targets = targets, verbose = TRUE, force = TRUE)
  cat("Reading IDAT files is completed!\n")
  
  # Preprocessing
  if (norm.method == "preprocessRaw") {
    MSet <- preprocessRaw(RGset)
  } else if (norm.method == "preprocessIllumina") {
    MSet <- preprocessIllumina(RGset)
  } else if (norm.method == "preprocessSWAN") {
    MSet <- preprocessSWAN(RGset)
  } else if (norm.method == "preprocessFunnorm") {
    MSet <- preprocessFunnorm(RGset)
  } else if (norm.method == "preprocessQuantile") {
    MSet <- preprocessQuantile(RGset)
  } else if (norm.method == "preprocessNoob") {
    MSet <- preprocessNoob(RGset)
  }
  cat("Normalization process is completed!\n")
  
  # Get array type
  ima <- annotation(MSet)[['array']]
  array_type <- gsub("IlluminaHumanMethylation", "", ima)
  array_type <- gsub("v[0-9]", "", array_type)
  cat(paste0("array type is : ",array_type, "\n"))
  # Get beta and M values
  RatioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
  saveRDS(RatioSet, file = file.path(outdir, "RatioSet.rds"))
  Gset <- mapToGenome(RatioSet)
  saveRDS(Gset, file = file.path(outdir, "Gset.rds"))
  Bval <- getBeta(Gset)
  saveRDS(Bval, file = file.path(outdir, "Bval.rds"))
  Mval <- getM(Gset)
  saveRDS(Mval, file = file.path(outdir, "Mval.rds"))
  CN <- getCN(Gset)
  cat("Beta and M-value conversion process is completed!\n")
  
  # Probe filtering (if required)
  if (prob.filter) {
    message("Performing probe filtering...", Sys.time())
    amb.filter <- read.table(file.path("mnp_training-master", "filter", "amb_3965probes.vh20151030.txt"), header = FALSE)
    epic.filter <- read.table(file.path("mnp_training-master", "filter", "epicV1B2_32260probes.vh20160325.txt"), header = FALSE)
    snp.filter <- read.table(file.path("mnp_training-master", "filter", "snp_7998probes.vh20151030.txt"), header = FALSE)
    xy.filter <- read.table(file.path("mnp_training-master", "filter", "xy_11551probes.vh20151030.txt"), header = FALSE)
    rs.filter <- grep("rs", rownames(MSet))
    ch.filter <- grep("ch", rownames(MSet))
    
    # Filter CpG probes
    remove <- unique(c(
      match(amb.filter[, 1], rownames(MSet)),
      match(epic.filter[, 1], rownames(MSet)),
      match(snp.filter[, 1], rownames(MSet)),
      match(xy.filter[, 1], rownames(MSet)),
      rs.filter,
      ch.filter
    ))
    
    MSet_filtered <- MSet[-remove, ]
    RatioSet_filtered <- ratioConvert(MSet_filtered, what = "both", keepCN = TRUE)
    Gset_filtered <- mapToGenome(RatioSet_filtered)
    Bval_filtered <- getBeta(Gset_filtered)
    Mval_filtered <- getM(Gset_filtered)
  }
  
  # Annotation file
  if (array_type == "450K") {
    anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  } else if (array_type == "EPIC") {
    anno <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  }
  
  # Calculate gene-level beta values
  if (!is_empty(levels_beta)) {
    group <- levels_beta
    annoG <- anno
    if (length(strsplit(group[2], ";")[[1]]) >= 2) {
      grouplist <- c(group[2], sapply(strsplit(group[2], ";"), "[", 1))
      annoG <- annoG[annoG$UCSC_RefGene_Group %in% grouplist, ]
    } else {
      annoG <- annoG[annoG$UCSC_RefGene_Group == group[2], ]
    }
    annoG <- annoG[, c("Name", "UCSC_RefGene_Name")]
    annoG <- annoG[annoG$UCSC_RefGene_Name != "", ]
    annoG$UCSC_RefGene_Name <- sapply(strsplit(annoG$UCSC_RefGene_Name, ";"), "[", 1)
    rownames(annoG) <- annoG$Name
    annoG <- annoG[, "UCSC_RefGene_Name", drop = FALSE]
    
    # Merge with beta values
    beta_gene <- merge(annoG, Bval, by = 0, all = FALSE)
    beta_gene <- beta_gene[, -1]
    beta_gene <- aggregate(. ~ UCSC_RefGene_Name, data = beta_gene, FUN = median)
    write.table(beta_gene, file.path(outdir, "Gene_level_Beta_values.txt"), row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  # Differential methylation analysis
  if (dm.Anls) {
    group <- factor(targets$Sample_Group)
    design <- model.matrix(~0 + group, data = targets)
    colnames(design) <- levels(group)
    fit <- lmFit(Mval_filtered, design)
    fit2 <- eBayes(fit)
    for (coef in colnames(design)) {
      DMPs <- topTable(fit2, coef = coef, num = Inf)
      write.csv(DMPs, file.path(outdir, paste0("DMPs_", coef, ".csv")), row.names = TRUE)
    }
  }
}

sample<-read.csv("MB_samples.csv",header = T)
for (plate in unique(sample$Platform_methy)) {
  # Select targets for the current plate
  targets <- sample[sample$Platform_methy == plate, ]
  
  # Print number of samples for the current plate
  cat(paste0("Number of samples for ", plate, ": ", nrow(targets), "\n"))
  
  # Select required columns
  targets <- targets[, c("Sentrix_ID", "Sentrix_position", "idat_filename")]
  targets$Basename <- paste0("./", targets$idat_filename)
  
  # Check if all files exist
  if (!all(file.exists(targets$Basename))) {
    cat("Warning: One or more files do not exist in targets\n")
  }
  
  # Create output directory if it does not exist
  outdir <- paste0(getwd(),"/Methylation_Data_processed/", plate)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Run methylation processing and handle potential errors
  tryCatch({
    methyl <- run_methyl(targets = targets, levels_beta = NULL, outdir = outdir)
    cat(paste0(plate, " Processed files saved in ", outdir, "\n"))
  }, error = function(e) {
    cat(paste0("Error in processing plate ", plate, ": ", e$message, "\n"))
  })
}
