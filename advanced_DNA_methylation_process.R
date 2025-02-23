
library(easypackages)
libraries("minfi","GEOquery","limma","tidyverse","openxlsx","Gviz","DMRcate","stringr","missMethyl","psych", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICanno.ilm10b4.hg19", "IlluminaHumanMethylationEPICv2anno.20a1.hg38")

                   #### 1.1 Detection p-value Filtering ####

#' Filter RGChannelSet by detection p-values
#'
#' @param RGset An RGChannelSet object (from minfi).
#' @param pval.threshold A detection p-value threshold (default 0.01).
#' @param max.sample.fail Maximum fraction of probes allowed to fail per sample
#'        before that sample is removed (default 0.05 => 5%).
#' @param max.probe.fail Maximum fraction of samples allowed to fail per probe
#'        before that probe is removed (default 0.10 => 10%).
#' @return A filtered RGChannelSet.
filter_detection_pvalues <- function(RGset,
                                     pval.threshold = 0.01,
                                     max.sample.fail = 0.05,
                                     max.probe.fail = 0.10) {
  message("Computing detection p-values...")
  detP <- detectionP(RGset)
  
  # 1) Remove poor-quality samples
  sampleFailRate <- colMeans(detP > pval.threshold)
  keepSamples <- sampleFailRate < max.sample.fail
  if (any(!keepSamples)) {
    message("Removing ", sum(!keepSamples), 
            " samples with >", max.sample.fail * 100, 
            "% probes failing detection.")
    RGset <- RGset[, keepSamples]
    detP  <- detP[, keepSamples]
  }
  
  # 2) Remove probes failing in too many remaining samples
  probeFailRate <- rowMeans(detP > pval.threshold)
  keepProbes <- probeFailRate < max.probe.fail
  if (any(!keepProbes)) {
    message("Removing ", sum(!keepProbes), 
            " probes failing in >", max.probe.fail * 100, 
            "% of samples.")
    RGset <- RGset[keepProbes, ]
  }
  
  return(RGset)
}


                  #### 1.2 Cross-reactive, SNP, XY Probe Filtering ####

#' Filter MethylationSet by ambiguity (multi-hit), cross-reactive, SNP, XY, etc.
#'
#' This function removes various categories of "problematic" probes:
#' - Probes listed in an "ambiguity" (multi-hit/cross-reactive) file (often from Chen or Price).
#' - Probes in an EPIC-specific list of known problematic sites (if you have such a file).
#' - Probes known to overlap SNPs (SNP probe list).
#' - Probes on the X or Y chromosome (from a file listing them).
#' - Probes whose names start with "rs" (dbSNP probes) or "ch" (non-CpG controls).
#'
#' @param MSet A Methylation object (e.g. MSet, after normalization).
#' @param ambig.file Path to a file listing "ambiguity" (multi-hit/cross-reactive) probes.
#'                   Typically from Chen et al. or Price et al.
#' @param epic.file Path to an EPIC-specific file of known problematic probes (optional).
#' @param snp.file Path to a file of SNP-containing probes (optional).
#' @param xy.file Path to a file listing probes on chrX or chrY (optional).
#' @param remove.rs Logical; remove probes whose IDs begin with "rs" (SNP probes).
#' @param remove.ch Logical; remove probes whose IDs begin with "ch" (non-CpG or control).
#'
#' @return A filtered Methylation object (the same class as `MSet`), 
#'         with the specified "problematic" probes removed.
#'
#' @examples
#' # Suppose you have these text files of probe IDs (one probe ID per line):
#' # ambig_probes.txt - from Chen's cross-reactive or "multi-hit" list
#' # epic_probes.txt  - known epic-specific problematic probes
#' # snp_probes.txt   - known SNP probes
#' # xy_probes.txt    - probes mapped to chrX or chrY
#' #
#' # Then call:
#' # MSet_filtered <- filter_problematic_probes(
#' #   MSet,
#' #   ambig.file = "ambig_probes.txt",
#' #   epic.file  = "epic_probes.txt",
#' #   snp.file   = "snp_probes.txt",
#' #   xy.file    = "xy_probes.txt"
#' # )
filter_problematic_probes <- function(GSet,
                                      probe.filter.folder=NULL,
                                      ambig.filter = TRUE,
                                      epic.filter  = TRUE,
                                      snp.filter   = TRUE,
                                      xy.filter    = TRUE,
                                      remove.rs  = TRUE,
                                      remove.ch  = TRUE) {
  message("Filtering ambiguous (multi-hit), cross-reactive, SNP, XY, rs*, ch* probes...")

  # Build default file paths
  if (is.null(probe.filter.folder)){
    folder<-"/data/asanigaris/Complete_Deconvolution/MB_Methylation/mnp_training-master/filter"
    
  } else {
    folder<-probe.filter.folder
  }
  
  message("filter probe file are in:", folder)

  ambig.file <- file.path(folder, "amb_3965probes.vh20151030.txt")
  epic.file  <- file.path(folder, "epicV1B2_32260probes.vh20160325.txt")
  snp.file   <- file.path(folder, "snp_7998probes.vh20151030.txt")
  xy.file    <- file.path(folder, "xy_11551probes.vh20151030.txt")

  # Collect all probe indices to remove in this vector
  remove_idxs <- integer(0)
  
  # 1) Ambiguity / cross-reactive file (e.g., Chen multi-hit list)
  if (!is.null(ambig.filter) && file.exists(ambig.file)) {
    ambig_list <- read.table(ambig.file, header = FALSE, stringsAsFactors = FALSE)
    ambig_idx  <- match(ambig_list[, 1], rownames(GSet))
    remove_idxs <- c(remove_idxs, ambig_idx)
  }
  
  # 2) EPIC-specific problematic probes
  if (!is.null(epic.filter) && file.exists(epic.file)) {
    epic_list <- read.table(epic.file, header = FALSE, stringsAsFactors = FALSE)
    epic_idx  <- match(epic_list[, 1], rownames(GSet))
    remove_idxs <- c(remove_idxs, epic_idx)
  }
  
  # 3) SNP-overlapping probes
  if (!is.null(snp.filter) && file.exists(snp.file)) {
    snp_list <- read.table(snp.file, header = FALSE, stringsAsFactors = FALSE)
    snp_idx  <- match(snp_list[, 1], rownames(GSet))
    remove_idxs <- c(remove_idxs, snp_idx)
  }
  
  # 4) XY-chromosome probes
  if (!is.null(xy.filter) && file.exists(xy.file)) {
    xy_list <- read.table(xy.file, header = FALSE, stringsAsFactors = FALSE)
    xy_idx  <- match(xy_list[, 1], rownames(GSet))
    remove_idxs <- c(remove_idxs, xy_idx)
  }
  
  # 5) Probes starting with "rs" or "ch"
  #    "rs" = dbSNP, "ch" = non-CpG control or other special cases
  if (remove.rs) {
    rs_idx <- grep("^rs", rownames(MSet))
    remove_idxs <- c(remove_idxs, rs_idx)
  }
  if (remove.ch) {
    ch_idx <- grep("^ch", rownames(MSet))
    remove_idxs <- c(remove_idxs, ch_idx)
  }
  
  # Clean up indices (remove NA, duplicates)
  remove_idxs <- unique(remove_idxs[!is.na(remove_idxs)])
  
  # Filter out the unwanted probes
  if (length(remove_idxs) > 0) {
    message("Removing ", length(remove_idxs), 
            " probes (ambiguity/cross-reactive/SNP/XY/rs/ch).")
    Gset_filtered  <- GSet[!remove_idxs, ]
  } else {
    message("No matching probes found in GSet based on the provided lists.")
  }
  
  return(Gset_filtered)
}


               ####   1.3 Gene-Level Beta Calculation (Optional Helper)   ####

#' Calculate gene-level beta values
#'
#' @param Bval A beta value matrix (rows = probes, cols = samples).
#' @param anno An annotation data.frame (with rownames = probe IDs).
#' @param levels_beta A character vector specifying annotation column and 
#'        grouping. Example: c("UCSC_RefGene_Name", "TSS1500;TSS1500").
#' @param out.path File path to save the resulting gene-level beta table.
#' @return A data.frame with aggregated gene-level beta values.
calculate_gene_level_beta <- function(Bval, 
                                      anno, 
                                      levels_beta = c("UCSC_RefGene_Name", "TSS1500;TSS1500"),
                                      out.path = NULL) {
  
  if (length(levels_beta) < 2) {
    stop("levels_beta should contain [annotation column, group info], e.g., 
         c('UCSC_RefGene_Name','TSS1500;TSS1500').")
  }
  
  group_col <- levels_beta[1]
  group_spec <- levels_beta[2]
  
  # Copy annotation so we can subset
  annoG <- anno
  
  # The group_spec might contain semicolon-delimited categories
  if (length(strsplit(group_spec, ";")[[1]]) >= 2) {
    # e.g. TSS1500;TSS200
    grouping_list <- c(group_spec, sapply(strsplit(group_spec, ";"), `[`, 1))
    annoG <- annoG[annoG$UCSC_RefGene_Group %in% grouping_list, ]
  } else {
    annoG <- annoG[annoG$UCSC_RefGene_Group == group_spec, ]
  }
  
  # Keep relevant columns
  annoG <- annoG[, c("Name", group_col)]
  annoG <- annoG[annoG[[group_col]] != "", ]
  annoG[[group_col]] <- sapply(strsplit(annoG[[group_col]], ";"), `[`, 1)
  
  rownames(annoG) <- annoG$Name
  annoG <- annoG[, group_col, drop = FALSE]
  
  # Merge with beta values
  commonProbes <- intersect(rownames(annoG), rownames(Bval))
  beta_gene <- merge(annoG[commonProbes, , drop=FALSE],
                     Bval[commonProbes, , drop=FALSE],
                     by = 0,
                     all = FALSE)
  
  # Drop rownames col
  beta_gene <- beta_gene[, -1]
  
  # Summarize by median (or mean if you prefer)
  beta_gene <- aggregate(. ~ get(group_col), data = beta_gene, FUN = median)
  colnames(beta_gene)[1] <- group_col  # rename the aggregated column
  
  # Save if desired
  if (!is.null(out.path)) {
    write.table(beta_gene, out.path, 
                row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  return(beta_gene)
}



                ####  2. Main Pipeline Function  ####

run_methyl <- function(targets, 
                       norm.method = "preprocessFunnorm", 
                        # -- Quality contol arguments --
                       Quality.control = FALSE,
                       pval.threshold = 0.01,
                       max.sample.fail = 0.05,
                       max.probe.fail = 0.10,
                       # -- Filter arguments --
                       probe.filter = TRUE, 
                       probe.filter.folder = NULL,
                       ambig.filter = TRUE,
                       epic.filter  = TRUE,
                       snp.filter   = TRUE,
                       xy.filter    = TRUE,
                       # -- other options --
                       levels_beta = c("UCSC_RefGene_Name", "TSS1500;TSS1500"), 
                       outdir = "MethylResult", 
                       dm.Anls = FALSE) {
  # Create output directory
  if (!file.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }
  
  # 1) Read IDATs -> RGChannelSet
  RGset <- read.metharray.exp(targets = targets, verbose = TRUE, force = TRUE)
  cat("Reading IDAT files is completed!\n")
  
  # 2) (Optional) Detection p-value filtering on the RGset
  if (Quality.control) {
    RGset <- filter_detection_pvalues(RGset,
                                      pval.threshold = pval.threshold,
                                      max.sample.fail = max.sample.fail,
                                      max.probe.fail = max.probe.fail)
  }
  
  # Save these intermediate objects
  saveRDS(RGset, file = file.path(outdir, "RGSet.rds"))

  # 3) Preprocessing / normalization
  if (norm.method == "preprocessRaw") {
    MSet <- preprocessRaw(RGset)
  } else if (norm.method == "preprocessIllumina") {
    MSet <- preprocessIllumina(RGset)
  } else if (norm.method == "preprocessSWAN") {
    MSet <- preprocessSWAN(RGset)
  } else if (norm.method == "preprocessNoob") {
    MSet <- preprocessNoob(RGset)
  } else if (norm.method == "preprocessFunnorm") {
    GSet <- preprocessFunnorm(RGset)
  } else if (norm.method == "preprocessQuantile") {
    GSet <- preprocessQuantile(RGset)
  }
  cat("Normalization process is completed!\n")

  
  # 4) Convert to RatioSet, then map to genome => Gset
  #    SKIP if norm.method == "preprocessFunnorm" OR "preprocessQuantile"

  if (norm.method %in% c("preprocessFunnorm","preprocessQuantile")) {
    cat("Skipping ratioConvert() and mapToGenome() for Funnorm/Quantile data.\n")
    # In fact, preprocessFunnorm() and preprocessQuantile() 
    # already return a GenomicRatioSet, so just rename:
  
  } else {
    cat("Converting to RatioSet and mapping to genome...\n")
    RatioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
    Gset <- mapToGenome(RatioSet)
  }
  
  # 5) Further filtering for cross-reactive, XY, etc.
  if (probe.filter) {
    GSet <- filter_problematic_probes(GSet,
                                      probe.filter.folder = probe.filter.folder,
                                      ambig.filter = ambig.filter ,
                                      epic.filter  = epic.filter,
                                      snp.filter   = snp.filter,
                                      xy.filter   = xy.filter,
                                      remove.rs  = TRUE,
                                      remove.ch  = TRUE)
  }
  
  
  saveRDS(Gset, file = file.path(outdir, "Gset.rds"))
  
  # 6) Get final Beta / M-values
  Bval <- getBeta(Gset)
  Mval <- getM(Gset)
  CN <- getCN(Gset)
  saveRDS(Bval, file = file.path(outdir, "Bval.rds"))
  saveRDS(Mval, file = file.path(outdir, "Mval.rds"))
  cat("Beta and M-value conversion process is completed!\n")
  
  # 7) Annotation (for 450K or EPIC)
  if (annotation(MSet)[["array"]] %in% c("IlluminaHumanMethylation450k")) {
    anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  } else {
    # Assume EPIC
    anno <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  }
  
  # 8) Gene-level beta calculation
  if (!rlang::is_empty(levels_beta)) {
    gene_beta <- calculate_gene_level_beta(
      Bval, 
      anno,
      levels_beta = levels_beta,
      out.path = file.path(outdir, "Gene_level_Beta_values.txt")
    )
  }
  
  # 9) (Optional) Differential Methylation Analysis
  if (dm.Anls) {
    # Example: group is in targets$Sample_Group
    group <- factor(targets$Sample_Group)
    design <- model.matrix(~0 + group, data = targets)
    colnames(design) <- levels(group)
    
    # Fit linear model using the final, filtered M-values
    fit <- limma::lmFit(Mval, design)
    fit2 <- limma::eBayes(fit)
    
    # Loop over contrasts in design
    for (coef in colnames(design)) {
      DMPs <- limma::topTable(fit2, coef = coef, number = Inf)
      write.csv(DMPs, file.path(outdir, paste0("DMPs_", coef, ".csv")), row.names = TRUE)
    }
    cat("Differential methylation analysis (DMP) is completed!\n")
  }
  
  cat("run_methyl pipeline completed.\n")
}


