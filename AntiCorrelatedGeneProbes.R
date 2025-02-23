
library(minfi)
library(psych)
library(tidyverse)
library(dplyr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)


# 1) Correlation between methylated probes and Gene expression values 

# Outside the function, please load required packages and annotation:
# library(minfi)
# library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# library(dplyr)

probeGeneCor <- function(
  betaData,
  expData,
  array.type       = "450K",
  RefGene_Group    = "TSS1500",   # e.g. "TSS1500", "TSS200", "5'UTR"
  match_mode       = c("exact", "multi", "partial"),
  genelist         = NULL,
  corr.method      = "pearson"
) {
  # 1) Load the correct annotation
  if (array.type == "450K") {
    anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } else if (array.type == "EPIC") {
    anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  } else {
    stop("Unsupported array type. Choose '450K' or 'EPIC'.")
  }
  
  annoG <- as.data.frame(anno)
  
  # 2) Determine how to filter 'UCSC_RefGene_Group' based on match_mode
  match_mode <- match.arg(match_mode)  # ensures valid choice
  
  if (match_mode == "exact") {
    # Keep rows that are exactly the 'RefGene_Group'
    keep_idx <- (annoG$UCSC_RefGene_Group == RefGene_Group)
    
  } else if (match_mode == "multi") {
    # Pattern: ^RefGene_Group(;RefGene_Group)*$
    # e.g. if RefGene_Group="TSS1500",
    # then "TSS1500;TSS1500;TSS1500" matches,
    # but "TSS1500;TSS200" does NOT.
    
    # Escape any regex special chars in RefGene_Group, just in case
    refVal_escaped <- gsub("([\\.\\^\\$\\|\\(\\)\\[\\]\\{\\}\\*\\+\\?\\\\])", "\\\\\\1", RefGene_Group)
    regex_pattern  <- paste0("^", refVal_escaped, "(;", refVal_escaped, ")*$")
    
    keep_idx <- grepl(regex_pattern, annoG$UCSC_RefGene_Group)
    
  } else if (match_mode == "partial") {
    # match_mode == "partial"
    # Keep rows that contain the RefGene_Groupe anywhere
    # e.g. "TSS1500;TSS200", "TSS200;5'UTR"
    
    refVal_escaped <- gsub("([\\.\\^\\$\\|\\(\\)\\[\\]\\{\\}\\*\\+\\?\\\\])", "\\\\\\1", RefGene_Group)
    
    keep_idx <- grepl(refVal_escaped, annoG$UCSC_RefGene_Group)
  }
  
  annoG <- annoG[keep_idx, ]
  
  # 3) If a row has multiple gene names (e.g. "BRCA1;TP53"), keep only the first
  annoG$UCSC_RefGene_Name <- sapply(strsplit(annoG$UCSC_RefGene_Name, ";"), `[`, 1)
  
  # 4) If a genelist is provided, keep only those genes
  if (!is.null(genelist)) {
    annoG <- annoG[annoG$UCSC_RefGene_Name %in% genelist, ]
  }
  
  # 5) Keep only CpGs that exist in 'betaData' row names
  #    'Name' is the Illumina probe ID
  annoG <- annoG[annoG$Name %in% rownames(betaData), c("Name", "UCSC_RefGene_Name")]
  colnames(annoG) <- c("cpg", "gene")
  
  # 6) Check which genes exist in expression data
  commonGenes <- intersect(unique(annoG$gene), rownames(expData))
  if (length(commonGenes) == 0) {
    message("No overlap between filtered annotation genes and expData. Returning empty result.")
    return(data.frame(gene=character(), cpg=character(), pval=numeric(), cor=numeric()))
  }
  
  # 7) Identify common samples
  commonSamples <- intersect(colnames(betaData), colnames(expData))
  if (length(commonSamples) < length(commonSamples)/2) {
    stop("Fewer than 50% common samples between betaData and expData. Cannot compute correlation.")
  }
  
  # 8) Loop over gene-CpG pairs to calculate correlation
  result_list <- vector("list", length=100000)
  k <- 0
  
  for (g in commonGenes) {
    # cpgs mapped to this gene
    cpg_idx <- which(annoG$gene == g)
    if (length(cpg_idx) == 0) next
    
    # expression vector
    exp_vec <- as.numeric(expData[g, commonSamples])
    
    for (i in seq_along(cpg_idx)) {
      cpg_id  <- annoG$cpg[cpg_idx[i]]
      beta_vec <- as.numeric(betaData[cpg_id, commonSamples])
      
      # correlation
      ct <- cor.test(exp_vec, beta_vec, method = corr.method)
      
      k <- k + 1
      result_list[[k]] <- data.frame(
        gene = g,
        cpg  = cpg_id,
        pval = ct$p.value,
        cor  = unname(ct$estimate)
      )
    }
  }
  
  result_list <- result_list[seq_len(k)]  # trim
  results_df <- do.call(rbind, result_list)
  rownames(results_df) <- NULL
  
  cat("Correlation is completed!\n",
      "RefGene filter mode:", match_mode,
      "| RefGene value:", RefGene_Group, "\n",
      "Total pairs:", nrow(results_df), "\n")
  
  return(results_df)
}



# 2) Function: Identify Anti-Correlated Probes

#' Identify anti-correlated probes and save CSV files
#'
#' This function reads:
#'   1) a methylation RDS file (cbtn_embl_beta),
#'   2) an RDS file with cell-type expression data,
#'   3) CSV files for Malignant / Myeloid DEGs,
#' and calculates probe–gene correlations using `probeGeneCor()`.
#' It filters for logFC > 0.1, pval < 0.05, and correlation < -0.1,
#' then saves anti-correlated probes to CSV for each cell type + gene cutoff.
#'
#' @param beta_file Path to the RDS file with beta matrix (e.g. "cbtn_embl_betaValues.rds").
#' @param cell_data_file      Path to the RDS file with expression data (named list by cell type).
#' @param deg_file_list      list of CSV with DEGs (e.g. "MalignantClusterNamedProteinCodingDEGs.csv").
#' @param GeneNum             Integer vector of gene cutoffs to test (e.g. c(100,300,500,...)).
#' @param Cells               Character vector of cell types (e.g. c("Malignant","Myeloid")).
#' @param cor_cutoff          Numeric, correlation threshold for "anti-correlation" (default -0.1).
#'
#' @return (Invisibly) a list of data frames (the last computed anti-correlated result per cell type),
#'         but primarily writes CSV files to disk.
#'
#' @examples
#' identify_anticorr_probes(
#'   beta_file = "cbtn_embl_betaValues.rds",
#'   cell_data_file      = "allDattaSixCelltypeExp.rds",
#'   deg_file_list      = "list(malig,myeloid)",
#'   GeneNum             = c(100,300,500),
#'   Cells               = c("Malignant","Myeloid"),
#'   cor_cutoff          = -0.1
#')
identify_anticorr_probes <- function(
  beta_file,
  cell_data_file,
  deg_files_list,
  GeneNum    = c(500, 1000, 2000, 3000, 5000),
  Cells      = c("Malignant", "Myeloid"),
  cor_cutoff = -0.1,
  outdir     = "AntiCorrelatedGeneProbes"
) {

  # Create output directory
  if (!file.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }
  #----------------------------------------------------------------
  # 1) LOAD BETA DATA (either from RDS file path or in-memory)
  #----------------------------------------------------------------
  if (is.character(beta_file)) {
    # User passed a path
    message("Loading beta data from file: ", beta_file)
    beta_df <- readRDS(beta_file)
  } else {
    # Assume it's already an object in memory
    message("Using provided beta data object directly.")
    beta_df <- beta_file
  }
  
  #----------------------------------------------------------------
  # 2) LOAD CELL-TYPE EXPRESSION DATA
  #----------------------------------------------------------------
  if (is.character(cell_data_file)) {
    message("Loading cell-type expression data from file: ", cell_data_file)
    cell_data <- readRDS(cell_data_file)
  } else {
    message("Using provided cell_data object directly.")
    cell_data <- cell_data_file
  }
  
  #----------------------------------------------------------------
  # 3) LOAD DEGs LIST
  #----------------------------------------------------------------
  # 'deg_files_list' can also be a path to an RDS or a list object
  if (is.character(deg_files_list)) {
    message("Loading DEG list from file: ", deg_files_list)
    deg_list <- readRDS(deg_files_list)
  } else {
    message("Using provided deg_files_list object directly.")
    deg_list <- deg_files_list
  }
  
  #----------------------------------------------------------------
  # 4) ITERATE OVER CELL TYPES AND GENE NUM CUTOFFS
  #----------------------------------------------------------------
  results_list <- list()  # to store results if needed
  
  for (Cell in Cells) {
    
    # Check if the DEGs for this cell exist in deg_list
    if (! Cell %in% names(deg_list)) {
      warning("Cell type '", Cell, "' not found in deg_list. Skipping...")
      next
    }
    
    # 4A: Subset the correct DEG data frame
    deg_df <- deg_list[[Cell]]
    
    for (Num in GeneNum) {
      
      # 4B: Filter top genes (logFC>0.1, pval<0.05, top 'Num' by logFC per group)
      topGenes <- deg_df %>%
        dplyr::filter(logFC > 0.1 & pval < 0.05) %>%
        dplyr::group_by(group) %>%
        dplyr::arrange(desc(logFC), .by_group = TRUE) %>%
        dplyr::slice_max(order_by = logFC, n = Num)
      
      # 4C: Get matching expression data for this cell type
      if (! Cell %in% names(cell_data)) {
        warning("Cell type '", Cell, "' not found in cell_data. Skipping...")
        next
      }
      exp <- cell_data[[Cell]]
      
      # 4D: Match columns in beta_df
      common_samples <- intersect(colnames(exp), colnames(beta_df))
      exp  <- exp[, common_samples, drop = FALSE]
      beta <- beta_df[, common_samples, drop = FALSE]
      
      # 4E: Calculate correlation (using your custom probeGeneCor function)
      allprobCorr_anti <- probeGeneCor(
        betaData      = beta,
        expData       = exp,
        genelist      = topGenes$feature,
        RefGene_Group = "TSS1500",
        match_mode    = "exact",
        array.type    = "450K"   # or "EPIC" if needed
      )
      
      # 4F: Filter out NA, keep cor < cor_cutoff
      allprobCorr_anti <- stats::na.omit(allprobCorr_anti)
      allprobCorr_anti <- allprobCorr_anti[allprobCorr_anti$cor < cor_cutoff, ]
      
      # 4G: Merge with topGenes
      merged <- merge(topGenes, allprobCorr_anti, 
                      by.x = "feature", by.y = "gene", all = FALSE)
      merged <- merged[!duplicated(merged$cpg), ]
      
      # 4H: Write results to CSV
      #     Make a subdir for each cell
      cell_outdir <- file.path(outdir, Cell)
      dir.create(cell_outdir, showWarnings = FALSE, recursive = TRUE)
      
      out_file <- file.path(cell_outdir,
                            paste0(Cell, "Cellstate", Num, 
                                   "Gene_AnticorProbesNeg.csv"))
      message("Writing file: ", out_file)
      write.csv(merged, out_file, row.names = FALSE)
      
      # Store result in list if needed
      results_list[[paste(Cell, Num, sep = "_")]] <- merged
    }
  }
  
  #----------------------------------------------------------------
  message("Finished identify_anticorr_probes()!")
  saveRDS(results_list, file=file.path(outdir,"AntiCorGeneprobeList.rds"))
  invisible(results_list)  # return (invisibly) the list of results
}



# 2) Function: Generate Gene-Level Beta Values
#' Generate gene-level beta values from anti-correlated probes
#'
#' This function:
#'   - Loops over the same `Cells` and `GeneNum` as in `identify_anticorr_probes()`.
#'   - Reads the previously saved CSV files for anti-correlated probes (one per cell type + gene cutoff).
#'   - Filters those probes by stricter correlation thresholds (`cor_values`).
#'   - Merges with the `NCI_beta_new` data, aggregates (median) by "feature" (gene),
#'   - Saves final CSV files in `NCI_Methylation/NCI_Genelevel_beta/<Cell>_CBTN_EMBL_450K_TSS1500/`.
#'
#' @param NCI_beta_file Path to the RDS containing a beta matrix object named "NCI_beta_new" in your original code.
#' @param Cells         Character vector (e.g., c("Malignant","Myeloid")), matching what was used previously.
#' @param GeneNum       Same integer vector used previously for gene cutoffs.
#' @param cor_values    Vector of correlation thresholds to filter further (e.g. c(-0.1, -0.2, -0.3)).
#'
#' @return No return value; writes final gene-level beta CSV files to disk.
#'
#' @examples
#' generate_gene_level_beta(
#'   NCI_beta_file = "NCI_beta_new.rds",
#'   Cells         = c("Malignant","Myeloid"),
#'   GeneNum       = c(100, 300, 500),
#'   cor_values    = c(-0.1, -0.2, -0.3, -0.4)
#' )
generate_gene_level_beta <- function(
  beta_file,
  Cells         = c("Malignant","Myeloid"),
  GeneNum       = c(100,300,500),
  cor_values    = c(-0.1, -0.2, -0.3, -0.4),
  outdir        = "Genelevel_beta",
  antiCorResults = NULL  # Optional list from identify_anticorr_probes()
) {
  # 0. Create output directory (top-level)
  if (!file.exists(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
  }

  # 1. Load the final beta data for merging (a matrix/data.frame with rownames = cpg IDs)
  if (is.character(beta_file)) {
    message("Loading beta data from file: ", beta_file)
    beta_df <- readRDS(beta_file)
    beta_df<-data.frame(beta_df)
  } else {
    message("Using provided beta data object directly.")
    beta_df <- beta_file
    beta_df<-data.frame(beta_df)
  }
  
  # 2. Loop over each cell type, gene number, and correlation threshold
  for (Cell in Cells) {
    
    # Create a subdirectory for each cell inside outdir
    cell_path <- file.path(outdir, Cell)
    dir.create(cell_path, recursive = TRUE, showWarnings = FALSE)
    
    for (n_Gene in GeneNum) {
      
      #-------------------------------------------------------
      # 2A. Get the anti-cor probe data (either from memory list or CSV)
      #-------------------------------------------------------
      antiProbe_df <- NULL
      
      if (is.list(antiCorResults)) {
        # If a list was provided, we try to retrieve the data frame directly
        list_key <- paste(Cell, n_Gene, sep = "_")
        if (!list_key %in% names(antiCorResults)) {
          message("No result for (Cell=", Cell, ", nGene=", n_Gene, ") in antiCorResults. Skipping...")
          next
        }
        antiProbe_df <- antiCorResults[[list_key]]
        
      } else if (is.character(antiCorResults)) {
        # Otherwise, read CSV from folder or directory
        anti_file <- file.path(antiCorResults,Cell, 
          paste0(Cell, "Cellstate", n_Gene, "Gene_AnticorProbesNeg.csv")
        )
        
        if (!file.exists(anti_file)) {
          message("File not found: ", anti_file, " (Skipping...)")
          next
        }
        
        antiProbe_df <- read.csv(anti_file, header = TRUE, stringsAsFactors = FALSE)
      } else {
        stop("antiCorResults must be either a directory path (character) or a list.")
      }
      
      #-------------------------------------------------------
      # 2B. For each correlation threshold, filter further
      #-------------------------------------------------------
      for (cor_val in cor_values) {
        
        # Keep only rows with cor < cor_val
        subDF <- antiProbe_df[antiProbe_df$cor < cor_val, ]
        if (nrow(subDF) == 0) {
          message("No probes pass cor<", cor_val, " for ", Cell, ", Gene=", n_Gene)
          next
        }
        
        # Keep only cpg + feature columns
        if (!all(c("cpg","feature") %in% colnames(subDF))) {
          warning("The data frame doesn't have 'cpg' or 'feature' columns. Skipping.")
          next
        }
        subDF <- subDF[, c("cpg","feature")]
        
        #-------------------------------------------------------
        # 2C. Merge with beta_df by cpg
        #     (beta_df should have rownames=cpg IDs)
        #-------------------------------------------------------
        GeneBval <- merge(subDF, beta_df, 
                          by.x = "cpg", by.y = 1,  
                          all = FALSE)
        
        # Remove the 'cpg' column
        GeneBval <- GeneBval[, -1, drop=FALSE]
        message("Merged data rows: ", nrow(GeneBval))
        #-------------------------------------------------------
        # 2D. Aggregate by gene (the 'feature' column)
        #     using median across all sample columns
        #-------------------------------------------------------
        GeneBval <- aggregate(. ~ feature, data = GeneBval, FUN = median)
        
        #-------------------------------------------------------
        # 2E. Write to CSV
        #-------------------------------------------------------
        cor_formatted <- gsub("\\.", "", formatC(abs(cor_val), format = "f", digits = 2))
        
        out_file <- file.path(cell_path, 
          paste0(Cell, "Cellstate", n_Gene,
                 "Gene_TSS1500_Median_betaValues_", cor_formatted, ".csv")
        )
        
        write.csv(GeneBval, out_file, row.names = FALSE)
        
        message("Gene-level beta for Cell=", Cell, ", Genes=", n_Gene,
                ", cor<", cor_val, " => written to:\n  ", out_file)
      } # End of cor_values loop
    } # End of GeneNum loop
  } # End of Cells loop
  
  message("Finished generate_gene_level_beta_v2()!")
}


