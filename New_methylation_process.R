library(minfi)
library(dplyr)
source("advanced_DNA_methylation_process.R")

sample<-read.csv("all_MBGP34_Samples.csv",header = T)
sample<-sample[!duplicated(sample$idat_filename),]
outdir<- paste0(getwd(),"/Methylation_Data_processed_with_Raw")
for (plate in unique(sample$Platform_methy)) {
  # Select targets for the current plate
  targets <- sample[sample$Platform_methy == plate, ]
  
  # Print number of samples for the current plate
  cat(paste0("Number of samples for ", plate, ": ", nrow(targets), "\n"))
  
  # Select required columns
  targets <- targets[, c("Sentrix_ID", "Sentrix_position", "idat_filename")]
  targets$Basename <- paste0("./IDAT_files/", targets$idat_filename)
  
  # Check if all files exist
  if (!all(file.exists(targets$Basename))) {
    cat("Warning: One or more files do not exist in targets\n")
  }
  
  # Create output directory if it does not exist
  plate_outdir <- file.path(outdir, plate)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  # Run methylation processing and handle potential errors
  tryCatch({
    methyl <- run_methyl(targets = targets,norm.method = "preprocessRaw",probe.filter = FALSE,  levels_beta = NULL, outdir = plate_outdir)
    cat(paste0(plate, " Processed files saved in ", outdir, "\n"))
  }, error = function(e) {
    cat(paste0("Error in processing plate ", plate, ": ", e$message, "\n"))
  })
}


## merge all data sets
library(tidyverse)
library(dplyr)
library(purrr)
Beta_files<-list.files(path = outdir,pattern = "Bval.rds", recursive = TRUE,)
Beta_data_list<-list()

for (file in Beta_files) {
  file_name <- sapply(strsplit(file, "HumanMethylation"), `[`, 2)
  file_name <- sapply(strsplit(file_name, "/"), `[`, 1)
  
  data <- data.frame(readRDS(file))
  
  if (file_name=="450"){
    print("skiping row processing..1")
    data$rownames<-rownames(data)
  }else{
    data$rownames <- sapply(strsplit(rownames(data), "_"), `[`, 1)
    data<- data[!duplicated(data$rownames),]
  }

  Beta_data_list[[file_name]] <- data
  
}



# Then use purrr::reduce with a join:
merged_df <- reduce(Beta_data_list, function(x, y) {
  inner_join(x, y, by = "rownames")  # or full_join, left_join, etc.
})

merged_df<-merged_df %>% dplyr::select(rownames, everything())
saveRDS(merged_df,file="NCI_combined_plate_BetaValue.rds")
write.csv(merged_df, file.path(outdir,"NCI_combined_plate_BetaValue.csv"), row.names = F)


