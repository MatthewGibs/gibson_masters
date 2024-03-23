# This script generates the metadata object required by DESeq2, and saves it to file.
# NOTE: Should you be using a different dataset than the default for this pipeline then you MUST change these parameters appropriately.

generate_metadata <- function(mainDir) {
  # Create a data frame containing metadata for the txi object. This is essential for DESeq2
  donor_number <- c("One", "Two", "Three", "One", "Two", "Three", "One", "Two", "Three")
  condition <- c("untreated", "untreated", "untreated", "treated", "treated", "treated", "treated", "treated", "treated")
  treatment <- c("none", "none", "none", "MCSF", "MCSF", "MCSF", "GMCSF", "GMCSF", "GMCSF")
  
  # Join the variables to create a data frame
  sample_metadata <- data.frame(donor_number, condition, treatment, stringsAsFactors = TRUE)
  sample_metadata$donor_number <- factor(sample_metadata$donor_number, levels = c("One", "Two", "Three"))
  rownames(sample_metadata) <- sample_names
  
  # Whichever column is to be used for the DEA must have its factor levels reordered.
  # specify the 'ref' as your control condition (in this case the samples that have received no treatment).
  sample_metadata$treatment <- relevel(sample_metadata$treatment, ref = "none")
  
  setwd(mainDir)
  
  output_dir <- "Annotation"
  if (!dir.exists(output_dir)){
    dir.create(output_dir)
    print("Created 'Work/Annotation' directory.")
  } else {}
  
  saveRDS(object = sample_metadata,
          file = sprintf("%s/Annotation/sample_metadata.rds", mainDir))
  print(sprintf("The metadata file has been generated. It can be found in the '%s/Annotation' directory", mainDir))
}

