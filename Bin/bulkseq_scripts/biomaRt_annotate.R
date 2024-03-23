# This script contains the function to acquire annotation infomation from the ensembl database. 
# If this information has already been saved to file, then it reads in said annotation object.

generate_ensembl <- function(mainDir)
{
  # required packages
  annotate_packages <- c("biomaRt",
                         "dplyr")
  lapply(annotate_packages, require, character.only = TRUE)
  
  if (file.exists(sprintf("%s/Annotation/ensembl_annotation.rds", mainDir))) {
    
    print("The 'ensembl_annotation.rds' file already exists, importing file.")
    # Read in the annotation file
    ensembl_annotation <- readRDS(sprintf("%s/Annotation/ensembl_annotation.rds", mainDir))

  } else {
    
    print("Acquiring annotation data from the ensembl database. Please ensure that you have a stable internet connection.")
    ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
    
    setwd(mainDir)
    output_dir <- "Annotation"
    if (!dir.exists(output_dir)){
      dir.create(output_dir)
      print("Created 'Work/Annotation' directory.")
    } else {}
    
    print("writing 'ensembl_annotation.rds' to file")
    saveRDS(object = ensembl, 
            file = sprintf("%s/Annotation/ensembl_annotation.rds", mainDir))
  }
  
}

# The 'convert_ID' function matches two different ID types. 
# Please note: in cases where IDs do not match perfectly, only the first match to the final ID will be retained.
convert_ID <- function(initial_ID, # This is the starting ID that you wish to change to another type.
                       final_ID, # The ID type you want in your final object.
                       object # character vector containing the 'initial ID' type that you want to convert
                       )
  {
  # required packages
  annotate_packages <- c("biomaRt",
                         "dplyr")
  lapply(annotate_packages, require, character.only = TRUE)
  
  # Read in the annotation file
  ensembl_annotation <- readRDS(sprintf("%s/Annotation/ensembl_annotation.rds", mainDir))
  
  annotation_info <- getBM(attributes = c(initial_ID, final_ID),
                           mart = ensembl_annotation) # Produces a data frame matching the two ID types specified by the user.
  colnames(annotation_info) <- c(initial_ID, final_ID)
  
  # Certain ID types are merely numbers. These should be converted to characters.
  if (!is.character(annotation_info[,2])) {
    annotation_info[,2] <- as.character(annotation_info[,2])
  }
  
  object <- as.data.frame(object)
  names(object) <- initial_ID
  
  # The gene names must be matched to the gene IDs of the genes that 
  # are being expressed (ie the gene IDs from the count matrix)
  object_final <- left_join(object, annotation_info)
  
  # Keep only the first match for each initial ID
  object_final <- object_final %>% group_by(across(all_of(initial_ID))) %>% dplyr::slice(1)
  
  # Ensure uniqueness of the resulting data frame
  object_final <- unique(object_final)
  
  # double check that the initial and final IDs have the same number of rows.
  if (identical(nrow(object), nrow(object_final)) == TRUE){
    print(sprintf("%s and %s have been successfully matched.", initial_ID, final_ID))
  } else {
    print(sprintf("Warning: There are an unequal number of %s and %s. Annotation results will be unreliable.", initial_ID, final_ID))
  }
  
  # There may be 2 types of missing data within the 'object_final', 'NA' or '' (an empty cell).
  # Handle missing or empty cells
  missing_cells <- which(is.na(object_final[, final_ID]) | object_final[, final_ID] == "")
  
  # Report the number of 'NA' or empty cells in 'object_final'.
  missing_ID2 <- length(missing_cells)
  if (missing_ID2 == 0) {
    print("No missing gene names were detected.")
  } else {
    print(sprintf("%s target IDs had missing information. These fields will retain their original ID type.", missing_ID2))
  }
  
  # Replace missing and na values with the original ID
  object_final[missing_cells, final_ID] <- object_final[missing_cells, initial_ID]
  
  # Ensure that the 'object_final' is written to the global environment
  assign("object_final", object_final, envir = globalenv())
  
  print(sprintf("Final matchings of %s and %s are complete and can be located in 'object_final'.", initial_ID, final_ID))
  print("Conversion complete")
  print("------------------------------")
}

# The 'gather_attributes' function takes in a vector of a specified ID type and gathers information on specified attributes.
gather_attributes <- function(initial_ID, # This is the starting ID type contained in the 'object'.
                              target_attributes, # A vector listing all the target attributes you wish to gather information for.
                              object # character vector containing the 'initial ID' type that you want to convert
)
{
  # required packages
  annotate_packages <- c("biomaRt",
                         "dplyr")
  lapply(annotate_packages, require, character.only = TRUE)
  
  # Read in the annotation file
  ensembl_annotation <- readRDS(sprintf("%s/Annotation/ensembl_annotation.rds", mainDir))
  
  print("Gathering information on all the specified attributes.")
  print("Please note: this may take some time depending on the type and number of specified attributes.")
  annotation_info <- getBM(attributes = c(initial_ID, target_attributes),
                           filters = initial_ID,
                           values = object,
                           mart = ensembl_annotation) # Produces a data frame matching the two ID types specified by the user.
  colnames(annotation_info) <- c(initial_ID, target_attributes)
  
  object <- as.data.frame(object)
  names(object) <- initial_ID
  
  # The gene names must be matched to the gene IDs of the genes that 
  # are being expressed (ie the gene IDs from the count matrix)
  object_final <- left_join(object, annotation_info)
  
  # Ensure uniqueness of the resulting data frame
  object_final <- unique(object_final)
  
  # Ensure that the 'object_final' is written to the global environment
  assign("object_final", object_final, envir = globalenv())
  
  print("Annotation complete")
  print("------------------------------")
}