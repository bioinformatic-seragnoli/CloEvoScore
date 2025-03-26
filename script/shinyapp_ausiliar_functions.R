# === Purpose: Utility functions to validate and load input files in a Shiny app ===

# === Function 1: Validate uploaded input file structure ===
validate_input <- function(file, required_cols) {
  
  # Read the uploaded file
  data <- fread(file$datapath)
  
  # Check if all required columns are present in the uploaded data
  missing_cols <- setdiff(required_cols, names(data))
  
  # If any required columns are missing, stop and return an error message
  if (length(missing_cols) > 0) {
    stop(paste("Input file is missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # If all checks pass, return the loaded data
  return(data)
}

# === Function 2: Load user-uploaded file or fall back to default ===
load_data <- function(input_file, default_path, required_cols = NULL) {
  
  # If a user-provided file is uploaded
  if (!is.null(input_file)) {
    cat("Loading user-provided file...\n")
    data <- fread(input_file$datapath)
  } else {
    
    # If no file uploaded, load the default file from the www/ folder
    cat("Loading default file from", default_path, "...\n")
    data <- fread(file.path("www", default_path))
  }
  
  # If required column names are specified, check for their presence
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(data))
    
    # Stop execution and report missing columns if any are absent
    if (length(missing_cols) > 0) {
      stop(paste("File is missing required columns:", paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Return the validated and loaded dataset
  return(data)
}
