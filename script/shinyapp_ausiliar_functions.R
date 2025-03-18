

validate_input <- function(file, required_cols) {
  data <- fread(file$datapath)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Input file is missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  return(data)
}

load_data <- function(input_file, default_path, required_cols = NULL) {
  if (!is.null(input_file)) {
    cat("Loading user-provided file...\n")
    data <- fread(input_file$datapath)
  } else {
    cat("Loading default file from", default_path, "...\n")
    data <- fread(file.path("www", default_path))
  }
  
  if (!is.null(required_cols)) {
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      stop(paste("File is missing required columns:", paste(missing_cols, collapse = ", ")))
    }
  }
  
  return(data)
}
