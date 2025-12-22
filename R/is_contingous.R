  is_contiguous <- function(x) {
    
    # checks whether values of x appear contiguously without being broken by other values or NA's.
     
    r <- rle(x)
    vals <- r$values
    
    # Indices of runs that are not NA
    idx_non_na <- which(!is.na(vals))
    
    # For each non-NA run, record its value and its run index
    df <- data.frame(
      value = vals[idx_non_na],
      run   = idx_non_na
    )
    
    # For each distinct non-NA value, count how many *separate runs* it has
    runs_per_value <- tapply(df$run, df$value, length)
    
    # Every non-NA value must appear in exactly one run
    all(runs_per_value == 1)
  }