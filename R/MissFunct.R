# ============================================================
# BIOFLOW - PHENOTYPE QA/QC
# Missingness analysis functions
# ============================================================


# ------------------------------------------------------------
# 1. Design variable missingness
# ------------------------------------------------------------

check_design_missingness <- function(data,
                                     design_vars = c(
                                       "genotype",
                                       "designation",
                                       "environment",
                                       "management",
                                       "replicate",
                                       "block",
                                       "row",
                                       "column",
                                       "plot"
                                     )) {
  
  vars_present <- intersect(design_vars, names(data))
  
  if (length(vars_present) == 0) {
    return(data.frame())
  }
  
  result <- lapply(vars_present, function(v) {
    
    n_missing <- sum(is.na(data[[v]]) | data[[v]] == "", na.rm = TRUE)
    
    data.frame(
      variable = v,
      n_total = nrow(data),
      n_missing = n_missing,
      pct_missing = round(100 * n_missing / nrow(data),4),
      stringsAsFactors = FALSE
    )
  })
  
  result <- do.call(rbind, result)
  
  result$status <- ifelse(
    result$pct_missing == 0,
    "OK",
    ifelse(
      result$pct_missing < 10,
      "WARNING",
      "HIGH"
    )
  )
  
  result
}


# ------------------------------------------------------------
# 2. Identify rows with missing design information
# ------------------------------------------------------------

get_design_missing_rows <- function(data,
                                    design_vars) {
  
  vars_present <- intersect(design_vars, names(data))
  
  if (length(vars_present) == 0) {
    return(data.frame())
  }
  
  missing_matrix <- sapply(
    data[vars_present],
    function(x) is.na(x) | x == ""
  )
  
  if (is.null(dim(missing_matrix))) {
    missing_matrix <- matrix(
      missing_matrix,
      ncol = 1
    )
  }
  
  affected <- apply(
    missing_matrix,
    1,
    any
  )
  
  result <- data[affected, , drop = FALSE]
  
  result$missing_design_variables <-
    apply(
      missing_matrix[affected, , drop = FALSE],
      1,
      function(x) paste(
        vars_present[x],
        collapse = ", "
      )
    )
  
  result
}


# ------------------------------------------------------------
# 3. Identify possible unphenotyped individuals
# ------------------------------------------------------------

identify_unphenotyped <- function(data,
                                  genotype_col,
                                  trait_cols) {
  
  trait_cols <- intersect(
    trait_cols,
    names(data)
  )
  
  if (length(trait_cols) == 0) {
    return(data.frame())
  }
  
  trait_missing <- rowSums(
    is.na(data[trait_cols])
  )
  
  result <- data[
    !is.na(data[[genotype_col]]) &
      trait_missing == length(trait_cols),
    ,
    drop = FALSE
  ]
  
  if (nrow(result) > 0) {
    
    result$missing_all_traits <- TRUE
    
    result$missingness_action <-
      "REVIEW_FOR_GENOMIC_PREDICTION"
    
  }
  
  result
}


# ------------------------------------------------------------
# 4. Trait missingness by environment
# ------------------------------------------------------------

missing_by_environment_trait <- function(data,
                                         environment_col,
                                         trait_cols) {
  
  trait_cols <- intersect(
    trait_cols,
    names(data)
  )
  
  if (!environment_col %in% names(data)) {
    return(data.frame())
  }
  
  result <- lapply(
    trait_cols,
    function(trait) {
      
      tmp <- aggregate(
        is.na(data[[trait]]) ~ data[[environment_col]],
        FUN = mean
      )
      
      names(tmp) <- c(
        "environment",
        "pct_missing"
      )
      
      tmp$trait <- trait
      
      tmp
    }
  )
  
  result <- do.call(
    rbind,
    result
  )
  
  result$pct_missing <-
    result$pct_missing * 100
  
  result[, c(
    "environment",
    "trait",
    "pct_missing"
  )]
}


# ------------------------------------------------------------
# 5. Missingness by environment / replicate / block
# ------------------------------------------------------------

missing_by_design <- function(data,
                              environment_col,
                              grouping_cols,
                              trait_cols) {
  
  grouping_cols <- intersect(
    grouping_cols,
    names(data)
  )
  
  trait_cols <- intersect(
    trait_cols,
    names(data)
  )
  
  if (length(grouping_cols) == 0) {
    return(data.frame())
  }
  
  results <- list()
  
  counter <- 1
  
  for (trait in trait_cols) {
    
    tmp <- data
    
    tmp$.missing <- is.na(
      tmp[[trait]]
    )
    
    summary <- aggregate(
      tmp$.missing,
      by = tmp[grouping_cols],
      FUN = mean
    )
    
    names(summary)[
      ncol(summary)
    ] <- "pct_missing"
    
    summary$pct_missing <-
      summary$pct_missing * 100
    
    summary$trait <- trait
    
    results[[counter]] <- summary
    
    counter <- counter + 1
  }
  
  result <- do.call(
    rbind,
    results
  )
  
  result
}


# ------------------------------------------------------------
# 6. Detect high missingness
# ------------------------------------------------------------

detect_high_missingness <- function(
    env_trait_table,
    warning_threshold = 10,
    critical_threshold = 30) {
  
  result <- env_trait_table
  
  result$status <- ifelse(
    result$pct_missing >= critical_threshold,
    "HIGH",
    ifelse(
      result$pct_missing >= warning_threshold,
      "WARNING",
      "OK"
    )
  )
  
  result
}


# ------------------------------------------------------------
# 7. Detect environment missingness
# ------------------------------------------------------------

check_environment_missingness <- function(
    data,
    environment_col) {
  
  if (!environment_col %in% names(data)) {
    
    return(
      list(
        available = FALSE,
        missing = TRUE,
        n_missing = NA_integer_,
        pct_missing = NA_real_,
        allow_gxe = FALSE
      )
    )
  }
  
  n_missing <- sum(
    is.na(data[[environment_col]]) |
      data[[environment_col]] == ""
  )
  
  pct_missing <-
    round(100 * n_missing / nrow(data),4)
  
  list(
    available = TRUE,
    missing = n_missing > 0,
    n_missing = n_missing,
    pct_missing = pct_missing,
    allow_gxe = n_missing == 0
  )
}


# ------------------------------------------------------------
# 8. Row x Column missingness grid
# ------------------------------------------------------------

missing_grid <- function(data,
                         environment_col,
                         row_col,
                         column_col,
                         trait) {
  
  required <- c(
    environment_col,
    row_col,
    column_col,
    trait
  )
  
  if (!all(required %in% names(data))) {
    return(NULL)
  }
  
  result <- data[, required]
  
  names(result) <- c(
    "environment",
    "row",
    "column",
    "missing"
  )
  
  result$missing <- is.na(
    result$missing
  )
  
  result
}


# ------------------------------------------------------------
# 9. Overall trait missingness
# ------------------------------------------------------------

trait_missingness <- function(data,
                              trait_cols) {
  
  trait_cols <- intersect(
    trait_cols,
    names(data)
  )
  
  result <- data.frame(
    trait = trait_cols,
    n_total = nrow(data),
    n_missing = sapply(
      data[trait_cols],
      function(x) sum(is.na(x))
    )
  )
  
  result$pct_missing <-
    100 *
    result$n_missing /
    result$n_total
  
  result$status <- ifelse(
    result$pct_missing >= 30,
    "HIGH",
    ifelse(
      result$pct_missing >= 10,
      "WARNING",
      "OK"
    )
  )
  
  result
}


# ------------------------------------------------------------
# 10. Create complete missingness object
# ------------------------------------------------------------

run_missingness_analysis <- function(
    data,
    genotype_col,
    environment_col,
    trait_cols,
    design_vars,
    replicate_col = NULL,
    block_col = NULL,
    row_col = NULL,
    column_col = NULL,
    warning_threshold = 10,
    critical_threshold = 30) {
  
  # ------------------------------------------
  # Design
  # ------------------------------------------
  
  design <- check_design_missingness(
    data,
    design_vars
  )
  
  design_rows <- get_design_missing_rows(
    data,
    design_vars
  )
  
  # ------------------------------------------
  # Unphenotyped
  # ------------------------------------------
  
  unphenotyped <- identify_unphenotyped(
    data,
    genotype_col,
    trait_cols
  )
  
  # ------------------------------------------
  # Environment
  # ------------------------------------------
  
  environment <- check_environment_missingness(
    data,
    environment_col
  )
  
  # ------------------------------------------
  # Trait missingness
  # ------------------------------------------
  
  traits <- trait_missingness(
    data,
    trait_cols
  )
  
  # ------------------------------------------
  # Environment × Trait
  # ------------------------------------------
  
  env_trait <- missing_by_environment_trait(
    data,
    environment_col,
    trait_cols
  )
  
  env_trait <- detect_high_missingness(
    env_trait,
    warning_threshold,
    critical_threshold
  )
  
  # ------------------------------------------
  # Replicate / Block
  # ------------------------------------------
  
  grouping_cols <- c(
    environment_col,
    replicate_col,
    block_col
  )
  
  grouping_cols <- grouping_cols[
    !is.na(grouping_cols)
  ]
  
  grouping_cols <- grouping_cols[
    grouping_cols != ""
  ]
  
  design_trait <- missing_by_design(
    data,
    environment_col,
    grouping_cols,
    trait_cols
  )
  
  # ------------------------------------------
  # Pipeline rules
  # ------------------------------------------
  
  pipeline <- list(
    
    allow_STA = TRUE,
    
    allow_MTA = TRUE,
    
    allow_GEBV = nrow(unphenotyped) > 0,
    
    allow_GxE = environment$allow_gxe
  )
  
  # ------------------------------------------
  # Return
  # ------------------------------------------
  
  list(
    
    data = data,
    
    design = design,
    
    design_rows = design_rows,
    
    unphenotyped = unphenotyped,
    
    environment = environment,
    
    traits = traits,
    
    env_trait = env_trait,
    
    design_trait = design_trait,
    
    pipeline = pipeline,
    
    decisions = list(
      
      genomic_prediction = character(),
      
      remove_rows = character(),
      
      keep_environments = unique(
        data[[environment_col]]
      ),
      
      keep_traits = trait_cols
    )
  )
}
