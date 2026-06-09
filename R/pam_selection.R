coalesce_chr <- function(...) {
  x <- list(...)
  out <- as.character(x[[1]])
  
  if (length(x) > 1) {
    for (i in 2:length(x)) {
      xi <- as.character(x[[i]])
      replace_idx <- is.na(out) | out == ""
      out[replace_idx] <- xi[replace_idx]
    }
  }
  
  out
}


#' Check trait quality before running initial selection
#'
#' Returns a list of flagged traits with reasons. Called by UI to show
#' confirmation dialog before proceeding.
#' @export
checkTraitQuality <- function(args, dt_object) {
  
  check_entry_type_value <- args$checkEntryTypeValue
  if (is.null(check_entry_type_value) || !nzchar(check_entry_type_value)) {
    check_entry_type_value <- NULL
  }
  
  preds <- dt_object$predictions
  mta_preds <- preds[preds$analysisId == args$mtaStamp & preds$effectType == "designation", ]
  mta_preds <- mta_preds[mta_preds$trait %in% args$traitsToEvaluate, ]
  
  # Deduplicate (same as runInitialProdAdv)
  dup_key <- paste(mta_preds$designation, mta_preds$trait, sep = "|||")
  if (any(duplicated(dup_key))) {
    agg_cols <- c("predictedValue")
    if ("reliability" %in% colnames(mta_preds)) agg_cols <- c(agg_cols, "reliability")
    if ("stdError" %in% colnames(mta_preds)) agg_cols <- c(agg_cols, "stdError")
    mta_preds <- do.call(rbind, lapply(split(mta_preds, dup_key), function(x) {
      row1 <- x[1, , drop = FALSE]
      for (col in agg_cols) {
        if (col %in% colnames(x)) row1[[col]] <- mean(x[[col]], na.rm = TRUE)
      }
      row1
    }))
    rownames(mta_preds) <- NULL
  }
  
  # Filter to candidates only (same as runInitialProdAdv)
  if (!is.null(args$selectedCandidates)) {
    if (is.null(check_entry_type_value)) {
      mta_preds <- mta_preds[mta_preds$designation %in% args$selectedCandidates, ]
    } else {
      mta_preds <- mta_preds[
        mta_preds$designation %in% args$selectedCandidates |
          mta_preds$entryType == check_entry_type_value, ]
    }
  }
  
  if (is.null(check_entry_type_value)) {
    candidate_preds <- mta_preds
  } else {
    candidate_preds <- mta_preds[mta_preds$entryType != check_entry_type_value, ]
  }
  
  reliability_threshold <- 0.20
  reliability_pct_threshold <- 0.60
  rsr_threshold <- 1.0
  
  flagged_traits <- character(0)
  flag_reasons <- list()
  
  for (t in args$traitsToEvaluate) {
    trait_cand <- candidate_preds[candidate_preds$trait == t, ]
    reasons <- character(0)
    
    if ("reliability" %in% colnames(trait_cand)) {
      trait_rel <- trait_cand$reliability[!is.na(trait_cand$reliability)]
      if (length(trait_rel) > 0) {
        pct_low <- mean(trait_rel < reliability_threshold)
        if (pct_low > reliability_pct_threshold) {
          reasons <- c(reasons, sprintf("low reliability (%.0f%% below %.2f)", pct_low * 100, reliability_threshold))
        }
      }
    }
    
    if ("stdError" %in% colnames(trait_cand)) {
      blups <- trait_cand$predictedValue[!is.na(trait_cand$predictedValue)]
      ses <- trait_cand$stdError[!is.na(trait_cand$stdError)]
      if (length(blups) > 2 && length(ses) > 0 && mean(ses) > 0) {
        rsr <- sd(blups) / mean(ses)
        if (rsr < rsr_threshold) {
          reasons <- c(reasons, sprintf("ranking not separable (RSR=%.2f)", rsr))
        }
      }
    }
    
    if (length(reasons) > 0) {
      flagged_traits <- c(flagged_traits, t)
      flag_reasons[[t]] <- paste(reasons, collapse = "; ")
    }
  }
  
  list(
    flagged_traits = flagged_traits,
    flag_reasons = flag_reasons,
    has_issues = length(flagged_traits) > 0
  )
}


runInitialProdAdv <- function(analysisId = as.numeric(Sys.time()),
                             analysisIdName,
                             args,
                             dt_object){

  check_entry_type_value <- args$checkEntryTypeValue
  if (is.null(check_entry_type_value) || !nzchar(check_entry_type_value)) {
    check_entry_type_value <- NULL
  }

  #Filter MTA predictions
  preds <- dt_object$predictions
  mta_preds <- preds[preds$analysisId == args$mtaStamp & preds$effectType == "designation",]
  mta_preds <- mta_preds[mta_preds$trait %in% args$traitsToEvaluate,]
  mta_preds <- mta_preds[order(mta_preds$designation), ]

  # Deduplicate: keep one row per designation × trait (average if duplicated)
  # This handles cases where MTA stores multiple predictions per designation
  dup_key <- paste(mta_preds$designation, mta_preds$trait, sep = "|||")
  if (any(duplicated(dup_key))) {
    agg_cols <- c("predictedValue")
    if ("reliability" %in% colnames(mta_preds)) agg_cols <- c(agg_cols, "reliability")
    if ("stdError" %in% colnames(mta_preds)) agg_cols <- c(agg_cols, "stdError")
    
    # Keep first occurrence for non-numeric columns, average for numeric
    mta_preds_dedup <- do.call(rbind, lapply(split(mta_preds, dup_key), function(x) {
      row1 <- x[1, , drop = FALSE]
      for (col in agg_cols) {
        if (col %in% colnames(x)) {
          row1[[col]] <- mean(x[[col]], na.rm = TRUE)
        }
      }
      row1
    }))
    rownames(mta_preds_dedup) <- NULL
    mta_preds <- mta_preds_dedup[order(mta_preds_dedup$designation), ]
  }

  #Filter treatments according to pre-selection
  if (is.null(check_entry_type_value)) {
    mta_preds <- mta_preds[
      mta_preds$designation %in% args$selectedCandidates,
    ]
  } else {
    mta_preds <- mta_preds[
      mta_preds$designation %in% args$selectedCandidates |
        mta_preds$entryType == check_entry_type_value,
    ]
  }

  # ---- Trait quality assessment ----
  # Compute reliability and RSR metrics per trait (informational — no exclusion)
  # Flagged traits are reported back to UI for user decision
  reliability_threshold <- 0.20
  reliability_pct_threshold <- 0.60
  rsr_threshold <- 1.0
  flagged_traits <- character(0)
  flag_reasons <- list()

  # Identify candidates for quality check
  if (is.null(check_entry_type_value)) {
    candidate_preds <- mta_preds
  } else {
    candidate_preds <- mta_preds[mta_preds$entryType != check_entry_type_value, ]
  }

  for (t in args$traitsToEvaluate) {
    trait_cand <- candidate_preds[candidate_preds$trait == t, ]
    reasons <- character(0)

    # Check 1: Reliability
    if ("reliability" %in% colnames(trait_cand)) {
      trait_rel <- trait_cand$reliability[!is.na(trait_cand$reliability)]
      if (length(trait_rel) > 0) {
        pct_low <- mean(trait_rel < reliability_threshold)
        if (pct_low > reliability_pct_threshold) {
          reasons <- c(reasons, sprintf("low reliability (%.0f%% of candidates below %.2f)", pct_low * 100, reliability_threshold))
        }
      }
    }

    # Check 2: Ranking Separability Ratio (RSR)
    if ("stdError" %in% colnames(trait_cand)) {
      blups <- trait_cand$predictedValue[!is.na(trait_cand$predictedValue)]
      ses <- trait_cand$stdError[!is.na(trait_cand$stdError)]
      if (length(blups) > 2 && length(ses) > 0 && mean(ses) > 0) {
        rsr <- sd(blups) / mean(ses)
        if (rsr < rsr_threshold) {
          reasons <- c(reasons, sprintf("ranking not separable (RSR=%.2f, threshold=%.1f)", rsr, rsr_threshold))
        }
      }
    }

    if (length(reasons) > 0) {
      flagged_traits <- c(flagged_traits, t)
      flag_reasons[[t]] <- paste(reasons, collapse = "; ")
    }
  }

  # Use only the traits the user chose to keep (after seeing the flags)
  # args$traitsToUse is set by the UI: either all traits or with flagged ones removed
  traits_to_use <- args$traitsToEvaluate
  if (!is.null(args$excludeTraits) && length(args$excludeTraits) > 0) {
    traits_to_use <- setdiff(traits_to_use, args$excludeTraits)
  }
  if (length(traits_to_use) == 0) {
    stop("No traits remaining after exclusion. Cannot proceed with selection.")
  }

  # ---- Build prediction and reliability matrices ----
  all_designations <- unique(mta_preds$designation)
  n_des <- length(all_designations)

  mta_preds_matrix <- matrix(nrow = n_des, ncol = length(traits_to_use))
  rownames(mta_preds_matrix) <- all_designations
  colnames(mta_preds_matrix) <- traits_to_use

  # Reliability matrix (per-individual, per-trait) for index weighting
  rel_matrix <- matrix(1, nrow = n_des, ncol = length(traits_to_use))
  rownames(rel_matrix) <- all_designations
  colnames(rel_matrix) <- traits_to_use

  for (idx in seq_along(traits_to_use)) {
    t <- traits_to_use[idx]
    trait_data <- mta_preds[mta_preds$trait == t, ]
    mta_preds_matrix[, idx] <- trait_data$predictedValue[match(all_designations, trait_data$designation)]

    # Populate reliability matrix if available
    if ("reliability" %in% colnames(trait_data)) {
      rel_vals <- trait_data$reliability[match(all_designations, trait_data$designation)]
      rel_vals[is.na(rel_vals)] <- 0
      # Clamp to [0, 1]
      rel_vals <- pmax(0, pmin(1, rel_vals))
      rel_matrix[, idx] <- rel_vals
    }
  }

  # ---- Compute reliability-weighted selection index ----
  # Strategy: scale traits, then weight each individual's scaled BLUP by sqrt(reliability)
  # This penalizes poorly-estimated individuals: their contribution to the index is dampened
  # sqrt is used rather than raw reliability to avoid being too aggressive (sqrt(0.5) ≈ 0.71)
  index_weights <- args$customWeights[traits_to_use]
  if (is.null(index_weights) || length(index_weights) != length(traits_to_use)) {
    index_weights <- rep(1, length(traits_to_use))
    names(index_weights) <- traits_to_use
  }

  scaled_matrix <- scale(mta_preds_matrix)
  # Apply reliability penalty: element-wise multiplication with sqrt(reliability)
  reliability_penalized <- scaled_matrix * sqrt(rel_matrix)
  index_preds <- reliability_penalized %*% index_weights

  # Build result data frame
  mta_preds_matrix_full <- cbind(mta_preds_matrix, index_preds)
  colnames(mta_preds_matrix_full) <- c(traits_to_use, "index_value")

  mta_preds_long <- as.data.frame(mta_preds_matrix_full)
  mta_preds_long$designation <- rownames(mta_preds_matrix_full)
  # Get entryType for each designation via match (handles duplicates safely)
  entry_type_lookup <- mta_preds[!duplicated(mta_preds$designation), c("designation", "entryType"), drop = FALSE]
  mta_preds_long$entryType <- entry_type_lookup$entryType[match(mta_preds_long$designation, entry_type_lookup$designation)]

  # ---- Apply index selection (% or number) ----
  if (is.null(check_entry_type_value)) {
    n_candidates <- nrow(mta_preds_long)
    candidate_indices <- seq_len(nrow(mta_preds_long))
  } else {
    # Case-insensitive comparison for entryType
    candidate_indices <- which(toupper(trimws(mta_preds_long$entryType)) != toupper(trimws(check_entry_type_value)))
    n_candidates <- length(candidate_indices)
  }

  if (n_candidates == 0) {
    # Fallback: treat all as candidates if no match found
    n_candidates <- nrow(mta_preds_long)
    candidate_indices <- seq_len(nrow(mta_preds_long))
  }

  if (!is.null(args$nSelected) && args$nSelected > 0) {
    n_selected <- min(args$nSelected, n_candidates)
  } else {
    top_pct <- if (!is.null(args$topPctSelected)) args$topPctSelected else 20
    n_selected <- as.integer(n_candidates * (top_pct / 100))
  }
  n_selected <- max(c(1, n_selected))
  
  message(sprintf("[PAM] n_candidates=%d, n_selected=%d, topPct=%s, nSelected_arg=%s, check_entry_type='%s'",
                  n_candidates, n_selected, 
                  as.character(args$topPctSelected), as.character(args$nSelected),
                  as.character(check_entry_type_value)))

  index_order <- order(mta_preds_long$index_value[candidate_indices], decreasing = TRUE)
  selected_positions <- candidate_indices[index_order[1:n_selected]]
  selected_designations <- mta_preds_long$designation[selected_positions]

  mta_preds_long$index_decision <- "NOT SELECTED"
  mta_preds_long$index_decision[mta_preds_long$designation %in% selected_designations] <- "SELECTED"
  if (!is.null(check_entry_type_value)) {
    mta_preds_long$index_decision[toupper(trimws(mta_preds_long$entryType)) == toupper(trimws(check_entry_type_value))] <- "CHECK"
  }

  # ---- Apply trait-specific thresholds ----
  trait_decision <- matrix(nrow = nrow(mta_preds_long), ncol = length(traits_to_use))
  colnames(trait_decision) <- paste0(traits_to_use, "_trait_decision")

  for (t in traits_to_use) {
    trait_rules <- args$traitRules[[t]]
    trait_value <- mta_preds_long[, t]

    # Default: no threshold set (all pass)
    if (is.null(trait_rules) || is.null(trait_rules$ruleType) || trait_rules$ruleType == "None") {
      decision <- rep("SELECTED", length(trait_value))
    } else if (trait_rules$ruleType == "Threshold") {
      if (trait_rules$direction == "Higher is better") {
        decision <- ifelse(trait_value >= trait_rules$threshold, "SELECTED", "NOT SELECTED")
      } else {
        decision <- ifelse(trait_value <= trait_rules$threshold, "SELECTED", "NOT SELECTED")
      }
    } else if (trait_rules$ruleType == "Acceptable range") {
      decision <- ifelse(trait_value >= trait_rules$minValue & trait_value <= trait_rules$maxValue,
                         "SELECTED", "NOT SELECTED")
    } else if (trait_rules$ruleType == "% over check") {
      if (is.null(check_entry_type_value)) {
        stop("Checks are required when using a '% over check' trait threshold.")
      }
      check_value <- mean(trait_value[mta_preds_long$entryType == check_entry_type_value], na.rm = TRUE)
      if (trait_rules$direction == "Higher is better") {
        target <- check_value + check_value * (trait_rules$threshold / 100)
        decision <- ifelse(trait_value >= target, "SELECTED", "NOT SELECTED")
      } else {
        target <- check_value - check_value * (trait_rules$threshold / 100)
        decision <- ifelse(trait_value <= target, "SELECTED", "NOT SELECTED")
      }
    } else if (trait_rules$ruleType == "% over mean") {
      if (is.null(check_entry_type_value)) {
        mean_value <- mean(trait_value, na.rm = TRUE)
      } else {
        mean_value <- mean(trait_value[mta_preds_long$entryType != check_entry_type_value], na.rm = TRUE)
      }
      if (trait_rules$direction == "Higher is better") {
        target <- mean_value + mean_value * (trait_rules$threshold / 100)
        decision <- ifelse(trait_value >= target, "SELECTED", "NOT SELECTED")
      } else {
        target <- mean_value - mean_value * (trait_rules$threshold / 100)
        decision <- ifelse(trait_value <= target, "SELECTED", "NOT SELECTED")
      }
    } else {
      decision <- rep("SELECTED", length(trait_value))
    }

    if (!is.null(check_entry_type_value)) {
      decision[toupper(trimws(mta_preds_long$entryType)) == toupper(trimws(check_entry_type_value))] <- "CHECK"
    }

    trait_decision[, paste0(t, "_trait_decision")] <- decision
  }

  mta_preds_long <- cbind(mta_preds_long, trait_decision)

  # ---- Combine index + threshold decisions ----
  # Final initial decision: must pass index selection AND all trait thresholds
  decision_cols <- grep("_decision", colnames(mta_preds_long), value = TRUE)
  decision_set <- mta_preds_long[, decision_cols, drop = FALSE]
  initial_decision <- apply(decision_set, 1, function(x) {
    if (any(x == "CHECK")) return("CHECK")
    if (any(x == "NOT SELECTED")) return("NOT SELECTED")
    return("SELECTED")
  })

  mta_preds_long$initial_decision <- initial_decision

  # ---- Create modifications table ----
  modifications <- data.frame(
    module = "Init_prodAdv",
    analysisId = analysisId,
    designation = mta_preds_long$designation,
    reason = "initial_selection",
    value = mta_preds_long$initial_decision,
    stringsAsFactors = FALSE
  )

  # ---- Create modeling table ----
  .null_or_name <- function(x, name) {
    if (is.null(x)) NULL else name
  }

  modeling <- NULL
  for (t in traits_to_use) {
    weight_value <- index_weights[t]
    trait_rules_t <- args$traitRules[[t]]

    arg_values <- c(
      args$mtaStamp,
      args$staStamp,
      "Weighted index",
      as.character(args$topPctSelected),
      as.character(args$nSelected),
      as.character(weight_value),
      "TRUE",
      if (!is.null(trait_rules_t$ruleType)) trait_rules_t$ruleType else "None",
      .null_or_name(trait_rules_t$direction, trait_rules_t$direction),
      .null_or_name(trait_rules_t$threshold, as.character(trait_rules_t$threshold)),
      .null_or_name(trait_rules_t$referenceCheck, trait_rules_t$referenceCheck),
      .null_or_name(trait_rules_t$minValue, as.character(trait_rules_t$minValue)),
      .null_or_name(trait_rules_t$maxValue, as.character(trait_rules_t$maxValue))
    )

    arg_param <- c(
      "mta_stamp",
      "sta_stamp",
      "decision_logic",
      "top_perc_selected",
      .null_or_name(args$nSelected, "n_selected"),
      "index_weight",
      "scaled_traits",
      "trait_rule_type",
      .null_or_name(trait_rules_t$direction, "direction"),
      .null_or_name(trait_rules_t$threshold, "threshold"),
      .null_or_name(trait_rules_t$referenceCheck, "reference_check"),
      .null_or_name(trait_rules_t$minValue, "min_value"),
      .null_or_name(trait_rules_t$maxValue, "max_value")
    )

    # Remove NULLs (pairs must match)
    keep <- !sapply(arg_param, is.null)
    arg_param <- unlist(arg_param[keep])
    arg_values <- unlist(arg_values[keep])

    trait_model <- data.frame(
      module = "Init_prodAdv",
      analysisId = analysisId,
      trait = t,
      environment = NA,
      parameter = arg_param,
      value = arg_values,
      stringsAsFactors = FALSE
    )

    modeling <- if (is.null(modeling)) trait_model else rbind(modeling, trait_model)
  }

  # Record which traits were actually excluded by user choice
  if (!is.null(args$excludeTraits) && length(args$excludeTraits) > 0) {
    excl_rows <- data.frame(
      module = "Init_prodAdv",
      analysisId = analysisId,
      trait = args$excludeTraits,
      environment = NA,
      parameter = "user_excluded_trait",
      value = "TRUE",
      stringsAsFactors = FALSE
    )
    modeling <- rbind(modeling, excl_rows)
  }

  # Record flagged traits that were KEPT (not excluded) — informational only
  kept_flagged <- setdiff(flagged_traits, if (!is.null(args$excludeTraits)) args$excludeTraits else character(0))
  if (length(kept_flagged) > 0) {
    flag_rows <- data.frame(
      module = "Init_prodAdv",
      analysisId = analysisId,
      trait = kept_flagged,
      environment = NA,
      parameter = "flagged_low_quality",
      value = sapply(kept_flagged, function(t) flag_reasons[[t]]),
      stringsAsFactors = FALSE
    )
    modeling <- rbind(modeling, flag_rows)
  }

  # ---- Create status table ----
  status <- data.frame(
    module = "Init_prodAdv",
    analysisId = analysisId,
    analysisIdName = analysisIdName,
    stringsAsFactors = FALSE
  )

  # ---- Update data object ----
  if (is.null(dt_object$modifications$selection)) {
    dt_object$modifications$selection <- modifications
  } else {
    dt_object$modifications$selection <- rbind(dt_object$modifications$selection, modifications)
  }

  dt_object$modeling <- rbind(dt_object$modeling, modeling)
  dt_object$status <- rbind(dt_object$status, status)

  # Attach flagged traits info as attribute for UI warning
  attr(dt_object, "pam_flagged_traits") <- flagged_traits
  attr(dt_object, "pam_flag_reasons") <- flag_reasons

  dt_object

}

savePlotProdAdvSelection <- function(
    analysisId = as.numeric(Sys.time()),
    analysisIdName = NULL,
    initialSelectionStamp,
    plotSelectionStamp = NULL,
    manual_decisions,
    dt_object
) {
  
  if (!is.list(dt_object)) {
    stop("'dt_object' must be a data object list.")
  }
  
  if (is.null(dt_object$modifications$selection)) {
    stop("No existing selection modifications table found in dt_object.")
  }
  
  if (is.null(dt_object$status)) {
    stop("No status table found in dt_object.")
  }
  
  if (!is.data.frame(manual_decisions) || nrow(manual_decisions) == 0) {
    stop("No manual decisions were provided.")
  }
  
  required_cols <- c("designation", "plot_decision")
  missing_cols <- setdiff(required_cols, colnames(manual_decisions))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in manual_decisions: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  valid_values <- c("SELECTED", "NOT SELECTED")
  if (!all(manual_decisions$plot_decision %in% valid_values)) {
    stop("manual_decisions$plot_decision must contain only 'SELECTED' or 'NOT SELECTED'.")
  }
  
  base_selection <- dt_object$modifications$selection
  base_selection <- base_selection[
    base_selection$analysisId %in% initialSelectionStamp &
      base_selection$module == "Init_prodAdv" &
      base_selection$reason == "initial_selection",
    ,
    drop = FALSE
  ]
  
  if (nrow(base_selection) == 0) {
    stop("No initial selection records found for the provided initialSelectionStamp.")
  }
  
  names(base_selection)[names(base_selection) == "value"] <- "base_decision"
  base_selection <- unique(base_selection[, c("designation", "base_decision"), drop = FALSE])
  
  if (!is.null(plotSelectionStamp) && nzchar(plotSelectionStamp)) {
    previous_plot <- dt_object$modifications$selection
    previous_plot <- previous_plot[
      previous_plot$analysisId %in% plotSelectionStamp &
        previous_plot$module == "Plot_prodAdv" &
        previous_plot$reason == "manual_plot_selection",
      ,
      drop = FALSE
    ]
    
    if (nrow(previous_plot) > 0) {
      names(previous_plot)[names(previous_plot) == "value"] <- "base_decision"
      previous_plot <- unique(previous_plot[, c("designation", "base_decision"), drop = FALSE])
      base_selection <- previous_plot
    }
  }
  
  merged_decisions <- merge(
    base_selection,
    manual_decisions,
    by = "designation",
    all.y = TRUE
  )
  
  if (nrow(merged_decisions) == 0) {
    stop("None of the manual decisions matched the candidate designations.")
  }
  
  modifications <- data.frame(
    module = "Plot_prodAdv",
    analysisId = analysisId,
    designation = merged_decisions$designation,
    reason = "manual_plot_selection",
    value = merged_decisions$plot_decision,
    stringsAsFactors = FALSE
  )
  
  status <- data.frame(
    module = "Plot_prodAdv",
    analysisId = analysisId,
    analysisIdName = analysisIdName,
    stringsAsFactors = FALSE
  )
  
  dt_object$modifications$selection <- rbind(
    dt_object$modifications$selection,
    modifications
  )
  
  dt_object$status <- rbind(
    dt_object$status,
    status
  )
  
  dt_object
}



build_prodadv_decision_table_data <- function(dt,
                                              initial_stamp,
                                              plot_stamp = "__none__",
                                              final_stamp = "__none__",
                                              final_overrides = NULL) {
  
  if (is.null(dt)) {
    stop("dt is required.")
  }
  
  if (is.null(initial_stamp) || !nzchar(initial_stamp)) {
    stop("initial_stamp is required.")
  }
  
  if (is.null(dt$modeling)) {
    stop("dt$modeling is missing.")
  }
  
  if (is.null(dt$predictions)) {
    stop("dt$predictions is missing.")
  }
  
  if (is.null(dt$modifications) ||
      is.null(dt$modifications$selection)) {
    stop("dt$modifications$selection is missing.")
  }
  
  modeling_init <- dt$modeling[
    dt$modeling$analysisId %in% initial_stamp &
      dt$modeling$module == "Init_prodAdv",
    ,
    drop = FALSE
  ]
  
  if (nrow(modeling_init) == 0) {
    stop("No Init_prodAdv modeling rows found for the selected initial stamp.")
  }
  
  selected_traits <- unique(modeling_init$trait)
  selected_traits <- selected_traits[!is.na(selected_traits) & nzchar(selected_traits)]
  
  # Exclude traits marked as user_excluded_trait
  excluded_trait_rows <- modeling_init[modeling_init$parameter == "user_excluded_trait", , drop = FALSE]
  if (nrow(excluded_trait_rows) > 0) {
    selected_traits <- setdiff(selected_traits, excluded_trait_rows$trait)
  }
  
  if (length(selected_traits) == 0) {
    stop("No selected traits found in Init_prodAdv modeling rows.")
  }
  
  mta_stamp <- modeling_init$value[modeling_init$parameter == "mta_stamp"][1]
  
  if (is.na(mta_stamp) || !nzchar(mta_stamp)) {
    stop("No mta_stamp found in Init_prodAdv modeling rows.")
  }
  
  preds <- dt$predictions[
    dt$predictions$analysisId %in% mta_stamp &
      dt$predictions$trait %in% selected_traits,
    ,
    drop = FALSE
  ]
  
  if (nrow(preds) == 0) {
    stop("No predictions found for the selected MTA stamp and traits.")
  }
  
  pred_wide <- reshape(
    preds[, c("designation", "trait", "predictedValue"), drop = FALSE],
    idvar = "designation",
    timevar = "trait",
    direction = "wide"
  )
  
  names(pred_wide) <- sub("^predictedValue\\.", "", names(pred_wide))
  
  # Compute index_value from stored weights in modeling table (with reliability weighting)
  # IMPORTANT: only use designations in the initial selection for scaling (same population as runInitialProdAdv)
  weight_rows <- modeling_init[modeling_init$parameter == "index_weight", , drop = FALSE]
  
  init_sel <- dt$modifications$selection[
    dt$modifications$selection$analysisId %in% initial_stamp &
      dt$modifications$selection$module == "Init_prodAdv" &
      dt$modifications$selection$reason == "initial_selection",
    ,
    drop = FALSE
  ]
  
  if (nrow(init_sel) == 0) {
    stop("No initial selection rows found for the selected initial stamp.")
  }
  
  init_sel <- unique(init_sel[, c("designation", "value"), drop = FALSE])
  names(init_sel)[names(init_sel) == "value"] <- "initial_decision"
  review_designations <- unique(init_sel$designation)
  
  # Filter pred_wide to only review designations BEFORE computing index
  pred_wide <- pred_wide[pred_wide$designation %in% review_designations, , drop = FALSE]
  
  if (nrow(weight_rows) > 0) {
    index_weights <- as.numeric(weight_rows$value)
    names(index_weights) <- weight_rows$trait
    avail_traits <- intersect(names(index_weights), colnames(pred_wide))
    if (length(avail_traits) > 0) {
      trait_matrix <- as.matrix(pred_wide[, avail_traits, drop = FALSE])
      scaled_matrix <- scale(trait_matrix)
      w <- index_weights[avail_traits]
      
      # Apply reliability weighting (same as runInitialProdAdv)
      if ("reliability" %in% colnames(preds)) {
        preds_review <- preds[preds$designation %in% review_designations, , drop = FALSE]
        rel_data <- reshape(
          preds_review[, c("designation", "trait", "reliability"), drop = FALSE],
          idvar = "designation", timevar = "trait", direction = "wide"
        )
        names(rel_data) <- sub("^reliability\\.", "", names(rel_data))
        rel_avail <- intersect(avail_traits, colnames(rel_data))
        if (length(rel_avail) == length(avail_traits)) {
          rel_matrix <- as.matrix(rel_data[match(pred_wide$designation, rel_data$designation), avail_traits, drop = FALSE])
          rel_matrix[is.na(rel_matrix)] <- 0
          rel_matrix <- pmax(0, pmin(1, rel_matrix))
          reliability_penalized <- scaled_matrix * sqrt(rel_matrix)
          pred_wide$index_value <- as.numeric(reliability_penalized %*% w)
        } else {
          pred_wide$index_value <- as.numeric(scaled_matrix %*% w)
        }
      } else {
        pred_wide$index_value <- as.numeric(scaled_matrix %*% w)
      }
    } else {
      pred_wide$index_value <- NA_real_
    }
  } else {
    pred_wide$index_value <- NA_real_
  }
  
  # Filter preds to review designations (init_sel and review_designations already defined above)
  preds <- preds[
    preds$designation %in% review_designations,
    ,
    drop = FALSE
  ]
  
  entry_type <- unique(preds[, c("designation", "entryType"), drop = FALSE])
  
  base_df <- merge(
    pred_wide,
    entry_type,
    by = "designation",
    all = FALSE
  )
  
  base_df <- merge(
    base_df,
    init_sel,
    by = "designation",
    all = FALSE
  )
  
  check_entry_type_value <- unique(
    preds$entryType[
      init_sel$initial_decision[
        match(preds$designation, init_sel$designation)
      ] == "CHECK"
    ]
  )
  
  check_entry_type_value <- check_entry_type_value[!is.na(check_entry_type_value)]
  
  if (length(check_entry_type_value) == 0) {
    check_entry_type_value <- NULL
  } else {
    check_entry_type_value <- check_entry_type_value[[1]]
  }
  
  trait_decision_df <- build_prodadv_trait_decision_table(
    preds = preds[, c("designation", "trait", "predictedValue", "entryType"), drop = FALSE],
    modeling_init = modeling_init,
    check_entry_type_value = check_entry_type_value
  )
  
  base_df <- merge(
    base_df,
    trait_decision_df,
    by = "designation",
    all.x = TRUE
  )
  
  if (!is.null(plot_stamp) &&
      nzchar(plot_stamp) &&
      plot_stamp != "__none__") {
    
    plot_manual <- dt$modifications$selection[
      dt$modifications$selection$analysisId %in% plot_stamp &
        dt$modifications$selection$module == "Plot_prodAdv" &
        dt$modifications$selection$reason == "manual_plot_selection",
      ,
      drop = FALSE
    ]
    
    if (nrow(plot_manual) > 0) {
      plot_manual <- unique(plot_manual[, c("designation", "value"), drop = FALSE])
      names(plot_manual)[names(plot_manual) == "value"] <- "plot_decision"
      
      base_df <- merge(
        base_df,
        plot_manual,
        by = "designation",
        all.x = TRUE
      )
    } else {
      base_df$plot_decision <- NA_character_
    }
    
  } else {
    base_df$plot_decision <- NA_character_
  }
  
  if (!is.null(final_stamp) &&
      nzchar(final_stamp) &&
      final_stamp != "__none__") {
    
    final_manual <- dt$modifications$selection[
      dt$modifications$selection$analysisId %in% final_stamp &
        dt$modifications$selection$module == "Final_prodAdv" &
        dt$modifications$selection$reason == "final_manual_decision",
      ,
      drop = FALSE
    ]
    
    if (nrow(final_manual) > 0) {
      final_manual <- unique(final_manual[, c("designation", "value"), drop = FALSE])
      names(final_manual)[names(final_manual) == "value"] <- "final_manual_decision"
      
      base_df <- merge(
        base_df,
        final_manual,
        by = "designation",
        all.x = TRUE
      )
    } else {
      base_df$final_manual_decision <- NA_character_
    }
    
  } else {
    base_df$final_manual_decision <- NA_character_
  }
  
  if (is.null(final_overrides)) {
    final_overrides <- data.frame(
      designation = character(),
      final_decision = character(),
      stringsAsFactors = FALSE
    )
  }
  
  if (nrow(final_overrides) > 0) {
    base_df <- merge(
      base_df,
      final_overrides,
      by = "designation",
      all.x = TRUE
    )
  } else {
    base_df$final_decision <- NA_character_
  }
  
  base_df$final_decision_initial <- coalesce_chr(
    base_df$final_decision,
    base_df$final_manual_decision,
    base_df$plot_decision,
    base_df$initial_decision
  )
  
  # Sort by index_value descending (checks participate in ordering)
  if ("index_value" %in% colnames(base_df)) {
    base_df <- base_df[order(-as.numeric(base_df$index_value)), , drop = FALSE]
  } else {
    base_df <- base_df[order(base_df$designation), , drop = FALSE]
  }
  
  base_df
}

build_prodadv_review_plot_data <- function(dt,
                                           initial_stamp,
                                           plot_stamp = "__none__",
                                           final_stamp = "__none__") {
  
  if (is.null(dt)) {
    stop("dt is required.")
  }
  
  if (is.null(initial_stamp) || !nzchar(initial_stamp)) {
    stop("initial_stamp is required.")
  }
  
  if (is.null(dt$predictions)) {
    stop("dt$predictions is missing.")
  }
  
  if (is.null(dt$modifications)) {
    stop("dt$modifications is missing.")
  }
  
  if (is.null(dt$status)) {
    stop("dt$status is missing.")
  }
  
  if (is.null(dt$modifications$selection)) {
    stop("dt$modifications$selection is missing.")
  }
  
  if (is.null(dt$modeling)) {
    stop("dt$modeling is missing.")
  }
  
  init_sel <- dt$modifications$selection
  init_sel <- init_sel[
    init_sel$analysisId %in% initial_stamp &
      init_sel$module == "Init_prodAdv" &
      init_sel$reason == "initial_selection",
    ,
    drop = FALSE
  ]
  
  if (nrow(init_sel) == 0) {
    stop("No initial selection records found for the selected stamp.")
  }
  
  init_sel <- unique(init_sel[, c("designation", "value"), drop = FALSE])
  names(init_sel)[names(init_sel) == "value"] <- "initial_decision"
  review_designations <- unique(init_sel$designation)
  
  
  plot_sel <- NULL
  
  if (!is.null(plot_stamp) &&
      nzchar(plot_stamp) &&
      plot_stamp != "__none__") {
    
    plot_sel <- dt$modifications$selection
    plot_sel <- plot_sel[
      plot_sel$analysisId %in% plot_stamp &
        plot_sel$module == "Plot_prodAdv" &
        plot_sel$reason == "manual_plot_selection",
      ,
      drop = FALSE
    ]
    
    if (nrow(plot_sel) > 0) {
      plot_sel <- unique(plot_sel[, c("designation", "value"), drop = FALSE])
      names(plot_sel)[names(plot_sel) == "value"] <- "plot_decision"
    } else {
      plot_sel <- NULL
    }
  }
  
  final_sel <- NULL
  
  if (!is.null(final_stamp) &&
      nzchar(final_stamp) &&
      final_stamp != "__none__") {
    
    final_sel <- dt$modifications$selection
    final_sel <- final_sel[
      final_sel$analysisId %in% final_stamp &
        final_sel$module == "Final_prodAdv" &
        final_sel$reason == "final_manual_decision",
      ,
      drop = FALSE
    ]
    
    if (nrow(final_sel) > 0) {
      final_sel <- unique(final_sel[, c("designation", "value"), drop = FALSE])
      names(final_sel)[names(final_sel) == "value"] <- "final_decision"
    } else {
      final_sel <- NULL
    }
  }
  
  modeling_init <- dt$modeling
  modeling_init <- modeling_init[
    modeling_init$analysisId %in% initial_stamp &
      modeling_init$module == "Init_prodAdv",
    ,
    drop = FALSE
  ]
  
  if (nrow(modeling_init) == 0) {
    stop("No modeling records found for the selected initial selection stamp.")
  }
  
  selected_traits <- unique(modeling_init$trait[!is.na(modeling_init$trait)])
  selected_traits <- selected_traits[nzchar(selected_traits)]
  
  if (length(selected_traits) == 0) {
    stop("No traits found in the modeling table for this selection stamp.")
  }
  
  mta_stamp_used <- unique(modeling_init$value[modeling_init$parameter == "mta_stamp"])
  sta_stamp_used <- unique(modeling_init$value[modeling_init$parameter == "sta_stamp"])
  
  mta_stamp_used <- mta_stamp_used[!is.na(mta_stamp_used) & nzchar(mta_stamp_used)]
  sta_stamp_used <- sta_stamp_used[!is.na(sta_stamp_used) & nzchar(sta_stamp_used)]
  
  if (length(mta_stamp_used) == 0) {
    stop("No MTA stamp recorded for this initial selection.")
  }
  
  if (length(sta_stamp_used) == 0) {
    stop("No STA stamp recorded for this initial selection.")
  }
  
  preds <- dt$predictions
  
  mta_preds <- preds[
    preds$analysisId %in% mta_stamp_used &
      preds$effectType == "designation" &
      preds$trait %in% selected_traits,
    ,
    drop = FALSE
  ]
  
  sta_preds <- preds[
    preds$analysisId %in% sta_stamp_used &
      preds$trait %in% selected_traits,
    ,
    drop = FALSE
  ]
  
  mta_preds <- mta_preds[
    mta_preds$designation %in% review_designations,
    ,
    drop = FALSE
  ]
  
  sta_preds <- sta_preds[
    sta_preds$designation %in% review_designations,
    ,
    drop = FALSE
  ]
  
  if (nrow(mta_preds) == 0) {
    stop("No MTA predictions found for the selected initial selection.")
  }
  
  if (nrow(sta_preds) == 0) {
    stop("No STA predictions found for the selected initial selection.")
  }
  
  mta_cols <- c("designation", "trait", "predictedValue")
  
  for (extra_col in c("reliability", "stdError")) {
    if (extra_col %in% colnames(mta_preds)) {
      mta_cols <- c(mta_cols, extra_col)
    }
  }
  
  mta_preds <- mta_preds[, mta_cols, drop = FALSE]
  
  for (extra_col in c("reliability", "stdError")) {
    if (!extra_col %in% colnames(mta_preds)) {
      mta_preds[[extra_col]] <- NA_real_
    }
  }
  
  sta_cols <- c("designation", "environment", "trait", "predictedValue")
  
  for (extra_col in c("reliability", "stdError")) {
    if (extra_col %in% colnames(sta_preds)) {
      sta_cols <- c(sta_cols, extra_col)
    }
  }
  
  sta_preds <- sta_preds[, sta_cols, drop = FALSE]
  
  for (extra_col in c("reliability", "stdError")) {
    if (!extra_col %in% colnames(sta_preds)) {
      sta_preds[[extra_col]] <- NA_real_
    }
  }
  
  pred_wide <- reshape(
    mta_preds[, c("designation", "trait", "predictedValue"), drop = FALSE],
    idvar = "designation",
    timevar = "trait",
    direction = "wide"
  )
  
  names(pred_wide) <- sub("^predictedValue\\.", "", names(pred_wide))
  
  rel_wide <- reshape(
    mta_preds[, c("designation", "trait", "reliability"), drop = FALSE],
    idvar = "designation",
    timevar = "trait",
    direction = "wide"
  )
  
  names(rel_wide) <- sub("^reliability\\.", "reliability_", names(rel_wide))
  
  se_wide <- reshape(
    mta_preds[, c("designation", "trait", "stdError"), drop = FALSE],
    idvar = "designation",
    timevar = "trait",
    direction = "wide"
  )
  
  names(se_wide) <- sub("^stdError\\.", "stdError_", names(se_wide))
  
  review_df <- merge(
    pred_wide,
    init_sel,
    by = "designation",
    all.x = TRUE
  )
  
  review_df <- merge(
    review_df,
    rel_wide,
    by = "designation",
    all.x = TRUE
  )
  
  review_df <- merge(
    review_df,
    se_wide,
    by = "designation",
    all.x = TRUE
  )
  
  if (!is.null(plot_sel)) {
    review_df <- merge(
      review_df,
      plot_sel,
      by = "designation",
      all.x = TRUE
    )
  } else {
    review_df$plot_decision <- NA_character_
  }
  
  if (!is.null(final_sel)) {
    review_df <- merge(
      review_df,
      final_sel,
      by = "designation",
      all.x = TRUE
    )
  } else {
    review_df$final_decision <- NA_character_
  }
  
  review_df$plot_status <- coalesce_chr(
    review_df$final_decision,
    review_df$plot_decision,
    review_df$initial_decision
  )
  
  list(
    review_df = review_df,
    mta_long = mta_preds,
    sta_long = sta_preds,
    traits = selected_traits,
    initialSelectionStamp = initial_stamp,
    plotSelectionStamp = plot_stamp,
    finalSelectionStamp = final_stamp
  )
}


## Helpers for decision table

get_prodadv_trait_rules <- function(modeling_init) {
  stopifnot(is.data.frame(modeling_init))
  
  modeling_init <- modeling_init[modeling_init$module == "Init_prodAdv", , drop = FALSE]
  
  traits <- unique(modeling_init$trait)
  traits <- traits[!is.na(traits) & nzchar(traits)]
  
  out <- lapply(traits, function(tr) {
    rows <- modeling_init[modeling_init$trait == tr, , drop = FALSE]
    
    get_param <- function(param_name) {
      vals <- rows$value[rows$parameter == param_name]
      vals <- vals[!is.na(vals)]
      if (length(vals) == 0) return(NULL)
      vals[[1]]
    }
    
    as_num_or_null <- function(x) {
      if (is.null(x) || identical(x, "") || is.na(x)) return(NULL)
      as.numeric(x)
    }
    
    list(
      trait = tr,
      ruleType = get_param("trait_rule_type"),
      direction = get_param("direction"),
      threshold = as_num_or_null(get_param("threshold")),
      referenceCheck = get_param("reference_check"),
      minValue = as_num_or_null(get_param("min_value")),
      maxValue = as_num_or_null(get_param("max_value"))
    )
  })
  
  names(out) <- traits
  out
}


apply_prodadv_trait_rule <- function(trait_df, rule_def, check_entry_type_value = NULL) {
  stopifnot(is.data.frame(trait_df))
  stopifnot(all(c("designation", "predictedValue", "entryType") %in% colnames(trait_df)))
  
  trait_value <- trait_df$predictedValue
  decision <- rep(NA_character_, nrow(trait_df))
  
  # Handle NULL or "None" rule type — all pass
  if (is.null(rule_def$ruleType) || identical(rule_def$ruleType, "None") || identical(rule_def$ruleType, "")) {
    decision <- rep("SELECTED", nrow(trait_df))
    if (!is.null(check_entry_type_value)) {
      decision[trait_df$entryType == check_entry_type_value] <- "CHECK"
    }
    return(data.frame(
      designation = trait_df$designation,
      decision = decision,
      stringsAsFactors = FALSE
    ))
  }
  
  if (identical(rule_def$ruleType, "Threshold")) {
    if (identical(rule_def$direction, "Higher is better")) {
      decision <- ifelse(trait_value >= rule_def$threshold, "SELECTED", "NOT SELECTED")
    } else if (identical(rule_def$direction, "Lower is better")) {
      decision <- ifelse(trait_value <= rule_def$threshold, "SELECTED", "NOT SELECTED")
    } else {
      stop("Unknown direction for Threshold rule.")
    }
  }
  
  else if (identical(rule_def$ruleType, "Acceptable range")) {
    decision <- ifelse(
      trait_value >= rule_def$minValue & trait_value <= rule_def$maxValue,
      "SELECTED",
      "NOT SELECTED"
    )
  }
  
  else if (identical(rule_def$ruleType, "% over check")) {
    
    ref_check <- rule_def$referenceCheck
    
    check_value <- mean(
      trait_value[trait_df$designation == ref_check],
      na.rm = TRUE
    )
    
    if (!is.finite(check_value)) {
      decision <- rep("NOT SELECTED", nrow(trait_df))
    } else {
      if (identical(rule_def$direction, "Higher is better")) {
        if (!is.null(rule_def$threshold) && rule_def$threshold > 0) {
          cutoff <- check_value + check_value * (rule_def$threshold / 100)
          decision <- ifelse(trait_value >= cutoff, "SELECTED", "NOT SELECTED")
        } else {
          decision <- ifelse(trait_value >= min(trait_value, na.rm = TRUE), "SELECTED", "NOT SELECTED")
        }
      } else if (identical(rule_def$direction, "Lower is better")) {
        if (!is.null(rule_def$threshold) && rule_def$threshold > 0) {
          cutoff <- check_value - check_value * (rule_def$threshold / 100)
          decision <- ifelse(trait_value <= cutoff, "SELECTED", "NOT SELECTED")
        } else {
          decision <- ifelse(trait_value <= max(trait_value, na.rm = TRUE), "SELECTED", "NOT SELECTED")
        }
      } else {
        stop("Unknown direction for % over check rule.")
      }
    }
  }
  
  else if (identical(rule_def$ruleType, "% over mean")) {
    if (is.null(check_entry_type_value)) {
      mean_value <- mean(trait_value, na.rm = TRUE)
    } else {
      mean_value <- mean(
        trait_value[mta_preds_long$entryType != check_entry_type_value],
        na.rm = TRUE
      )
    }
    
    if (!is.finite(mean_value)) {
      decision <- rep("NOT SELECTED", nrow(trait_df))
    } else {
      if (identical(rule_def$direction, "Higher is better")) {
        cutoff <- mean_value + mean_value * (rule_def$threshold / 100)
        decision <- ifelse(trait_value >= cutoff, "SELECTED", "NOT SELECTED")
      } else if (identical(rule_def$direction, "Lower is better")) {
        cutoff <- mean_value - mean_value * (rule_def$threshold / 100)
        decision <- ifelse(trait_value <= cutoff, "SELECTED", "NOT SELECTED")
      } else {
        stop("Unknown direction for % over mean rule.")
      }
    }
  }
  
  else if (is.null(rule_def$ruleType) || identical(rule_def$ruleType, "None") || identical(rule_def$ruleType, "")) {
    # No threshold set — all candidates pass
    decision <- rep("SELECTED", nrow(trait_df))
  }
  
  else {
    stop(paste("Unsupported rule type:", rule_def$ruleType))
  }
  
  if (!is.null(check_entry_type_value)) {
    decision[trait_df$entryType == check_entry_type_value] <- "CHECK"
  }
  
  data.frame(
    designation = trait_df$designation,
    decision = decision,
    stringsAsFactors = FALSE
  )
}


build_prodadv_trait_decision_table <- function(preds, modeling_init, check_entry_type_value = NULL) {
  stopifnot(is.data.frame(preds))
  stopifnot(is.data.frame(modeling_init))
  stopifnot(all(c("designation", "trait", "predictedValue", "entryType") %in% colnames(preds)))
  
  # Exclude traits marked as user_excluded_trait
  excl_traits <- modeling_init$trait[modeling_init$parameter == "user_excluded_trait"]
  modeling_filtered <- modeling_init[!modeling_init$trait %in% excl_traits, , drop = FALSE]
  
  rule_list <- get_prodadv_trait_rules(modeling_filtered)
  selected_traits <- names(rule_list)
  
  decision_tables <- lapply(selected_traits, function(tr) {
    trait_df <- preds[preds$trait == tr, c("designation", "predictedValue", "entryType"), drop = FALSE]
    
    out <- apply_prodadv_trait_rule(
      trait_df = trait_df,
      rule_def = rule_list[[tr]],
      check_entry_type_value = check_entry_type_value
    )
    
    names(out)[names(out) == "decision"] <- paste0(tr, "_trait_decision")
    out
  })
  
  trait_decisions <- Reduce(
    function(x, y) merge(x, y, by = "designation", all = TRUE),
    decision_tables
  )
  
  trait_decisions
}



#Helpers for TPE plot

idw_surface <- function(df_points,
                        value_col = "env_value",
                        lon_col = "LON",
                        lat_col = "LAT",
                        n_grid = 80,
                        power = 2,
                        lon_min = NULL,
                        lon_max = NULL,
                        lat_min = NULL,
                        lat_max = NULL) {
  
  if (!is.data.frame(df_points)) {
    stop("df_points must be a data.frame.")
  }
  
  required_cols <- c(lon_col, lat_col, value_col)
  if (!all(required_cols %in% colnames(df_points))) {
    stop("df_points must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  df_points <- df_points[
    stats::complete.cases(df_points[, required_cols]),
    ,
    drop = FALSE
  ]
  
  if (nrow(df_points) < 2) {
    stop("At least two environments with complete coordinates are needed to interpolate a surface.")
  }
  
  if (is.null(lon_min)) lon_min <- min(df_points[[lon_col]], na.rm = TRUE)
  if (is.null(lon_max)) lon_max <- max(df_points[[lon_col]], na.rm = TRUE)
  if (is.null(lat_min)) lat_min <- min(df_points[[lat_col]], na.rm = TRUE)
  if (is.null(lat_max)) lat_max <- max(df_points[[lat_col]], na.rm = TRUE)
  
  if (!all(is.finite(c(lon_min, lon_max, lat_min, lat_max)))) {
    stop("Invalid coordinate range for interpolation.")
  }
  
  if (lon_min == lon_max) {
    stop("Longitude values are identical; cannot interpolate a 2D surface.")
  }
  
  if (lat_min == lat_max) {
    stop("Latitude values are identical; cannot interpolate a 2D surface.")
  }
  
  lon_seq <- seq(lon_min, lon_max, length.out = n_grid)
  lat_seq <- seq(lat_min, lat_max, length.out = n_grid)
  
  grid <- expand.grid(
    LON = lon_seq,
    LAT = lat_seq
  )
  
  pred_vals <- vapply(seq_len(nrow(grid)), function(i) {
    dx <- df_points[[lon_col]] - grid$LON[i]
    dy <- df_points[[lat_col]] - grid$LAT[i]
    d <- sqrt(dx^2 + dy^2)
    
    if (any(d == 0, na.rm = TRUE)) {
      return(df_points[[value_col]][which.min(d)])
    }
    
    w <- 1 / (d^power)
    sum(w * df_points[[value_col]], na.rm = TRUE) / sum(w, na.rm = TRUE)
  }, numeric(1))
  
  grid$z <- pred_vals
  grid
}


get_tpe_basemap <- function(df_points, buffer_deg = 3, square_expand_factor = 1.05) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for TPE map plotting.")
  }
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    stop("Package 'rnaturalearth' is required for TPE map plotting.")
  }
  
  raw_bbox <- c(
    xmin = min(df_points$LON, na.rm = TRUE) - buffer_deg,
    xmax = max(df_points$LON, na.rm = TRUE) + buffer_deg,
    ymin = min(df_points$LAT, na.rm = TRUE) - buffer_deg,
    ymax = max(df_points$LAT, na.rm = TRUE) + buffer_deg
  )
  
  plot_bbox <- make_squareish_bbox(raw_bbox, expand_factor = square_expand_factor)
  
  bbox_obj <- sf::st_bbox(
    plot_bbox,
    crs = sf::st_crs(4326)
  )
  
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  cropped <- suppressWarnings(sf::st_crop(world, bbox_obj))
  
  list(
    basemap = cropped,
    bbox = bbox_obj
  )
}


sf_to_plotly_polygons <- function(sf_obj) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required.")
  }
  
  if (is.null(sf_obj) || nrow(sf_obj) == 0) {
    return(data.frame())
  }
  
  geom <- sf::st_geometry(sf_obj)
  geom <- sf::st_cast(geom, "MULTIPOLYGON", warn = FALSE)
  
  out_list <- vector("list", length(geom))
  idx <- 1L
  
  for (i in seq_along(geom)) {
    mp <- geom[i][[1]]
    
    poly_list <- vector("list", length(mp))
    
    for (j in seq_along(mp)) {
      coords <- mp[[j]][[1]]
      df <- data.frame(
        LON = coords[, 1],
        LAT = coords[, 2],
        group = paste0("g", i, "_", j),
        stringsAsFactors = FALSE
      )
      poly_list[[j]] <- df
    }
    
    out_list[[idx]] <- do.call(rbind, poly_list)
    idx <- idx + 1L
  }
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}


make_squareish_bbox <- function(bbox, expand_factor = 1.05) {
  xmin <- unname(bbox["xmin"])
  xmax <- unname(bbox["xmax"])
  ymin <- unname(bbox["ymin"])
  ymax <- unname(bbox["ymax"])
  
  xspan <- xmax - xmin
  yspan <- ymax - ymin
  
  if (!all(is.finite(c(xmin, xmax, ymin, ymax, xspan, yspan)))) {
    stop("Invalid bbox values.")
  }
  
  if (xspan <= 0 || yspan <= 0) {
    stop("bbox must have positive width and height.")
  }
  
  target_span <- max(xspan, yspan) * expand_factor
  
  xmid <- (xmin + xmax) / 2
  ymid <- (ymin + ymax) / 2
  
  new_bbox <- c(
    xmin = xmid - target_span / 2,
    xmax = xmid + target_span / 2,
    ymin = ymid - target_span / 2,
    ymax = ymid + target_span / 2
  )
  
  new_bbox
}

build_tpe_mask <- function(df_points, buffer_km = 120) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for interpolation masking.")
  }
  
  pts <- df_points[, c("LON", "LAT"), drop = FALSE]
  pts <- pts[stats::complete.cases(pts), , drop = FALSE]
  
  if (nrow(pts) < 1) {
    return(NULL)
  }
  
  pts_sf <- sf::st_as_sf(pts, coords = c("LON", "LAT"), crs = 4326)
  
  pts_proj <- sf::st_transform(pts_sf, 3857)
  
  mask_proj <- pts_proj |>
    sf::st_buffer(dist = buffer_km * 1000) |>
    sf::st_union()
  
  sf::st_transform(mask_proj, 4326)
}


mask_grid_to_polygon <- function(grid_df, mask_sf) {
  if (is.null(mask_sf)) {
    return(grid_df)
  }
  
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("Package 'sf' is required for interpolation masking.")
  }
  
  if (!all(c("LON", "LAT") %in% colnames(grid_df))) {
    stop("grid_df must contain LON and LAT.")
  }
  
  grid_sf <- sf::st_as_sf(grid_df, coords = c("LON", "LAT"), crs = 4326)
  
  inside <- sf::st_within(grid_sf, mask_sf, sparse = FALSE)[, 1]
  out <- grid_df[inside, , drop = FALSE]
  
  rownames(out) <- NULL
  out
}

buildProdAdvPerformanceProfile <- function(
    mta_long,
    review_df,
    modeling_df,
    performance_profile_scale = c("% over check", "% over mean"),
    plot_selection_overrides = NULL
) {
  performance_profile_scale <- match.arg(performance_profile_scale)
  
  if (!is.data.frame(mta_long)) {
    stop("'mta_long' must be a data.frame.")
  }
  
  if (!is.data.frame(review_df)) {
    stop("'review_df' must be a data.frame.")
  }
  
  if (!is.data.frame(modeling_df)) {
    stop("'modeling_df' must be a data.frame.")
  }
  
  required_mta_cols <- c("designation", "trait", "predictedValue")
  missing_mta_cols <- setdiff(required_mta_cols, colnames(mta_long))
  if (length(missing_mta_cols) > 0) {
    stop(
      "Missing required columns in 'mta_long': ",
      paste(missing_mta_cols, collapse = ", ")
    )
  }
  
  required_review_cols <- c("designation", "plot_status")
  missing_review_cols <- setdiff(required_review_cols, colnames(review_df))
  if (length(missing_review_cols) > 0) {
    stop(
      "Missing required columns in 'review_df': ",
      paste(missing_review_cols, collapse = ", ")
    )
  }
  
  required_modeling_cols <- c("trait", "parameter", "value")
  missing_modeling_cols <- setdiff(required_modeling_cols, colnames(modeling_df))
  if (length(missing_modeling_cols) > 0) {
    stop(
      "Missing required columns in 'modeling_df': ",
      paste(missing_modeling_cols, collapse = ", ")
    )
  }
  
  status_df <- unique(review_df[, c("designation", "plot_status"), drop = FALSE])
  
  df_plot <- merge(
    mta_long,
    status_df,
    by = "designation",
    all.x = TRUE
  )
  
  if (!is.null(plot_selection_overrides)) {
    if (!is.data.frame(plot_selection_overrides)) {
      stop("'plot_selection_overrides' must be NULL or a data.frame.")
    }
    
    if (nrow(plot_selection_overrides) > 0) {
      required_override_cols <- c("designation", "plot_decision")
      missing_override_cols <- setdiff(required_override_cols, colnames(plot_selection_overrides))
      
      if (length(missing_override_cols) > 0) {
        stop(
          "Missing required columns in 'plot_selection_overrides': ",
          paste(missing_override_cols, collapse = ", ")
        )
      }
      
      df_plot <- merge(
        df_plot,
        plot_selection_overrides,
        by = "designation",
        all.x = TRUE
      )
      
      df_plot$plot_status <- ifelse(
        !is.na(df_plot$plot_decision),
        df_plot$plot_decision,
        as.character(df_plot$plot_status)
      )
      
      df_plot$plot_decision <- NULL
    }
  }
  
  df_plot <- df_plot[
    stats::complete.cases(df_plot[, c("designation", "trait", "predictedValue", "plot_status")]),
    ,
    drop = FALSE
  ]
  
  if (nrow(df_plot) == 0) {
    stop("No complete observations available for performance profile plot.")
  }
  
  selected_traits <- unique(modeling_df$trait[!is.na(modeling_df$trait)])
  selected_traits <- selected_traits[nzchar(selected_traits)]
  
  if (length(selected_traits) == 0) {
    stop("No traits found in the modeling table for this selection stamp.")
  }
  
  df_plot <- df_plot[df_plot$trait %in% selected_traits, , drop = FALSE]
  
  candidate_df <- df_plot[df_plot$plot_status != "CHECK", , drop = FALSE]
  checks_df <- df_plot[df_plot$plot_status == "CHECK", , drop = FALSE]
  
  if (nrow(candidate_df) == 0) {
    stop("No non-check candidates available for performance profile plot.")
  }
  
  trait_means <- stats::aggregate(
    predictedValue ~ trait,
    data = candidate_df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  names(trait_means)[names(trait_means) == "predictedValue"] <- "trait_mean"
  
  ref_rows <- modeling_df[
    modeling_df$parameter %in% c("reference_check", "referenceCheck"),
    c("trait", "value"),
    drop = FALSE
  ]
  
  names(ref_rows)[names(ref_rows) == "value"] <- "reference_check"
  
  ref_rows <- ref_rows[
    !is.na(ref_rows$trait) &
      nzchar(ref_rows$trait) &
      !is.na(ref_rows$reference_check) &
      nzchar(ref_rows$reference_check),
    ,
    drop = FALSE
  ]
  
  ref_rows <- ref_rows[!duplicated(ref_rows$trait), , drop = FALSE]
  
  check_ref_list <- lapply(selected_traits, function(tr) {
    trait_checks <- checks_df[checks_df$trait == tr, , drop = FALSE]
    
    if (nrow(trait_checks) == 0) {
      return(data.frame(
        trait = tr,
        check_mean = NA_real_,
        reference_check_used = NA_character_,
        stringsAsFactors = FALSE
      ))
    }
    
    ref_check <- ref_rows$reference_check[match(tr, ref_rows$trait)]
    
    if (length(ref_check) == 1 && !is.na(ref_check) && nzchar(ref_check)) {
      specific_ref <- trait_checks[
        trait_checks$designation == ref_check,
        ,
        drop = FALSE
      ]
      
      if (nrow(specific_ref) > 0) {
        return(data.frame(
          trait = tr,
          check_mean = mean(specific_ref$predictedValue, na.rm = TRUE),
          reference_check_used = ref_check,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    data.frame(
      trait = tr,
      check_mean = mean(trait_checks$predictedValue, na.rm = TRUE),
      reference_check_used = "ALL_CHECKS_MEAN",
      stringsAsFactors = FALSE
    )
  })
  
  check_means <- do.call(rbind, check_ref_list)
  
  candidate_df <- merge(candidate_df, trait_means, by = "trait", all.x = TRUE)
  candidate_df <- merge(candidate_df, check_means, by = "trait", all.x = TRUE)
  
  if (identical(performance_profile_scale, "% over mean")) {
    candidate_df$profile_value <- ifelse(
      is.na(candidate_df$trait_mean) | candidate_df$trait_mean == 0,
      NA_real_,
      100 * (candidate_df$predictedValue - candidate_df$trait_mean) / abs(candidate_df$trait_mean)
    )
    x_lab <- "% over mean"
  } else {
    candidate_df$profile_value <- ifelse(
      is.na(candidate_df$check_mean) | candidate_df$check_mean == 0,
      NA_real_,
      100 * (candidate_df$predictedValue - candidate_df$check_mean) / abs(candidate_df$check_mean)
    )
    x_lab <- "% over check"
  }
  
  candidate_df <- candidate_df[is.finite(candidate_df$profile_value), , drop = FALSE]
  
  if (nrow(candidate_df) == 0) {
    stop(paste("No valid values available for", x_lab, "calculation."))
  }
  
  candidate_df$profile_value_clipped <- pmax(-100, pmin(100, candidate_df$profile_value))
  
  candidate_df <- candidate_df[, c(
    "designation",
    "trait",
    "predictedValue",
    "plot_status",
    "trait_mean",
    "check_mean",
    "reference_check_used",
    "profile_value",
    "profile_value_clipped"
  ), drop = FALSE]
  
  list(
    plot_df = candidate_df,
    x_lab = x_lab
  )
}


#' Detect the relatedness visualization mode based on available data
#'
#' Determines which visualization mode to use for the relatedness plot by
#' checking for the presence of pedigree and genomic data in a phenoDTfile object.
#'
#' @param phenoDTfile A list representing the phenotypic data object, expected to
#'   contain \code{$data$pedigree}, \code{$metadata$pedigree}, and
#'   \code{$data$geno_imp}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{mode}{Character string: one of \code{"pedigree-only"},
#'     \code{"pedigree-genomic"}, \code{"genomic-only"}, or \code{"none"}.}
#'   \item{has_pedigree}{Logical indicating whether valid pedigree data is available.}
#'   \item{has_genomic}{Logical indicating whether genomic data is available.}
#'   \item{message}{Character string with a validation message when mode is
#'     \code{"none"}, or \code{NULL} otherwise.}
#' }
#'
#' @export
detect_relatedness_mode <- function(phenoDTfile) {
  # Check pedigree availability:
  # Pedigree data must be non-null AND metadata$pedigree must have >= 3 mapped columns
  has_pedigree <- FALSE
  if (!is.null(phenoDTfile$data$pedigree)) {
    ped_meta <- phenoDTfile$metadata$pedigree
    if (!is.null(ped_meta) && is.data.frame(ped_meta)) {
      # Check that designation, mother, and father are all mapped
      required_params <- c("designation", "mother", "father")
      mapped_params <- ped_meta$parameter[!is.na(ped_meta$value) & nzchar(trimws(ped_meta$value))]
      if (all(required_params %in% mapped_params)) {
        has_pedigree <- TRUE
      }
    }
  }

  # Check genomic availability:
  # geno_imp must have at least one entry (length >= 1)
  has_genomic <- FALSE
  if (!is.null(phenoDTfile$data$geno_imp)) {
    if (length(phenoDTfile$data$geno_imp) >= 1) {
      has_genomic <- TRUE
    }
  }

  # Determine mode
  if (has_pedigree && has_genomic) {
    mode <- "pedigree-genomic"
  } else if (has_pedigree && !has_genomic) {
    mode <- "pedigree-only"
  } else if (!has_pedigree && has_genomic) {
    mode <- "genomic-only"
  } else {
    mode <- "none"
  }

  # Build message for "none" mode
  msg <- NULL
  if (mode == "none") {
    msg <- "Relatedness plot requires pedigree or genotypic information. Please upload pedigree or marker data in the Data Retrieval section."
  }

  list(
    mode = mode,
    has_pedigree = has_pedigree,
    has_genomic = has_genomic,
    message = msg
  )
}


#' Assign individuals to families based on pedigree cross information
#'
#' Each individual is assigned a family label derived from the unique combination
#' of its mother and father. Individuals with both parents missing are placed in
#' "Founders / Unknown Family". Individuals with one parent missing use
#' "KnownParent \u00d7 Unknown".
#'
#' @param pedigree_df A data.frame containing pedigree records.
#' @param designation_col Character string naming the column with individual identifiers.
#' @param mother_col Character string naming the column with mother identifiers.
#' @param father_col Character string naming the column with father identifiers.
#'
#' @return A named character vector mapping each designation to its family label.
#'
#' @export
assign_families <- function(pedigree_df, designation_col, mother_col, father_col) {
  designations <- as.character(pedigree_df[[designation_col]])
  mothers <- as.character(pedigree_df[[mother_col]])
  fathers <- as.character(pedigree_df[[father_col]])

  # Treat NA and empty string as missing
  mother_missing <- is.na(mothers) | mothers == ""
  father_missing <- is.na(fathers) | fathers == ""

  family_labels <- character(length(designations))

  for (i in seq_along(designations)) {
    m_miss <- mother_missing[i]
    f_miss <- father_missing[i]

    if (m_miss && f_miss) {
      family_labels[i] <- "Founders / Unknown Family"
    } else if (m_miss) {
      family_labels[i] <- paste0(fathers[i], " \u00d7 Unknown")
    } else if (f_miss) {
      family_labels[i] <- paste0(mothers[i], " \u00d7 Unknown")
    } else {
      family_labels[i] <- paste0(mothers[i], " \u00d7 ", fathers[i])
    }
  }

  names(family_labels) <- designations
  family_labels
}


#' Build HTML tooltip for the relatedness plot
#'
#' Constructs an HTML tooltip string for a single individual in the
#' relatedness plot hover display.
#'
#' @param designation Character scalar — individual identifier.
#' @param family_or_cluster Character scalar — family or cluster label.
#' @param plot_status Character scalar — current plot status (e.g. "SELECTED").
#' @param index_value Numeric scalar — Selection Index value (may be NA).
#' @param avg_relatedness Numeric scalar — average relatedness (may be NA).
#' @param is_diversity_candidate Logical scalar — whether this is a diversity candidate.
#'
#' @return A character string containing the HTML tooltip.
#'
#' @export
build_relatedness_tooltip <- function(designation, family_or_cluster, plot_status,
                                      index_value, avg_relatedness, is_diversity_candidate) {
  index_str <- if (is.na(index_value)) "N/A" else sprintf("%.2f", index_value)
  rel_str <- if (is.na(avg_relatedness)) "N/A" else sprintf("%.3f", avg_relatedness)

  tooltip <- paste0(
    "<b>", designation, "</b><br>",
    "Family/Cluster: ", family_or_cluster, "<br>",
    "Status: ", plot_status, "<br>",
    "Selection Index: ", index_str, "<br>",
    "Avg. Relatedness: ", rel_str
  )

  if (isTRUE(is_diversity_candidate)) {
    tooltip <- paste0(tooltip, "<br>\u2605 Diversity candidate")
  }

  tooltip
}















#' Compute summary bar text for the relatedness plot
#'
#' Returns a formatted character string summarizing the selected count, number of
#' groups (families or clusters), number of diversity candidates identified, and
#' the data mode label.
#'
#' @param selected_count Integer scalar — number of selected individuals.
#' @param group_count Integer scalar — number of families or clusters.
#' @param diversity_count Integer scalar — number of diversity candidates identified.
#' @param mode_string Character scalar — the operating mode label
#'   (e.g. "pedigree-only", "pedigree-genomic", "genomic-only").
#'
#' @return A formatted character string for the summary bar display.
#' @export
compute_summary_bar <- function(selected_count, group_count, diversity_count, mode_string) {
  group_label <- if (grepl("genomic-only", mode_string, fixed = TRUE)) {
    "clusters"
  } else {
    "families"
  }

  paste0(
    selected_count, " selected | ",
    group_count, " ", group_label, " | ",
    diversity_count, " diversity candidates | ",
    mode_string, " mode"
  )
}

#' Normalize trait values by direction so higher is always better
#'
#' For traits where the breeding goal is to increase the value, the original
#' values are returned unchanged. For traits where the goal is to decrease,
#' the values are negated so that higher normalized values still represent
#' better performance.
#'
#' @param values Numeric vector of trait values.
#' @param direction Character, either "increase" or "decrease".
#'
#' @return Numeric vector of the same length as \code{values}. Unchanged when
#'   \code{direction} is "increase"; negated when \code{direction} is "decrease".
#'
#' @export
normalize_trait_direction <- function(values, direction) {
  if (direction == "decrease") {
    return(-values)
  }
  values
}

#' Compute Genomic Relationship Matrix (GRM) from a genlight object
#'
#' Extracts the marker dosage matrix from a genlight object, column-mean centers
#' it, and computes the realized genomic relationship matrix as
#' \code{tcrossprod(M_centered) / ncol(M_centered)}.
#'
#' @param genlight_obj A genlight object (from the adegenet package) containing
#'   marker dosage data for a set of individuals.
#'
#' @return A symmetric numeric matrix of dimension n x n (where n is the number
#'   of individuals in the genlight object). Row and column names correspond to
#'   the individual names from the genlight object.
#'
#' @export
compute_grm <- function(genlight_obj) {
  # Extract marker dosage matrix; missing values replaced by column means

  M <- adegenet::tab(genlight_obj, NA.method = "mean")

  # Column-mean center the marker matrix
  col_means <- colMeans(M, na.rm = TRUE)
  M_centered <- sweep(M, 2, col_means, "-")

  # Compute GRM = M_centered %*% t(M_centered) / p, where p = number of markers
  GRM <- tcrossprod(M_centered) / ncol(M_centered)

  return(GRM)
}

#' Compute the Additive Relationship Matrix (A-matrix) from pedigree data
#'
#' Wraps \code{cgiarBase::nrm2()} using pedigree metadata column mappings
#' from the phenoDTfile object. Extracts the designation, mother, and father
#' column names from \code{phenoDTfile$metadata$pedigree} and passes them
#' along with the pedigree data frame to \code{nrm2()}.
#'
#' @param phenoDTfile A list object containing at minimum:
#'   \itemize{
#'     \item \code{data$pedigree}: A data frame with pedigree records.
#'     \item \code{metadata$pedigree}: A data frame with columns \code{parameter}
#'       and \code{value}, mapping "designation", "mother", and "father" to
#'       actual column names in the pedigree data frame.
#'   }
#'
#' @return A symmetric numeric matrix (the numerator relationship matrix) with
#'   row and column names set to the designation identifiers.
#'
#' @export
compute_a_matrix <- function(phenoDTfile) {
  paramsPed <- phenoDTfile$metadata$pedigree
  
  indivCol <- paramsPed[paramsPed$parameter == "designation", "value"]
  damCol <- paramsPed[paramsPed$parameter == "mother", "value"]
  sireCol <- paramsPed[paramsPed$parameter == "father", "value"]
  
  # Guard against missing metadata mappings

  if (length(indivCol) == 0 || length(damCol) == 0 || length(sireCol) == 0) {
    stop("Pedigree metadata mapping incomplete. Required: designation, mother, father.")
  }
  indivCol <- indivCol[1]
  damCol <- damCol[1]
  sireCol <- sireCol[1]

  N <- cgiarBase::nrm2(
    pedData = phenoDTfile$data$pedigree,
    indivCol = indivCol,
    damCol = damCol,
    sireCol = sireCol
  )

  return(N)
}

#' Classify non-selected individuals as diversity candidates
#'
#' A non-selected individual qualifies as a diversity candidate when:
#' \enumerate{
#'   \item Its Selection_Index (\code{index_value}) exceeds the 40th percentile
#'     of index values among the selected set.
#'   \item Its average relatedness to the selected set is below the median of
#'     pairwise relatedness coefficients among selected individuals.
#' }
#'
#' If fewer than 2 individuals are in \code{selected_set}, classification is
#' skipped and an empty character vector is returned. Individuals with missing
#' \code{index_value} or not found in the \code{similarity_matrix} are excluded.
#'
#' @param candidates_df A data.frame with at least columns \code{designation}
#'   (character) and \code{index_value} (numeric, i.e., the Selection_Index).
#' @param similarity_matrix A named symmetric matrix (A-matrix or GRM) with
#'   row and column names corresponding to designations.
#' @param selected_set Character vector of designations that are currently selected.
#'
#' @return Character vector of designations classified as diversity candidates.
#'
#' @export
classify_diversity_candidates <- function(candidates_df, similarity_matrix, selected_set) {
  # Guard: skip classification if fewer than 2 selected individuals

  if (length(selected_set) < 2) {
    return(character(0))
  }

  # Identify non-selected individuals from candidates_df
  non_selected_df <- candidates_df[!(candidates_df$designation %in% selected_set), , drop = FALSE]

  # Compute the 40th percentile of index_value among selected individuals
  selected_df <- candidates_df[candidates_df$designation %in% selected_set, , drop = FALSE]
  selected_index_values <- selected_df$index_value
  # Remove NAs from selected index values for percentile computation

  selected_index_values <- selected_index_values[!is.na(selected_index_values)]
  if (length(selected_index_values) == 0) {
    return(character(0))
  }
  index_threshold <- as.numeric(stats::quantile(selected_index_values, probs = 0.40))

  # Compute median of pairwise relatedness among selected individuals

  # Use only selected individuals present in the similarity matrix
  selected_in_matrix <- intersect(selected_set, rownames(similarity_matrix))
  if (length(selected_in_matrix) < 2) {
    return(character(0))
  }

  selected_submatrix <- similarity_matrix[selected_in_matrix, selected_in_matrix, drop = FALSE]
  # Extract upper triangle (pairwise relatedness among selected)
  pairwise_values <- selected_submatrix[upper.tri(selected_submatrix)]
  relatedness_threshold <- stats::median(pairwise_values, na.rm = TRUE)

  # Evaluate each non-selected individual
  diversity_candidates <- character(0)

  for (i in seq_len(nrow(non_selected_df))) {
    desig <- as.character(non_selected_df$designation[i])
    idx_val <- non_selected_df$index_value[i]

    # Exclude if missing index_value
    if (is.na(idx_val)) {
      next
    }

    # Exclude if not found in similarity_matrix
    if (!(desig %in% rownames(similarity_matrix))) {
      next
    }

    # Check condition (a): index_value > 40th percentile of selected
    if (idx_val <= index_threshold) {
      next
    }

    # Check condition (b): average relatedness to selected set < median pairwise among selected
    relatedness_to_selected <- similarity_matrix[desig, selected_in_matrix]
    avg_relatedness <- mean(relatedness_to_selected, na.rm = TRUE)

    if (avg_relatedness < relatedness_threshold) {
      diversity_candidates <- c(diversity_candidates, desig)
    }
  }

  diversity_candidates
}


#' Find top X unrelated non-selected individuals
#'
#' Identifies top X non-selected individuals meeting relatedness criteria,
#' operating in either pedigree or genomic mode.
#'
#' In pedigree mode: finds individuals from families with zero selected members,
#' sorted by descending index_value.
#'
#' In genomic mode: finds individuals whose average GRM coefficient to the
#' selected set is below the median of pairwise GRM coefficients among selected
#' individuals, sorted by descending index_value.
#'
#' @param candidates_df data.frame with columns: designation, index_value,
#'   and family (family only required for pedigree mode).
#' @param similarity_matrix Named symmetric matrix (A-matrix or GRM) with
#'   row and column names as designations.
#' @param selected_set Character vector of selected designations.
#' @param x Integer, number of individuals to return (0-20).
#' @param mode Character, either "pedigree" or "genomic".
#'
#' @return A data.frame of qualifying individuals ranked by descending
#'   index_value. Contains at minimum columns: designation, index_value
#'   (and family for pedigree mode). Returns empty data.frame if x = 0
#'   or no candidates qualify.
#'
#' @export
find_top_x_unrelated <- function(candidates_df, similarity_matrix, selected_set, x, mode) {


  # If x is 0, return empty data.frame immediately
  if (x == 0) {
    if (mode == "pedigree") {
      return(data.frame(designation = character(0),
                        index_value = numeric(0),
                        family = character(0),
                        stringsAsFactors = FALSE))
    } else {
      return(data.frame(designation = character(0),
                        index_value = numeric(0),
                        stringsAsFactors = FALSE))
    }
  }

  # Identify non-selected individuals
  non_selected_df <- candidates_df[!(candidates_df$designation %in% selected_set), , drop = FALSE]

  if (nrow(non_selected_df) == 0) {
    if (mode == "pedigree") {
      return(data.frame(designation = character(0),
                        index_value = numeric(0),
                        family = character(0),
                        stringsAsFactors = FALSE))
    } else {
      return(data.frame(designation = character(0),
                        index_value = numeric(0),
                        stringsAsFactors = FALSE))
    }
  }

  if (mode == "pedigree") {
    # --- Pedigree mode ---
    # Get families of selected individuals
    selected_df <- candidates_df[candidates_df$designation %in% selected_set, , drop = FALSE]
    selected_families <- unique(selected_df$family)

    # Find non-selected individuals from families NOT represented in the selected set
    qualifying <- non_selected_df[!(non_selected_df$family %in% selected_families), , drop = FALSE]

    # Remove rows with NA index_value
    qualifying <- qualifying[!is.na(qualifying$index_value), , drop = FALSE]

    if (nrow(qualifying) == 0) {
      return(data.frame(designation = character(0),
                        index_value = numeric(0),
                        family = character(0),
                        stringsAsFactors = FALSE))
    }

    # Sort by descending index_value
    qualifying <- qualifying[order(qualifying$index_value, decreasing = TRUE), , drop = FALSE]

    # Return the top X (or fewer if not enough candidates)
    n_return <- min(x, nrow(qualifying))
    result <- qualifying[seq_len(n_return), , drop = FALSE]

    # Ensure at minimum the required columns are present
    result <- result[, intersect(c("designation", "index_value", "family"), colnames(result)), drop = FALSE]
    rownames(result) <- NULL
    return(result)

  } else {
    # --- Genomic mode ---
    # Use only selected individuals present in the similarity matrix
    selected_in_matrix <- intersect(selected_set, rownames(similarity_matrix))

    if (length(selected_in_matrix) < 2) {
      # Cannot compute median pairwise if fewer than 2 selected in matrix
      return(data.frame(designation = character(0),
                        index_value = numeric(0),
                        stringsAsFactors = FALSE))
    }

    # Compute the median of pairwise relatedness among selected (upper triangle)
    selected_submatrix <- similarity_matrix[selected_in_matrix, selected_in_matrix, drop = FALSE]
    pairwise_values <- selected_submatrix[upper.tri(selected_submatrix)]
    median_pairwise <- stats::median(pairwise_values, na.rm = TRUE)

    # Evaluate each non-selected individual present in the similarity matrix
    qualifying_rows <- list()

    for (i in seq_len(nrow(non_selected_df))) {
      desig <- as.character(non_selected_df$designation[i])
      idx_val <- non_selected_df$index_value[i]

      # Skip if missing index_value
      if (is.na(idx_val)) next

      # Skip if not present in similarity matrix
      if (!(desig %in% rownames(similarity_matrix))) next

      # Compute average GRM coefficient to all selected individuals present in matrix
      relatedness_to_selected <- similarity_matrix[desig, selected_in_matrix]
      avg_relatedness <- mean(relatedness_to_selected, na.rm = TRUE)

      # Keep only those with average relatedness < median pairwise among selected
      if (avg_relatedness < median_pairwise) {
        qualifying_rows[[length(qualifying_rows) + 1]] <- data.frame(
          designation = desig,
          index_value = idx_val,
          stringsAsFactors = FALSE
        )
      }
    }

    if (length(qualifying_rows) == 0) {
      return(data.frame(designation = character(0),
                        index_value = numeric(0),
                        stringsAsFactors = FALSE))
    }

    qualifying <- do.call(rbind, qualifying_rows)

    # Sort by descending index_value
    qualifying <- qualifying[order(qualifying$index_value, decreasing = TRUE), , drop = FALSE]

    # Return the top X (or fewer if not enough qualify)
    n_return <- min(x, nrow(qualifying))
    result <- qualifying[seq_len(n_return), , drop = FALSE]
    rownames(result) <- NULL
    return(result)
  }
}


#' Cluster selected individuals using genomic relationship matrix
#'
#' Performs hierarchical clustering on the GRM for genomic-only mode.
#' Subsets the GRM to selected individuals, computes distances as 1 - GRM,
#' clusters with average linkage, and cuts at the median pairwise GRM value.
#' Within each cluster, individuals are ordered by descending Selection_Index.
#'
#' @param grm A named symmetric matrix (genomic relationship matrix).
#' @param selected_designations Character vector of selected individual names.
#' @param index_values Named numeric vector where names are designations and values
#'   are Selection_Index scores.
#'
#' @return A list with:
#'   \describe{
#'     \item{clusters}{Named integer vector where names are designations and values
#'       are cluster IDs, ordered within each cluster by descending index.}
#'     \item{dendrogram}{The hclust object, or NULL if fewer than 2 individuals present.}
#'     \item{missing}{Character vector of selected designations not found in GRM.}
#'   }
#'
#' @export
cluster_genomic_only <- function(grm, selected_designations, index_values) {
  # 1. Subset GRM to only those selected designations present in rownames(grm)
  present <- selected_designations[selected_designations %in% rownames(grm)]
  missing_desig <- selected_designations[!(selected_designations %in% rownames(grm))]

  # 2. If fewer than 2 selected are present in GRM, return early

  if (length(present) < 2) {
    empty_clusters <- integer(0)
    names(empty_clusters) <- character(0)
    return(list(clusters = empty_clusters, dendrogram = NULL, missing = missing_desig))
  }

  # Subset the GRM
  grm_subset <- grm[present, present, drop = FALSE]

  # 3. Convert GRM subset to distance: dist_mat = 1 - grm_subset (clamp negatives to 0)
  dist_values <- 1 - grm_subset
  dist_values[dist_values < 0] <- 0
  dist_mat <- stats::as.dist(dist_values)

  # 4. Hierarchical clustering with average linkage
  hc <- stats::hclust(dist_mat, method = "average")

  # 5. Compute median of pairwise GRM values among selected (upper triangle)
  pairwise_grm <- grm_subset[upper.tri(grm_subset)]
  median_grm <- stats::median(pairwise_grm, na.rm = TRUE)

  # 6. Cut the dendrogram at height = 1 - median_grm (distance space)
  cut_height <- 1 - median_grm
  clusters <- stats::cutree(hc, h = cut_height)

  # 7. Within each cluster, order individuals by descending index_values
  # Build a data.frame for ordering
  cluster_df <- data.frame(
    designation = names(clusters),
    cluster_id = as.integer(clusters),
    stringsAsFactors = FALSE
  )

  # Match index values to designations
  cluster_df$index_val <- index_values[cluster_df$designation]
  # For any missing index value, use -Inf so they sort last

  cluster_df$index_val[is.na(cluster_df$index_val)] <- -Inf

  # Order: first by cluster_id ascending, then by index_val descending
  cluster_df <- cluster_df[order(cluster_df$cluster_id, -cluster_df$index_val), ]

  # Build the named integer vector result (ordered within clusters by descending index)
  ordered_clusters <- cluster_df$cluster_id
  names(ordered_clusters) <- cluster_df$designation

  # 8. Return
  list(clusters = ordered_clusters, dendrogram = hc, missing = missing_desig)
}
