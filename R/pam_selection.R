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


  mta_preds_matrix <- matrix(nrow = length(unique(mta_preds$designation)), ncol = length(args$traitsToEvaluate))
  rownames(mta_preds_matrix) <- unique(mta_preds$designation)
  colnames(mta_preds_matrix) <- args$traitsToEvaluate

  idx <- 1
  for(t in args$traitsToEvaluate){
    mta_preds_matrix[,idx] <- mta_preds[mta_preds$trait == t,"predictedValue"]
    idx <- idx + 1
  }

  #Get index weights (if applicable)
  scaled_traits <- NULL

  if(args$decisionLogic == "Weighted index"){

    #Desired gains index selected
    if(args$weightSource == "Use desired gains weights"){
      if(is.null(args$indexStamp)){
        stop("Desired gains were selected. Please run the selection index module and provide an index analysis stamp")
      }else{
       index_model <- dt_object$modeling
       index_model <- index_model[index_model$analysisId == args$indexStamp,]

       if(any(grepl("_scaled",index_model$trait))){
         index_model$trait <- gsub("_scaled","",index_model$trait)
         scaled_traits <- T
       }else{
         scaled_traits <- F
       }

       index_model <- index_model[index_model$trait %in% args$traitsToEvaluate & index_model$parameter == "weight",]
       if (!all(args$traitsToEvaluate %in% index_model$trait)) {
         stop("Some of the traits selected are not present in the desired-gains index. Weights could not be retrieved")
       }

       index_model <- index_model[match(args$traitsToEvaluate,index_model$trait),]
       index_weights <- as.numeric(index_model$value)
       names(index_weights) <- args$traitsToEvaluate

       if(scaled_traits){
         index_preds <- scale(mta_preds_matrix) %*% index_weights
       }else{
         index_preds <- mta_preds_matrix %*% index_weights
       }

      }
    }

    #Custom weights selected
    if(args$weightSource == "Use custom weights"){
      if(is.null(args$customWeights)){
        stop("Custom weights needs to be provided")
      }

      index_preds <- mta_preds_matrix %*% args$customWeights

    }

    #Apply index selection according to top percentage selected
    mta_preds_matrix <- cbind(mta_preds_matrix,index_preds)
    colnames(mta_preds_matrix) <- c(args$traitsToEvaluate, "index")

    mta_preds_long <- as.data.frame(mta_preds_matrix)
    mta_preds_long$designation <- rownames(mta_preds_long)
    mta_preds_long$entryType <- unique(mta_preds[,c("designation","entryType")])[,"entryType"]

    n_selected <- sum(mta_preds_long$entryType != check_entry_type_value) * (args$topPctSelected / 100)
    n_selected <- as.integer(n_selected)
    n_selected <- max(c(1,n_selected))

    index_selection <- order(mta_preds_long$index[mta_preds_long$entryType != check_entry_type_value], decreasing = TRUE)[1:n_selected]
    index_selection <- mta_preds_long$designation[mta_preds_long$entryType != check_entry_type_value][index_selection]
    mta_preds_long$index_decision = NA
    if(length(index_selection) > 1){
      mta_preds_long$index_decision[mta_preds_long$designation %in% index_selection] <- "SELECTED"
      mta_preds_long$index_decision[!mta_preds_long$designation %in% index_selection] <- "NOT SELECTED"
    }else if (length(index_selection) == 1) {
      mta_preds_long$index_decision[mta_preds_long$designation == index_selection] <- "SELECTED"
      mta_preds_long$index_decision[mta_preds_long$designation != index_selection] <- "NOT SELECTED"
    }else{
      mta_preds_long$index_decision <- "NOT SELECTED"
    }
    
    if(!is.null(check_entry_type_value)){
      mta_preds_long$index_decision[mta_preds_long$entryType == check_entry_type_value] <- "CHECK"
    }
   
  }

  if (!exists("mta_preds_long")) {
  mta_preds_long <- as.data.frame(mta_preds_matrix)
  mta_preds_long$designation <- rownames(mta_preds_long)
  mta_preds_long$entryType <- unique(mta_preds[,c("designation","entryType")])[,"entryType"]
  }


  #Apply trait-specific rules

  trait_decision <- matrix(nrow = length(unique(mta_preds$designation)), ncol = length(args$traitsToEvaluate))
  colnames(trait_decision) <- paste0(args$traitsToEvaluate,"_trait_decision")

  for(t in args$traitsToEvaluate){
    trait_rules <- args$traitRules[[t]]
    trait_value <- mta_preds_long[,t]

    if(trait_rules$ruleType == "Threshold"){
      if(trait_rules$direction == "Higher is better"){
        decision <- trait_value >= trait_rules$threshold
      }else if(trait_rules$direction == "Lower is better"){
        decision <- trait_value <= trait_rules$threshold
      }

    }else if(trait_rules$ruleType == "Acceptable range"){
      decision <- trait_value >= trait_rules$minValue & trait_value <= trait_rules$maxValue
    }else if(trait_rules$ruleType == "% over check"){
      if (is.null(check_entry_type_value)) {
        stop("Checks are required when using a '% over check' trait rule.")
      }
      check_value <- mean(trait_value[mta_preds_long$designation == trait_rules$referenceCheck])
      if(trait_rules$direction == "Higher is better"){
        if(trait_rules$threshold > 0){
          decision <- trait_value >= (check_value + check_value*(trait_rules$threshold/100))
        }else{
          decision <- trait_value >= min(trait_value)
        }

      }else if(trait_rules$direction == "Lower is better"){
        if(trait_rules$threshold > 0){
          decision <- trait_value <= (check_value - check_value*(trait_rules$threshold/100))
        }else{
          decision <- trait_value <= max(trait_value)
        }
      }

    }else if(trait_rules$ruleType == "% over mean"){
      if (is.null(check_entry_type_value)) {
        mean_value <- mean(trait_value, na.rm = TRUE)
      } else {
        mean_value <- mean(
          trait_value[mta_preds_long$entryType != check_entry_type_value],
          na.rm = TRUE
        )
      }
      if(trait_rules$direction == "Higher is better"){
        decision <- trait_value >= (mean_value + mean_value*(trait_rules$threshold/100))
      }else if(trait_rules$direction == "Lower is better"){
        decision <- trait_value <= (mean_value - mean_value*(trait_rules$threshold/100))
      }
    }

    decision <- as.character(decision)
    decision[decision == "TRUE"] <- "SELECTED"
    decision[decision == "FALSE"] <- "NOT SELECTED"
    if(!is.null(check_entry_type_value)){
      decision[mta_preds_long$entryType == check_entry_type_value] <- "CHECK"
    }
    
    trait_decision[,paste0(t,"_trait_decision")] <- decision
  }

  mta_preds_long <- cbind(mta_preds_long,trait_decision)

  #Apply trait decision rule
  if(args$decisionLogic == "All traits must pass"){
    decision_set <- mta_preds_long[,grep("trait_decision",colnames(mta_preds_long))]
    initial_decision <- apply(decision_set, 1, function(x) ifelse(any(x == "NOT SELECTED"),"NOT SELECTED","SELECTED"))
    
    if(!is.null(check_entry_type_value)){
      initial_decision[mta_preds_long$entryType == check_entry_type_value] <- "CHECK"
    }
    
  }else if(args$decisionLogic == "At least N traits pass"){
    decision_set <- mta_preds_long[,grep("trait_decision",colnames(mta_preds_long))]
    count_selected <- apply(decision_set, 1, function(x) sum(x == "SELECTED"))
    initial_decision <- ifelse(count_selected >= args$minTraitsPass, "SELECTED","NOT SELECTED")
    if(!is.null(check_entry_type_value)){
      initial_decision[mta_preds_long$entryType == check_entry_type_value] <- "CHECK"
    }
  }else if(args$decisionLogic == "Weighted index"){
    decision_set <- mta_preds_long[,grep("_decision",colnames(mta_preds_long))]
    initial_decision <- apply(decision_set, 1, function(x) ifelse(any(x == "NOT SELECTED"),"NOT SELECTED","SELECTED"))
    if(!is.null(check_entry_type_value)){
      initial_decision[mta_preds_long$entryType == check_entry_type_value] <- "CHECK"
    }
  }

  mta_preds_long$initial_decision <- initial_decision


  #Create modifications table (new list: selection)
  modifications <- data.frame(module = "Init_prodAdv",
                             analysisId = analysisId,
                             designation = mta_preds_long$designation,
                             reason = "initial_selection",
                             value = mta_preds_long$initial_decision)

  #Create modeling table
  for(t in args$traitsToEvaluate){

    weight_value <- NULL

    if (identical(args$weightSource, "Use desired gains weights")) {
      weight_value <- index_weights[t]
    } else if (identical(args$weightSource, "Use custom weights")) {
      weight_value <- if (is.null(args$customWeights)) NULL else args$customWeights[t]
    }

    arg_values <- c(args$mtaStamp,
                   args$staStamp,
                   args$decisionLogic,
                   args$minTraitsPass,
                   args$weightSource,
                   args$topPctSelected,
                   weight_value,
                   scaled_traits,
                   args$traitRules[[t]]$ruleType,
                   args$traitRules[[t]]$direction,
                   args$traitRules[[t]]$threshold,
                   args$traitRules[[t]]$referenceCheck,
                   args$traitRules[[t]]$minValue,
                   args$traitRules[[t]]$maxValue)

    .null_or_name <- function(x, name) {
      if (is.null(x)) NULL else name
    }

    arg_param <- c(
      "mta_stamp",
      "sta_stamp",
      "decision_logic",
      .null_or_name(args$minTraitsPass, "min_traits_pass"),
      .null_or_name(args$weightSource, "weight_source"),
      .null_or_name(args$topPctSelected, "top_perc_selected"),
      .null_or_name(weight_value, "index_weight"),
      .null_or_name(scaled_traits, "scaled_traits"),
      "trait_rule_type",
      .null_or_name(args$traitRules[[t]]$direction, "direction"),
      .null_or_name(args$traitRules[[t]]$threshold, "threshold"),
      .null_or_name(args$traitRules[[t]]$referenceCheck, "reference_check"),
      .null_or_name(args$traitRules[[t]]$minValue, "min_value"),
      .null_or_name(args$traitRules[[t]]$maxValue, "max_value")
    )

    trait_model <- data.frame(module = "Init_prodAdv",
                             analysisId = analysisId,
                             trait = t,
                             environment = NA,
                             parameter = arg_param,
                             value = arg_values)

    if(t == args$traitsToEvaluate[[1]]){
      modeling <- trait_model
    }else{
      modeling <- rbind(modeling,trait_model)
    }

  }

  #Create status table
  status <- data.frame(module = "Init_prodAdv",
                      analysisId = analysisId,
                      analysisIdName = analysisIdName)

  #Update data object
  if(is.null(dt_object$modifications$selection)){
    dt_object$modifications$selection <- modifications
  }else{
    dt_object$modifications$selection <- rbind(dt_object$modifications$selection,modifications)
  }

  dt_object$modeling <- rbind(dt_object$modeling,modeling)
  dt_object$status <- rbind(dt_object$status,status)

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
  
  base_df$is_check_sort <- ifelse(base_df$initial_decision == "CHECK", 0, 1)
  base_df <- base_df[order(base_df$is_check_sort, base_df$designation), , drop = FALSE]
  
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
  
  rule_list <- get_prodadv_trait_rules(modeling_init)
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














