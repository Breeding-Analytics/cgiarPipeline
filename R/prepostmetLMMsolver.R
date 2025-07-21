premetLMMsolver <- function(phenoDTfile= NULL, fixedTerm= NULL, randomTerm=NULL){
  
  envUsed <- gsub("[^[:alnum:]]", "", unique(phenoDTfile$data$pheno[,phenoDTfile$metadata$pheno[which(phenoDTfile$metadata$pheno$parameter == "environment"),"value"]]))
  traitUsed <- unique(phenoDTfile$metadata$pheno[which( phenoDTfile$metadata$pheno$parameter == "trait"),"value"])
  
  gxeList <- list("gxeCS" = c("environment_designation", "designation_environment",
                              "environment:designation", "designation:environment"),
                  "gxeDMD" = c(paste0("env",envUsed,"_designation"), paste0("designation_env",envUsed),
                               paste0("env",envUsed,":designation"), paste0("designation:env",envUsed)),
                  "gxeFW" = c(paste0("value",traitUsed,"envIndex_designation"),paste0("designation_value",traitUsed,"envIndex"),
                              paste0("value",traitUsed,"envIndex:designation"),paste0("designation:value",traitUsed,"envIndex")))
  
  modelTerms <- c(unlist(lapply(fixedTerm, function(x){paste(x, collapse = ":")})),
                  unlist(lapply(randomTerm, function(x){paste(x, collapse = "_")})))
  
  gxeModelCount <- 0
  gxeModelNum <- 0
  for(i in 1:3){
    if(any(modelTerms %in% gxeList[[i]])){
      gxeModelCount <- gxeModelCount + 1
      gxeModelNum <- i
    }
  }
  
  if(gxeModelCount > 1){
    stop("Please ensure that model terms are not conflicting.", call. = FALSE)
  }
  
  if(gxeModelNum > 0){
    gxeTerms <- modelTerms[which(modelTerms %in% gxeList[[gxeModelNum]])]
  } else{
    gxeTerms <- NULL
  }
  
  return(list("gxeModelNum" = gxeModelNum, "gxeTerms" = gxeTerms))
}

postmetLMMsolver <- function(phenoDTfile= NULL, analysisId=NULL, 
                             gxeModelNum=NULL, gxeTerms=NULL){
  
  envUsed <- gsub("[^[:alnum:]]", "", unique(phenoDTfile$data$pheno[,phenoDTfile$metadata$pheno[which(phenoDTfile$metadata$pheno$parameter == "environment"),"value"]]))
  traitUsed <- gsub("[^[:alnum:]]", "", unique(phenoDTfile$metadata$pheno[which( phenoDTfile$metadata$pheno$parameter == "trait"),"value"]))
  
  pred <- phenoDTfile$predictions
  if(gxeModelNum !=0){
    if(gxeModelNum == 3){
      for (i in 1:length(traitUsed)){
        pred[which(pred$analysisId == analysisId & pred$effectType %in% gxeTerms & pred$trait == traitUsed[i]),"predictedValue"] <- pred[which(pred$analysisId == analysisId & pred$effectType %in% gxeTerms & pred$trait == traitUsed[i]),"predictedValue"] - pred[which(pred$analysisId == analysisId & pred$effectType == "(Intercept)" & pred$trait == traitUsed[i]),"predictedValue"]
      }
    }else {
      if(any(grepl(":designation|designation:", gxeTerms))){
        gxeTermsF <- gxeTerms[which(grepl(":designation|designation:", gxeTerms))]
        for(j in 1:length(envUsed)){
          pred[which(pred$analysisId == analysisId & pred$effectType %in% gxeTermsF & grepl(envUsed[j], pred$designation)),"environment"] <- envUsed[j]
        }
      }
    }
    phenoDTfile$predictions <- pred
  }
  
  return(phenoDTfile)
}