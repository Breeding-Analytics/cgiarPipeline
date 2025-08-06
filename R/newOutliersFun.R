
newOutliersFun <- function(myObject, trait, outlierCoefOutqPheno, traitLBOutqPheno = NULL, traitUBOutqPheno = NULL){
  
  if(is.null(outlierCoefOutqPheno)){outlierCoefOutqPheno <- NA}
  mydata <- myObject$data$pheno
  ### change column names for mapping
  paramsPheno <- myObject$metadata$pheno
  paramsPheno <- paramsPheno[which(paramsPheno$parameter != "trait"),]
  colnames(mydata) <- cgiarBase::replaceValues(colnames(mydata), Search = paramsPheno$value, Replace = paramsPheno$parameter )
  ###
  mydata$rowindex <- 1:nrow(mydata)
  mydata[, "environment"] <- as.factor(mydata[, "environment"])
  analysisId <- NA
  myoutliers <- myObject$modifications$pheno
  if(!is.null(myoutliers) & !is.null(nrow(myoutliers)) ){ # there's previous outliers for this trait
    traitold=gsub("-residual", "", trait)
    outsItrait <- which(myoutliers$trait %in% traitold)
    myoutliersReduced <- myoutliers[outsItrait,] # outliers for the trait in turn
  }else{
    myoutliersReduced <- data.frame(matrix(nrow=0, ncol=6))
    colnames(myoutliersReduced) <- c("module" ,"analysisId" ,"trait","reason","row" , "value" )
  }
  ## add new outliers
  outList <- list(); counter=1
  for (i in 1:nlevels(mydata[, "environment"])) {
    sampleDT <- mydata[which(mydata[, "environment"] == levels(mydata[, "environment"])[i]), ] # data for the ith environment
    if(!is.na(outlierCoefOutqPheno)){
      outlier <- grDevices::boxplot.stats(x=sampleDT[, trait],coef=outlierCoefOutqPheno )$out
      toSilence <- sampleDT[which(sampleDT[,trait] %in% outlier),"rowindex"]
      typeOut <- rep("outlierIQR",length(toSilence))
    }else{
      toSilence <- numeric()
      typeOut <- character()
    }
    if(!is.null(traitLBOutqPheno) & !is.null(traitUBOutqPheno)) {
      outOfBounds <- which((sampleDT[, trait] < traitLBOutqPheno) | (sampleDT[, trait] > traitUBOutqPheno ) )
      if(length(outOfBounds) > 0){toSilence <- c(toSilence, sampleDT[outOfBounds,"rowindex"]); typeOut <- c(typeOut, rep("outlierIQR",length(outOfBounds))) }
    }
    if(length(toSilence) > 0){
      outList[[counter]] <- data.frame(module="qaRaw",analysisId=analysisId,trait=trait,reason=typeOut,row=toSilence, value=NA);
      counter=counter+1
    }
  }# end for each trial
  if(length(outList) > 0){ # we found outliers
    myoutliersReduced2 <- unique(do.call(rbind, outList))
    if( !is.null(myoutliers) & !is.null(nrow(myoutliers)) ){
      myoutliersReduced2 <- unique(rbind(myoutliersReduced,myoutliersReduced2))
    }
  }else{ # we did not find outliers
    myoutliersReduced2 <- data.frame(module="qaRaw",analysisId=analysisId,trait=trait,reason="none",row=NA, value=NA);
    if(!is.null(myoutliers)){ # if there was already outliers in the data structure
      myoutliersReduced2 <- unique(rbind(myoutliersReduced,myoutliersReduced2))
    }else{}
  }
  ## reactive
  return(myoutliersReduced2)
}