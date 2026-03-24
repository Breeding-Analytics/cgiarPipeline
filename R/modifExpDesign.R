modifExpDesign <- function(object, df){
  # object is the data object from bioflow
  # df is a table of environments (rows) by factors (cols)
  analysisId <- as.numeric(Sys.time())
  mods <- as.data.frame(matrix(NA, nrow=0, ncol=6))# data.frame(module="qaRaw",analysisId=analysisId,trait=trait,reason=typeOut,row=toSilence, value=NA);
  colnames(mods) <- c("module","analysisId","trait","reason","row","value")
  traits <- object$metadata$pheno[which(object$metadata$pheno$parameter %in% "trait"),"value"]
  envCol <- object$metadata$pheno[which(object$metadata$pheno$parameter %in% "environment"),"value"]


  mydata <- object$data$pheno
  counter <- 1; myList <- list()
  for(i in 1:ncol(df)){ # for each experimental design factor
    for(j in 1:nrow(df)){ # for each environment
      if(df[j,i] == 0){ # if user wants to delete this information
        rowsToSilence <- which(mydata[,envCol] == rownames(df)[j]) # rows for this environments
        valueToSilence <- mydata[rowsToSilence,object$metadata$pheno[which(object$metadata$pheno$parameter == colnames(df)[i]),"value"]]
        myList[[counter]] <- data.frame(module="qaDesign",analysisId=analysisId,trait=colnames(df)[i],reason="badDesign",row=rowsToSilence, value=valueToSilence);
        counter <- counter + 1
      }
    }
  }
  myMods <- do.call(rbind, myList)
  if(!is.null(object$modifications$pheno)){
    object$modifications$pheno <- rbind(object$modifications$pheno, myMods)
  }else{
    object$modifications$pheno <- myMods
  }
  status <- data.frame(module="qaDesign", analysisId=analysisId, analysisIdName=NA)
  if(!is.null(object$status)){
    object$status <- rbind(object$status, status[,colnames(object$status)])
  }else{
    object$status <- status
  }
  return(object)
}
