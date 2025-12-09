# Function to calculate Base Index using economic weights

baseIndex <- function(
  phenoDTfile,       # input data
  analysisId,        # analysis from the predictions data  
  analysisIdName,    # analysis from the predictions data  
  trait,             # traits to include in the index
  weights            # vector of economic weights
) {

  if(is.null(phenoDTfile)){

    stop("Please provide a valid dataset", call. = FALSE)
  }

  idxAnalysisId <- as.numeric(Sys.time())

  weights <- as.numeric(weights)
  '%!in%' <- function(x,y)!('%in%'(x,y))
  #if("effectType" %!in% colnames(phenoDTfile$predictions) ){
  #  phenoDTfile$predictions$effectType <- "general"
  #}
  phenoDTfilePred <- phenoDTfile$predictions 
  phenoDTfilePred <- phenoDTfilePred[phenoDTfilePred$analysisId==analysisId,]
  phenoDTfilePred <- phenoDTfilePred[phenoDTfilePred$effectType=="designation",]
  phenoDTfilePred <- phenoDTfilePred[phenoDTfilePred$trait %in% trait,]

  Wide <- reshape(
    data=phenoDTfilePred[,c("designation","trait","predictedValue")],
    timevar = "trait",
    idvar = "designation",
    direction = "wide"
  )
  colnames(Wide)[-1] <- gsub("predictedValue.","", colnames(Wide)[-1])
  Wide.Mat <- as.matrix(Wide[,trait])

  Wide.Mat <- Wide.Mat[,sort(trait)]
  colnames(Wide.Mat) <- sort(trait)

  Wide.Mat <- apply(Wide.Mat,2,enhancer::imputev)

  WideScaled <- data.frame(scale(Wide.Mat, center = TRUE, scale = TRUE))
  WideScaled[which(is.na(WideScaled), arr.ind = TRUE)] <- 0
  WideScaled <- t(t(WideScaled) * weights)

  baseIndex <- rowSums(WideScaled)

  baseIndex <- data.frame(analysisId=idxAnalysisId,
                          trait="baseIndex",
                          designation=Wide[,1],
                          predictedValue=baseIndex,
                          stdError=1e-6,
                          reliability=1e-6,
                          environment="across"
  )
  ##########################################
  ## add timePoint of origin, stage and designation code
  entries <- unique(phenoDTfilePred[,"designation"])
  baseOrigin <- do.call(rbind, apply(data.frame(entries),1,function(x){
    out1 <- (sort(phenoDTfilePred[which(phenoDTfilePred$designation %in% x),"gid"], decreasing = FALSE))[1]
    out2 <- (sort(phenoDTfilePred[which(phenoDTfilePred$designation %in% x),"mother"], decreasing = FALSE))[1]
    out3 <- (sort(phenoDTfilePred[which(phenoDTfilePred$designation %in% x),"father"], decreasing = FALSE))[1]
    out4 <- paste(unique(sort(phenoDTfilePred[which(phenoDTfilePred$designation %in% x),"pipeline"], decreasing = FALSE)),collapse=", ")
    out5 <- paste(unique(sort(phenoDTfilePred[which(phenoDTfilePred$designation %in% x),"entryType"], decreasing = FALSE)),collapse=", ")
    out6 <- paste(unique(sort(phenoDTfilePred[which(phenoDTfilePred$designation %in% x),"effectType"], decreasing = FALSE)),collapse=", ")
    y <- data.frame(designation=x,gid=out1,mother=out2,father=out3,pipeline=out4, entryType=out5, effectType=out6)
    return(y)
  }))
  predictionsBind <- merge(baseIndex,baseOrigin, by="designation", all.x=TRUE)
  predictionsBind$module <- "indexB"
  if(length(which(predictionsBind$designation=="."))!=0){predictionsBind=predictionsBind[-which(predictionsBind$designation=="."),]}
  #########################################
  ## update databases
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,
                                   predictionsBind[,colnames(phenoDTfile$predictions)]
  )
  modeling <- data.frame(module="indexB",
                         analysisId=idxAnalysisId,
                         trait=paste0(trait,"_scaled"),
                         environment="across",
                         parameter=rep("weight",length(trait)),
                         value=weights
  )

  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[,colnames(phenoDTfile$modeling)])  
  newStatus <- data.frame(module="indexB", analysisId=idxAnalysisId, analysisIdName= analysisIdName)
  phenoDTfile$status <- rbind(phenoDTfile$status, newStatus[, colnames(phenoDTfile$status)] )
  modeling <- data.frame(module="indexB",
                         analysisId=idxAnalysisId,
                         trait="inputObject",
                         environment="general",
                         parameter= "analysisId",
                         value= analysisId)
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)

}




