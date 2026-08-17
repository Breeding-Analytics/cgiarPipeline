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
  ## vectorized baseOrigin computation (replaces row-by-row apply loop)
  .first_sorted <- function(x) (sort(x, decreasing = FALSE, na.last = TRUE))[1]
  .collapse_unique <- function(x) paste(unique(sort(x, decreasing=FALSE)), collapse=", ")
  boGid    <- aggregate(phenoDTfilePred[,"gid",       drop=FALSE], by=list(designation=phenoDTfilePred[,"designation"]), FUN=.first_sorted)
  boMother <- aggregate(phenoDTfilePred[,"mother",    drop=FALSE], by=list(designation=phenoDTfilePred[,"designation"]), FUN=.first_sorted)
  boFather <- aggregate(phenoDTfilePred[,"father",    drop=FALSE], by=list(designation=phenoDTfilePred[,"designation"]), FUN=.first_sorted)
  boPipe   <- aggregate(phenoDTfilePred[,"pipeline",  drop=FALSE], by=list(designation=phenoDTfilePred[,"designation"]), FUN=.collapse_unique)
  boEntry  <- aggregate(phenoDTfilePred[,"entryType", drop=FALSE], by=list(designation=phenoDTfilePred[,"designation"]), FUN=.collapse_unique)
  boEffect <- aggregate(phenoDTfilePred[,"effectType",drop=FALSE], by=list(designation=phenoDTfilePred[,"designation"]), FUN=.collapse_unique)
  baseOrigin <- data.frame(designation=boGid$designation, gid=boGid$gid,
                           mother=boMother$mother, father=boFather$father,
                           pipeline=boPipe$pipeline, entryType=boEntry$entryType,
                           effectType=boEffect$effectType, stringsAsFactors=FALSE)
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




