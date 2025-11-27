indexDesire <- function(
    phenoDTfile= NULL, # input data structure
    analysisId=NULL, # analysis to be picked from predictions database
    environmentToUse=NULL,
    entryTypeToUse=NULL,
    effectTypeToUse=NULL,
    trait= NULL, # traits to include in the index
    desirev = NULL, # vector of desired values
    scaled=TRUE, # whether predicted values should be scaled or not
    verbose=TRUE # should we print logs or not
){
  ## THIS FUNCTION CALCULATES A SELECTION INDEX FOR A SET OF TRAITS AND A VECTR OF DESIRED CHANGE
  ## IS USED IN THE BANAL APP UNDER THE GENETIC EVALUATION MODULES
  idxAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  moduleInput <- phenoDTfile$status[which(phenoDTfile$status$analysisId %in% analysisId),"module"]
  if(length(moduleInput)==0){stop("The file provided doesn't have the analysisId required.",call. = FALSE)}
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(any(moduleInput  %!in% c("mta","mtaFlex","mtaLmms","mas","mtaAsr") ) ){stop("Index can only be calculated on results from a MET analysis using across environment predictions",call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(length(trait) != length(desirev)){stop("The number of traits and desire values needs to be equal",call. = FALSE)}
  names(desirev) <- trait
  ############################
  # loading the dataset
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(!is.null(phenoDTfile$predictions)){
    if("effectType" %!in% colnames(phenoDTfile$predictions) ){
      phenoDTfile$predictions$effectType <- "general"
    }
  }
  mydata <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId),] # readRDS(file.path(wd,"predictions",paste0(phenoDTfile)))
  mydata <- mydata[which(mydata$trait %in% trait),]
  if(is.null(environmentToUse)){ environmentToUse <- names(sort(table(mydata$environment)))}
  mydata <- mydata[which(mydata$environment %in% environmentToUse),]
  if(!is.null(entryTypeToUse)){
    if(length(setdiff(entryTypeToUse,"")) > 0){
      mydata <- mydata[which(mydata$entryType %in% entryTypeToUse),]
    }
  }
  if(!is.null(effectTypeToUse)){
    if(length(setdiff(effectTypeToUse,"")) > 0){
      mydata <- mydata[which(mydata$effectType %in% effectTypeToUse),]
    }
  }
  trait <- intersect(trait, unique(mydata$trait))
  desirev <- desirev[trait]
  ############################
  # if the user provides two ids with same traits kill the job
  traitByIdCheck <- with(mydata, table(trait, analysisId))
  traitByIdCheck <- traitByIdCheck/traitByIdCheck;
  checkOnPreds <- apply(traitByIdCheck,1,sum, na.rm=TRUE)
  badIdSelection <- which( checkOnPreds > 1)
  if(length(badIdSelection) > 0){
    stop(paste( "You have selected multiple analysisId to be analyzed together but trait(s)",paste(names(checkOnPreds)[badIdSelection], collapse =", "),"has data in multiple files") , call. = FALSE)
  }
  ############################
  ## index calculation
  wide0 <- reshape(mydata[,c("designation","trait","predictedValue")], direction = "wide", idvar = "designation",
                   timevar = "trait", v.names = "predictedValue", sep= "_")
  wide <- as.matrix(wide0[,-1]); colnames(wide) <- gsub("predictedValue_","", colnames(wide0)[-1])#unique(mydata$trait)
  wide <- as.matrix(wide[,trait]); colnames(wide) <- trait # ensure order of the user so weights also match
  wide <- apply(wide,2,enhancer::imputev)
  if(scaled){
    if(verbose){cat(paste("scaled has been set to",scaled,"'desirev' values are expected to be the desired change in std. deviations \n"))}
    wide <- apply(wide,2,scale)
    wide[which(is.na(wide), arr.ind = TRUE)] <- 0
  }else{
    if(verbose){cat(paste("scaled has been set to",scaled,"'desirev' values are expected to be desired change in original units \n")) }
  }
  G <- cov(wide, use="pairwise.complete.obs")
  G[which(is.na(G), arr.ind = TRUE)] <- 0
  b <- MASS::ginv(G)%*%desirev # desired weights Ginv*d, equivalent to knowing w (economic weights)
  merit <- wide %*% b
  newped <- data.frame(analysisId=idxAnalysisId,trait="desireIndex",
                       designation=wide0[,1], predictedValue=merit,stdError=1e-6,reliability=1e-6,
                       environment=paste(environmentToUse, collapse="_") )
  ##########################################
  ## add timePoint of origin, stage and designation code
  entries <- unique(mydata[,"designation"])
  baseOrigin <- do.call(rbind, apply(data.frame(entries),1,function(x){
    out1 <- (sort(mydata[which(mydata$designation %in% x),"gid"], decreasing = FALSE))[1]
    out2 <- (sort(mydata[which(mydata$designation %in% x),"mother"], decreasing = FALSE))[1]
    out3 <- (sort(mydata[which(mydata$designation %in% x),"father"], decreasing = FALSE))[1]
    out4 <- paste(unique(sort(mydata[which(mydata$designation %in% x),"pipeline"], decreasing = FALSE)),collapse=", ")
    out5 <- paste(unique(sort(mydata[which(mydata$designation %in% x),"entryType"], decreasing = FALSE)),collapse=", ")
    out6 <- paste(unique(sort(mydata[which(mydata$designation %in% x),"effectType"], decreasing = FALSE)),collapse=", ")
    y <- data.frame(designation=x,gid=out1,mother=out2,father=out3,pipeline=out4, entryType=out5, effectType=out6)
    return(y)
  }))
  predictionsBind <- merge(newped,baseOrigin, by="designation", all.x=TRUE)
  predictionsBind$module <- "indexD"
  if(length(which(predictionsBind$designation=="."))!=0){predictionsBind=predictionsBind[-which(predictionsBind$designation=="."),]}
  #########################################
  ## update databases
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions, predictionsBind[,colnames(phenoDTfile$predictions)])
  modeling <- data.frame(module="indexD",analysisId=idxAnalysisId, trait=if(scaled){paste0(rep(trait,2),"_scaled")}else{rep(trait,2)},
                         environment=environmentToUse,parameter=c(rep("desire",length(trait)),rep("weight",length(trait))),value=c(desirev, b ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[,colnames(phenoDTfile$modeling)])
  newStatus <- data.frame(module="indexD", analysisId=idxAnalysisId,analysisIdName=NA)
  phenoDTfile$status <- rbind(phenoDTfile$status, newStatus[,colnames(phenoDTfile$status)])
  modeling1 <- data.frame(module="indexD",  analysisId=idxAnalysisId, trait=c("inputObject"), environment="general",
                          parameter= c("analysisId"), value= c(analysisId))
  modeling2 <- data.frame(module="indexD",  analysisId=idxAnalysisId, trait=c("general"), environment="general",
                          parameter= c("scaled", rep("entryTypeToUse",length(entryTypeToUse)) ), value= c(scaled, entryTypeToUse))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling1[, colnames(phenoDTfile$modeling)],modeling2[, colnames(phenoDTfile$modeling)])

  return(phenoDTfile)
}
