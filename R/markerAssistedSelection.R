
markerAssistedSelection <- function(
    object= NULL,
    analysisIdForGenoModifications= NULL,
    markersToBeUsed=NULL,
    positiveAlleles=NULL,
    desire=NULL, ploidy=2
){

  analysisId <- as.numeric(Sys.time())
  ############################
  # loading the dataset
  if (is.null(object)) stop("No input file specified.")
  if (is.null(analysisIdForGenoModifications)) stop("No geno clean file specified.")
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(!is.null(object$predictions)){
    if("effectType" %!in% colnames(object$predictions) ){
      object$predictions$effectType <- "general"
    }
  }

  # get markers
  if (class(object$data$geno)[1] == "genlight") {
        qas <- which(names(object$data$geno_imp) == analysisIdForGenoModifications)
        Markers <- as.data.frame(object$data$geno_imp[qas])
		colnames(Markers)<-object[["data"]][["geno_imp"]][[qas]]@loc.names
    }	
  Markers <- Markers[,markersToBeUsed]
  ## extract marker matrices and reference alleles
  names(positiveAlleles) <- markersToBeUsed
  X<- object[["data"]][["geno"]]@loc.all
  names(X)<-object[["data"]][["geno"]]@loc.names
  X<-X[names(X)%in%markersToBeUsed]
  X<-strsplit(X,"/")
  X<-data.frame(do.call(rbind,X))
  colnames(X)<-c("refAllele","altAllele")
    
  positiveAllelesN <- vector(mode = "numeric", length = length(positiveAlleles))
  for(j in 1:nrow(X)){
    positiveAllelesN[j] <- ifelse(X[j,c("refAllele")] == positiveAlleles[rownames(X)[j]], 2, 0 )
  }
  # if user doesn't express a desire we weight by 1 minus the frequency
  if(is.null(desire) ){
    p <- vector(mode="numeric",length = length(positiveAllelesN))
    for(k in 1:ncol(Markers)){
      ttb <- table(0:ploidy)-1 # table of zeros for dosages
      tto <- table(Markers[,k])
      ttb[names(tto)] <- ttb[names(tto)] + tto
      n <- sum(ttb)
      p[k] <- ( ifelse(positiveAllelesN[k]==0, ttb[1],ttb[length(ttb)]) + (ttb[2:(length(ttb)-1)] / factorial(2:(length(ttb)-1)) ) ) / n
    }
    desire <- (1 - p) # distance to fixation
  }else{ # we apply both, frequencies and weights from user
    # if(var(desire)==0){
    p <- vector(mode="numeric",length = length(positiveAllelesN))
    for(k in 1:ncol(Markers)){
      ttb <- table(0:ploidy)-1 # table of zeros for dosages
      tto <- table(Markers[,k])
      ttb[names(tto)] <- ttb[names(tto)] + tto
      n <- sum(ttb)
      p[k] <- ( ifelse(positiveAllelesN[k]==0, ttb[1],ttb[length(ttb)]) + (ttb[2:(length(ttb)-1)] / factorial(2:(length(ttb)-1)) ) ) / n
    }
    desire <- (1 - p)*desire # distance to fixation
    # }
  }

  G <- cov(Markers)
  Gi <- solve(G + diag(1e-6,ncol(G),ncol(G)))

  w <- Gi %*% desire

  merit <- as.matrix(Markers) %*% w

  ###############
  # other tables
  newStatus <- data.frame(module="mas", analysisId=analysisId, analysisIdName=NA)
  object$status <- rbind( object$status, newStatus[,colnames(object$status)])
  ## modeling
  modeling <- data.frame(module="mas",  analysisId=analysisId, trait=c(markersToBeUsed,markersToBeUsed,"inputObject"), environment="(Intercept)",
                         parameter= c(rep("markerUsed",length(markersToBeUsed)),rep("desire", length(markersToBeUsed)), "analysisId"  ),
                         value= c(markersToBeUsed,round(desire,6),ifelse(is.null(analysisIdForGenoModifications),NA,analysisIdForGenoModifications))
  )
  if(is.null(object$modeling)){
    object$modeling <-  modeling
  }else{
    object$modeling <- rbind(object$modeling, modeling[, colnames(object$modeling)])
  }
  ##############
  ## predictions
  preds <- data.frame(module="mas",  analysisId=analysisId, pipeline= "unknown",
                      trait="masMerit", gid=1:nrow(Markers), designation=rownames(Markers),
                      mother=NA,father=NA, entryType="test_entry", effectType="designation",
                      environment="(Intercept)", predictedValue=as.numeric(merit), stdError=NA,
                      reliability=NA
  )
  if(is.null(object$predictions)){
    object$predictions <-  preds
  }else{
    object$predictions <- rbind(object$predictions, preds[, colnames(object$predictions)])
  }

  return(object)

}




