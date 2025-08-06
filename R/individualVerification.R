
individualVerification <- function(
  object= NULL,
  analysisIdForGenoModifications= NULL,
  markersToBeUsed=NULL,
  colsForExpecGeno=NULL, ploidy=2,
  onlyMats=FALSE
){

  analysisId <- as.numeric(Sys.time())
  ############################
  # loading the dataset
  if (is.null(object)) stop("No input file specified.")
  if (is.null(analysisIdForGenoModifications)) stop("No geno clean file specified.")
  if (is.null(colsForExpecGeno)){colsForExpecGeno <- "designation"}

  '%!in%' <- function(x,y)!('%in%'(x,y))
  if("effectType" %!in% colnames(object$predictions) ){
    object$predictions$effectType <- "general"
  }
  # get markers
  #Markers <- object$data$geno
  if (class(object$data$geno)[1] == "genlight") {
        qas <- which(names(object$data$geno_imp) == analysisIdForGenoModifications)
        Markers <- as.data.frame(object$data$geno_imp[qas])
  }	 
  if(is.null(Markers)){stop("This function requires your object to have marker information.", call. = FALSE)}
  # apply modifications
  #if(!is.null(analysisIdForGenoModifications)){ # user didn't provide a modifications id
  #  modificationsMarkers <- object$modifications$geno
  #  theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdForGenoModifications)
  #  if(length(theresMatch) > 0){ # there's a modification file after matching the Id
  #    modificationsMarkers <- modificationsMarkers[theresMatch,]
  #    Markers <- cgiarBase::applyGenoModifications(M=Markers, modifications=modificationsMarkers)
  #  }
  #}
  if (is.null(markersToBeUsed)){markersToBeUsed <- 1:ncol(Markers)}else{markersToBeUsed <- intersect(colnames(Markers),markersToBeUsed)}
  Markers <- Markers[,markersToBeUsed]
  ## extract marker matrices and reference alleles
  ped <- object$data$pedigree
  metaPed <- object$metadata$pedigree
  colnames(ped) <- cgiarBase::replaceValues(colnames(ped), Search = metaPed$value, Replace = metaPed$parameter )
  colsForExpecGeno <- cgiarBase::replaceValues(colsForExpecGeno, Search = metaPed$value, Replace = metaPed$parameter )
  cross <- unique(ped[,c("mother","father","designation")]); colnames(cross) <- c("Var1","Var2","hybrid")
  cross <- cross[which(!duplicated(cross$hybrid)),]

  # controls
  if( any( c(sommer::propMissing(cross$Var1), sommer::propMissing(cross$Var2)) == 1 ) ){
    stop("You have too many fathers or mothers missing in the pedigree. Please correct your file", call. = FALSE)
  }
  mothersWithM <- intersect(unique(c(cross$Var1)), rownames(Markers))
  fathersWithM <- intersect(unique(c(cross$Var2)), rownames(Markers))
  hybsWithM <- intersect(unique(c(cross$hybrid)), rownames(Markers))
  if(any(c(length(mothersWithM), length(fathersWithM), length(hybsWithM) ) == 0) ){
    stop("Either your maternal, paternal or progeny genotypes have no marker information. Please correct.", call. = FALSE)
  }
  ## now build the marker matrices for mother, fathers and progeny
  Markers <- as(Markers, Class = "dgCMatrix")
  bad <- Matrix::Matrix(NA, nrow=1, ncol=ncol(Markers))
  rownames(bad) <- "bad"; colnames(bad) <- colnames(Markers)
  Markers <- rbind(Markers,bad)
  cross$Var1[which(is.na(cross$Var1))]="bad"
  cross$Var2[which(is.na(cross$Var2))]="bad"
  cross$hybrid[which(is.na(cross$hybrid))]="bad"
  # mother matrix
  Zfem <- Matrix::sparse.model.matrix(~Var1-1, data=cross)
  colnames(Zfem) <- gsub("Var1","",colnames(Zfem))
  Mfem <- Zfem[,mothersWithM] %*% Markers[mothersWithM,]
  # father matrix
  Zmal <- Matrix::sparse.model.matrix(~Var2-1, data=cross)
  colnames(Zmal) <- gsub("Var2","",colnames(Zmal))
  Mmal <- Zmal[,fathersWithM] %*% Markers[fathersWithM,]
  # progeny matrix
  Zpro <- Matrix::sparse.model.matrix(~hybrid-1, data=cross)
  colnames(Zpro) <- gsub("hybrid","",colnames(Zpro))
  Mpro <- Zpro[,hybsWithM] %*% Markers[hybsWithM,]
  rownames(Mpro) <- cross$hybrid
  # expected progeny matrix
  Mexpec <- Matrix::Matrix(0, nrow=nrow(cross), ncol=ncol(Markers))
  for(iGeno in colsForExpecGeno){
    if(iGeno == "mother"){J=Mfem}else{ if(iGeno == "father"){J=Mmal}else{J=Mpro} }
    Mexpec <- Mexpec + J
  }
  Mexpec <- Mexpec/length(colsForExpecGeno)
  ## calculate metrics

  res <- cgiarBase::crossVerification(Mf=Mfem,Mm=Mmal,Mp=Mpro,
                                Mexp=Mexpec,
                                ploidy=ploidy)
  if(onlyMats){
    return(res)
  }else{
  ###############
  # other tables
  newStatus <- data.frame(module="gVerif", analysisId=analysisId, analysisIdName=NA)
  object$status <- rbind( object$status, newStatus[,colnames(object$status)])
  ## modeling
  modeling <- data.frame(module="gVerif",  analysisId=analysisId, trait=c(rep("none",length(colsForExpecGeno)+length(markersToBeUsed) ),"inputObject"), environment="general",
                         parameter= c(rep("expectedGenoColumn",length(colsForExpecGeno)), rep("markerUsed",length(markersToBeUsed)), "analysisId"  ),
                         value= c(colsForExpecGeno, markersToBeUsed,analysisIdForGenoModifications ))
  if(is.null(object$modeling)){
    object$modeling <-  modeling
  }else{
    object$modeling <- rbind(object$modeling, modeling[, colnames(object$modeling)])
  }
  # metrics
  nMarkers = ncol(Mexpec)
  nInds = nrow(Mexpec)
  monomorphicMarkersN = length(which(apply(Mpro,2,var,na.rm=TRUE)==0))
  polymorphicMarkersN = nMarkers - monomorphicMarkersN
  indivMatchedN = length(which(res$metricsInd$probMatch == 1))
  indivUnmatchedN = nInds - indivMatchedN
  object$metrics <- rbind(object$metrics,
                               data.frame(module="gVerif",analysisId=analysisId, trait="none", environment="general",
                                          parameter= c("nMarkers","nInds", "monomorphicMarkersN","polymorphicMarkersN","indivMatchedN","indivUnmatchedN"),
                                          method= c("sum","sum","sum","sum","sum","sum"),
                                          value=c(nMarkers, nInds, monomorphicMarkersN,polymorphicMarkersN,indivMatchedN, indivUnmatchedN ),
                                          stdError= NA
                               )
  )
  ##############
  ## predictions
  pp <- reshape(res$metricsInd, idvar = "designation", varying = list(2:ncol(res$metricsInd)),
          v.names = "predictedValue", direction = "long", times = colnames(res$metricsInd)[-c(1)]
          )
  pp <- merge(pp,cross, by.x = "designation", by.y="hybrid", all.x = TRUE)

  preds <- data.frame(module="gVerif",  analysisId=analysisId, pipeline= "unknown",
             trait=pp$time, gid=1:nrow(pp), designation=pp$designation,
             mother=pp$Var1,father=pp$Var2, entryType="test", effectType="designation",
             environment="across", predictedValue=pp$predictedValue, stdError=NA,
             reliability=NA
  )
  if(is.null(object$predictions)){
    object$predictions <-  preds
  }else{
    object$predictions <- rbind(object$predictions, preds[, colnames(object$predictions)])
  }

  return(object)
  }
}


