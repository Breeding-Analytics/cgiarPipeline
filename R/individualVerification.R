
individualVerification <- function(
  object= NULL,
  analysisIdForGenoModifications= NULL,
  markersToBeUsed=NULL,
  colsForExpecGeno=NULL, 
  ploidy=2,
  sc_filter = NULL,
  het = NULL,
  matchThres = NULL,
  onlyMats=FALSE
){

  analysisId <- as.numeric(Sys.time())
  ############################
  # loading the dataset
  if (is.null(object)) stop("No input file specified.")
  if (is.null(analysisIdForGenoModifications)) stop("No geno clean file specified.")
  
  # get markers
  #Markers <- object$data$geno
  if (class(object$data$geno)[1] == "genlight") {
        qas <- which(names(object$data$geno_imp) == analysisIdForGenoModifications)
        Markers <- as.data.frame(object$data$geno_imp[[qas]])
        colnames(Markers) = adegenet::locNames(object$data$geno_imp[[qas]])
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
  if(length(hybsWithM)== 0){
    stop("None of your progeny genotypes have marker information. Please correct.", call. = FALSE)
  }
  if(length(mothersWithM)==0){
    stop("None of your mother genotypes have marker information. Please correct.", call. = FALSE)
  }
  if(length(fathersWithM) == 0){
    stop("None of your father genotypes have marker information. Please correct.", call. = FALSE)
  }
  if(any(is.na(cross$hybrid))){
    stop("Designation column in the pedigree file cannot have missing data. Please correct", call. = FALSE)
  } 
  
  
  #Check if pedigree data is fully available
  progeny_have_mother = !is.na(cross$Var1)
  progeny_have_father = !is.na(cross$Var1)
  progeny_have_parents = progeny_have_mother & progeny_have_father
  
  #Get failed output due to missing pedigree data
  fail_ped = cross[!progeny_have_parents,"hybrid"]
  
  #Check if marker data is fully available
  parents_have_geno = (cross$Var2 %in% rownames(Markers))&(cross$Var1 %in% rownames(Markers))
  progeny_have_geno = cross$hybrid %in% rownames(Markers)
  
  #Get failed output due to missing geno data
  fail_par = cross[!parents_have_geno,"hybrid"]
  fail_pro = cross[!progeny_have_geno,"hybrid"]
  
  #Create filtered df for crosses to be evaluated
  cross_filt = cross[(progeny_have_parents&parents_have_geno&progeny_have_geno),]
  cross_filt$Var1 = as.factor(cross_filt$Var1)
  cross_filt$Var2 = as.factor(cross_filt$Var2)
  cross_filt$hybrid = as.factor(cross_filt$hybrid)
  
  lev <- rownames(Markers)
  lev_fem <- lev[lev %in% cross_filt$Var1]
  lev_mal <- lev[lev %in% cross_filt$Var2]
  lev_hyb <- lev[lev %in% cross_filt$hybrid]

  # enforce that order in the design matrix
  cross_filt$Var1 <- factor(cross_filt$Var1, levels = lev_fem)
  cross_filt$Var2 <- factor(cross_filt$Var2, levels = lev_mal)
  cross_filt$hybrid <- factor(cross_filt$hybrid, levels = lev_hyb)
  
  ## now build the marker matrices for mother, fathers and progeny
  Markers <- Matrix::Matrix(data.matrix(Markers), sparse = TRUE)
  # mother matrix
  Zfem <- Matrix::sparse.model.matrix(~Var1-1, data=cross_filt)
  colnames(Zfem) <- gsub("Var1","",colnames(Zfem))
  Mfem <- Zfem %*% Markers[rownames(Markers) %in% unique(cross_filt$Var1),]
  # father matrix
  Zmal <- Matrix::sparse.model.matrix(~Var2-1, data=cross_filt)
  colnames(Zmal) <- gsub("Var2","",colnames(Zmal))
  Mmal <- Zmal %*% Markers[rownames(Markers) %in% unique(cross_filt$Var2),]
  # progeny matrix
  Zpro <- Matrix::sparse.model.matrix(~hybrid-1, data=cross_filt)
  colnames(Zpro) <- gsub("hybrid","",colnames(Zpro))
  Mpro <- Zpro %*% Markers[rownames(Markers) %in% unique(cross_filt$hybrid),]
  rownames(Mpro) <- cross_filt$hybrid
  # expected progeny matrix
  Mexpec <- Matrix::Matrix(0, nrow=nrow(cross_filt), ncol=ncol(Markers))
  for(iGeno in colsForExpecGeno){
    if(iGeno == "mother"){J=Mfem}else{ if(iGeno == "father"){J=Mmal}else{J=Mpro} }
    Mexpec <- Mexpec + J
  }
  Mexpec <- Mexpec/length(colsForExpecGeno)
  ## calculate metrics
  #Calculate markers informativeness score
  Mfem = as.matrix(Mfem)
  Mmal = as.matrix(Mmal)
  Mpro = as.matrix(Mpro)
  Mexpec = as.matrix(Mexpec)

  res <- cgiarBase::crossVerification(Mf=Mfem,Mm=Mmal,Mp=Mpro,
                                Mexp=Mexpec,
                                ploidy=ploidy,
                                sc_filter = sc_filter,
                                het=het)
  if(onlyMats){
    return(res)
  }else{
  ###############
  # other tables
  newStatus <- data.frame(module="gVerif", analysisId=analysisId, analysisIdName=NA)
  object$status <- rbind( object$status, newStatus[,colnames(object$status)])
  ## modeling
  modeling <- data.frame(module="gVerif",  analysisId=analysisId, trait=c(rep("none",length(colsForExpecGeno)+length(het)+length(matchThres)+length(markersToBeUsed) ),"inputObject"), environment="general",
                         parameter= c(rep("expectedGenoColumn",length(colsForExpecGeno)),"ParentMaxHetThresh", "UpperMatchProbThres","MidMatchProbThres","LowerMatchProbThres",rep("markerUsed",length(markersToBeUsed)), "analysisId"  ),
                         value= c(colsForExpecGeno,het,matchThres,markersToBeUsed,analysisIdForGenoModifications ))
  if(is.null(object$modeling)){
    object$modeling <-  modeling
  }else{
    object$modeling <- rbind(object$modeling, modeling[, colnames(object$modeling)])
  }
  # metrics
  nMarkers = ncol(Mexpec)
  nInds = nrow(cross)
  monomorphicMarkersN = length(which(apply(Mpro,2,var,na.rm=TRUE)==0))
  polymorphicMarkersN = nMarkers - monomorphicMarkersN
  
  #Number of individuals that received PASS:
  
  #Parental homozygosity check
  homo_Mf <- (Mfem == 0 | Mfem == ploidy)
  n_non_missing <- rowSums(!is.na(Mfem))
  n_homo <- rowSums(homo_Mf, na.rm = TRUE)
  het_pct_fem <- (pmax(n_non_missing - n_homo, 0) / pmax(n_non_missing, 1)) * 100
  
  homo_Mm <- (Mmal == 0 | Mmal == ploidy)
  n_non_missing <- rowSums(!is.na(Mmal))
  n_homo <- rowSums(homo_Mm, na.rm = TRUE)
  het_pct_mal <- (pmax(n_non_missing - n_homo, 0) / pmax(n_non_missing, 1)) * 100
  
  #Filters
  het_filter = (het_pct_fem < het & het_pct_mal < het)
  
  indivMatchedN = length(which(res$metricsInd$probMatch[het_filter] >= matchThres[1]))
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
  
  #Add information on missing pedigree and geno data
  res$metricsInd$HasPed = 1
  res$metricsInd$HasGeno = 1
  res$metricsInd$ParentHasGeno = 1
  
  if(length(fail_ped)!= 0){
    tmp = as.data.frame(matrix(NA, nrow = length(fail_ped), ncol = ncol(res$metricsInd)))
    colnames(tmp) = colnames(res$metricsInd)
    tmp$designation = fail_ped
    tmp$HasPed = 0
    tmp$ParentHasGeno = 0
    if(any(fail_ped %in%fail_pro)){
      tmp[fail_ped %in%fail_pro,"HasGeno"] = 0
      tmp[!fail_ped %in%fail_pro,"HasGeno"] = 1
    }else{
      tmp$HasGeno = 1
    }
    
    tmp$parHetFilter = 0
    
    #Update objects
    res$metricsInd = rbind(res$metricsInd,tmp)
    fail_pro = fail_pro[!fail_pro %in% fail_ped]
    fail_par = fail_par[!fail_par %in% fail_ped]
  }
  
  if(length(fail_par)!= 0){
    tmp = as.data.frame(matrix(NA, nrow = length(fail_par), ncol = ncol(res$metricsInd)))
    colnames(tmp) = colnames(res$metricsInd)
    tmp$designation = fail_par
    tmp$HasPed = 1
    tmp$ParentHasGeno = 0
    if(any(fail_par %in%fail_pro)){
      tmp[fail_par%in%fail_pro,"HasGeno"] = 0
      tmp[!fail_par%in%fail_pro,"HasGeno"] = 1
    }else{
      tmp$HasGeno = 1
    }
    
    tmp$parHetFilter = 0 
    
    #Update objects
    res$metricsInd = rbind(res$metricsInd,tmp)
    fail_pro = fail_pro[!fail_pro %in% fail_par]
  }
  
  if(length(fail_pro)!= 0){
    tmp = as.data.frame(matrix(NA, nrow = length(fail_pro), ncol = ncol(res$metricsInd)))
    colnames(tmp) = colnames(res$metricsInd)
    tmp$designation = fail_pro
    tmp$HasPed = 1
    tmp$ParentHasGeno = 1
    tmp$HasGeno = 0
    tmp$parHetFilter = 1
    
    #Update objects
    res$metricsInd = rbind(res$metricsInd,tmp)
  }
  
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


