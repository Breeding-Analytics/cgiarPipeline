
singleCrossMat <- function( # single cross matrix function
  object= NULL,
  analysisIdForGenoModifications= NULL,
  hybridBatch=1000,
  allHybrids=FALSE,
  verbose=TRUE,
  separator=":"
){

  build.HMM <- function(M1,M2, custom.hyb=NULL, return.combos.only=FALSE, separator=":"){
    # build hybrid marker matrix

    if(!is.null(custom.hyb)){
      pheno <- custom.hyb
      found <- length(which(colnames(pheno) %in% c("Var1","Var2","hybrid")))
      if(found != 3){
        stop("Column names Var1, Var2, hybrid need to be present when you provide \n       a data table to customize the hybrid genotypes to be build.\n", call. = FALSE)
      }
      return.combos.only=FALSE
    }else{
      a <- rownames(M1)
      b <- rownames(M2)
      pheno <- expand.grid(a,b)
      pheno <- pheno[!duplicated(t(apply(pheno, 1, sort))),]
      pheno$hybrid <- paste(pheno$Var1, pheno$Var2, sep=separator)
    }

    if(!return.combos.only){
      # check that marker matrices are in -1,0,1 format
      checkM1 <- c(length(which(M1 == -1)),length(which(M1 == 1)),length(which(M1 == 2)))
      checkM2 <- c(length(which(M2 == -1)),length(which(M2 == 1)),length(which(M2 == 2)))

      checkM1[which(checkM1 > 0)] <- 1
      checkM2[which(checkM2 > 0)] <- 1

      if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
        cat("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
      }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
        cat("Either -1 or 1 alleles not detected in M1, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
      }

      if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
      }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
        cat("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 1 instead of -1 and 1. We'll fix it.\n")
      }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
        cat("Either -1 or 1 alleles not detected in M2, we assume you have coded homozygotes \n       as 0 and 2 instead of -1 and 1. We'll fix it.\n")
      }

      ## add markers coming from parents M1
      parentsX1 <- unique(na.omit(pheno$Var1))
      if(length(parentsX1) > 1){ # more than one parent (usual)
        Z1 <- model.matrix(~Var1-1,pheno);dim(Z1);
        colnames(Z1) <- gsub("Var1","",colnames(Z1))
      }else{ # testcross case maybe?
        Z1 <- matrix(1,nrow=nrow(pheno));dim(Z1);
        colnames(Z1) <- parentsX1
      }
      M1 <- M1[colnames(Z1),, drop=FALSE]
      #M1[1:4,1:4]; Z1[1:4,1:4];
      ## add markers coming from parents M2
      parentsX2 <- unique(na.omit(pheno$Var2))
      if(length(parentsX2) > 1){ # more than one parent (usual)
        Z2 <- model.matrix(~Var2-1,pheno);dim(Z2);
        colnames(Z2) <- gsub("Var2","",colnames(Z2))
      }else{
        Z2 <- matrix(1,nrow=nrow(pheno));dim(Z2);
        colnames(Z2) <- parentsX2
      }
      M2 <- M2[colnames(Z2),, drop=FALSE]
      #M2[1:4,1:4]; Z2[1:4,1:4];

      ## create the
      # Z3 <- model.matrix(~hybrid-1,pheno);dim(Z3);
      # colnames(Z3) <- gsub("hybrid","",colnames(Z3))
      # hyb.names <- colnames(Z3)[as.vector(apply(Z3,1,function(x){which(x==1)}))] # names of hybrids
      hyb.names <- pheno$hybrid
      ## marker matrix for hybrids one for each parent
      cat(paste("Building hybrid marker matrix for",nrow(Z1),"hybrids\n"))

      # M1 <- as(M1, Class="dgCMatrix")
      # M2 <- as(M2, Class="dgCMatrix")
      # Z1 <- as(Z1, Class="dgCMatrix")
      # Z2 <- as(Z2, Class="dgCMatrix")

      cat("Extracting M1 contribution\n")
      if(all(checkM1 == c(1,1,0))){ # homo markers were coded correctly as -1,1
        Md <- Z1 %*% as.matrix(M1);  # was already converted to -1,1
      }else if(all(checkM1 == c(0,1,0)) | all(checkM1 == c(1,0,0))){ # homo markers were coded as 0 1
        Md <- 2*Z1 %*% as.matrix(M1) - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
      }else if(all(checkM1 == c(0,0,1))){ # homo markers were coded as 0 2
        Md <- Z1 %*% as.matrix(M1) - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
      }

      cat("Extracting M2 contribution\n")
      if(all(checkM2 == c(1,1,0))){ # homo markers were coded correctly as -1,1
        Mf <- Z2 %*% as.matrix(M2);  # was already converted to -1,1
      }else if(all(checkM2 == c(0,1,0)) | all(checkM2 == c(1,0,0))){ # homo markers were coded as 0 1
        Mf <- 2*Z2 %*% as.matrix(M2) - 1;  # 2*Z.dent %*% M.dent - 1   # convert to -1,1
      }else if(all(checkM2 == c(0,0,1))){ # homo markers were coded as 0 2
        Mf <- Z2 %*% as.matrix(M2) - 1;  # Z.dent %*% M.dent - 1   # convert to -1,1
      }

      ## marker matrix coded as additive -1,0,1
      Mdf <- (Md + Mf)*(1/2) # normal marker matrix for the hybrids
      rownames(Mdf) <- hyb.names
      #hist(Mdf)

      ## dominance matrix for hybrids (0,1 coded)
      Delta <- 1/2*(1 - Md * Mf) #performs element wise multiplication = Hadamard product
      rownames(Delta) <- hyb.names
      #hist(Delta)
      cat("Done!!\n")
      return(list(HMM.add=Mdf, HMM.dom=Delta, data.used=pheno))

    }else{
      return(list(HMM.add=NA, HMM.dom=NA, data.used=pheno))
    }
  }

  analysisId <- as.numeric(Sys.time())
  ############################
  # loading the dataset
  if (is.null(object)) stop("No input file specified.")

  ## extract marker matrices and reference alleles
  ped <- object$data$pedigree
  metaPed <- object$metadata$pedigree
  colnames(ped) <- cgiarBase::replaceValues(colnames(ped), Search = metaPed$value, Replace = metaPed$parameter )
  cross <- unique(ped[,c("mother","father","designation")]); colnames(cross) <- c("Var1","Var2","hybrid")

  # apply the QA modifications to the markers
  #Markers <- object$data$geno
  #if(is.null(analysisIdForGenoModifications)){ # user didn't provide a modifications id
  #  if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and you have not provided a modifications table to impute the genotype data. Please go to the 'Markers QA/QC' module before performing your analysis.", call. = FALSE)}
  #}else{ # user provided a modifications Id
  #  modificationsMarkers <- object$modifications$geno
  #  theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdForGenoModifications)
  #  if(length(theresMatch) > 0){ # there's a modification file after matching the Id
  #    modificationsMarkers <- modificationsMarkers[theresMatch,]
  #    Markers <- cgiarBase::applyGenoModifications(M=Markers, modifications=modificationsMarkers)
  # }else{ # there's no match of the modification file
  #    if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and your Id didn't have a match in the modifications table to impute the genotype data.", call. = FALSE)}
  #  }
  #}
  if (class(object$data$geno)[1] == "genlight") {
        qas <- which(names(object$data$geno_imp) == analysisIdForGenoModifications)
        Markers <- as.data.frame(object$data$geno_imp[qas])
  }	  
  # Markers <- Markers[,sample(1:min(c(ncol(Markers), nMarkersRRBLUP)))] # don't use all the markers if goes beyond nK
  # Markers <- Markers - 1 # center markers now

  M1 <- Markers[which(rownames(Markers) %in% unique(ped$mother)),]
  M2 <- Markers[which(rownames(Markers) %in% unique(ped$father)),]

  # M1 <- object$data$geno[which(rownames(object$data$geno) %in% unique(ped$mother)),]
  # M2 <- object$data$geno[which(rownames(object$data$geno) %in% unique(ped$father)),]

  ############################
  ## first get all possible hybrids
  if(allHybrids){ # reduce to only hybrids present in the designation column in the pedigree table
    res1 <- build.HMM(M1, M2,return.combos.only = TRUE)
  }else{
    possible <- apply(cross,1, function(x){
      ifelse(x[1]%in%rownames(M1),1,0) + ifelse(x[2]%in%rownames(M2),1,0)
    })
    possible <- which(possible == 2)
    res1 <- list(data.used=cross[possible,])
  }
  res1$data.used <- res1$data.used[which(!is.na(res1$data.used$Var1) & !is.na(res1$data.used$Var2)),]
  res1$data.used <- res1$data.used[sample(1:nrow(res1$data.used), nrow(res1$data.used)),]

  if(nrow(res1$data.used)>0){ # if there is hybrids to build
    ## build the marker matrix for batches of 1000 hybrids
    batches <- sort(rep(1:1000,min(c(nrow(res1$data.used),hybridBatch))))
    res1$data.used$batch <- batches[1:nrow(res1$data.used)]
    data.usedBatches <- split(res1$data.used, res1$data.used$batch)
    # start the loop
    for(i in 1:length(data.usedBatches)){
      prov <- build.HMM(M1=M1-1, M2=M2-1,
                                custom.hyb = data.usedBatches[[i]],
                                return.combos.only=FALSE,
                                separator=separator
      )
      if(i == 1){
        M <- prov$HMM.add +1
      }else{
        M <- rbind(M,prov$HMM.add+1)
      }
    }
  }else{
    M <- matrix(NA, nrow=0, ncol=ncol(M1)); colnames(M) <- colnames(M1)
  }

  # add missing markers to computed hybrid matrix
  missMarkers <- setdiff(colnames(as.data.frame(object$data$geno)),colnames(M))
  addMarkers <- data.frame(matrix(NA,nrow = nrow(M), ncol = length(missMarkers)),row.names = rownames(M))
  colnames(addMarkers) <- missMarkers
  M <- cbind(M, addMarkers)
  M <- M[order(rownames(M)),]

  #object$data$geno <- rbind(object$data$geno,M)
  # other tables
  newStatus <- data.frame(module="scm", analysisId=analysisId,analysisIdName=NA)
  object$status <- rbind( object$status, newStatus[,colnames(object$status)])
  ## add which data was used as input
  modeling <- data.frame(module="scm",  analysisId=analysisId, trait=c("none"), environment="general",
                         parameter= c("allHybrids"), value= c(allHybrids ))
  if(is.null(object$modeling)){
    object$modeling <-  modeling
  }else{
    object$modeling <- rbind(object$modeling, modeling[, colnames(object$modeling)])
  }

  ##
  object$metrics <- rbind(object$metrics,
                               data.frame(module="scm",analysisId=analysisId, trait= "none", environment="across",
                                          parameter=c("nHybrids"), method=c("(Md + Mf) * (1/2)"),
                                          value=c( nrow(M)  ),
                                          stdError=c(NA)
                               )
  )
  return(object)
}
