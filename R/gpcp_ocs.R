gpcp <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    analysisId=NULL,
    analysisIdgeno=NULL,
    relDTfile= NULL, # "nrm", "grm", "both"
    trait= NULL, # per trait
    environment="across",
    nCross=20,
    targetAngle=30, # in radians
    verbose=FALSE,
    maxRun=100,
    numberBest=100,
    entryType=NULL,
    effectType=NULL
){
  ## THIS FUNCTION CALCULATES THE OPTIMAL CROSS SELECTION USING A TRAIT AND A RELATIONSHIP MATRIX
  ## IS USED IN THE BANAL APP UNDER THE GENETIC EVALUATION MODULES
  gpcpAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the predictions", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the analysisId", call. = FALSE)}
  if(is.null(relDTfile)){stop("Please make sure that you have the marker data to calculate the grm relationship matrix", call. = FALSE)}
  if(nchar(relDTfile)==0){stop("Please make sure that you have the marker data to calculate the grm relationship matrix", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(length(trait) > 1){stop(paste0(" Only one trait can be used for genomic prediction of cross performance. We suggest using an index."), call. = FALSE)}
  if(length(environment) > 1){stop(paste0(" Only one environment can be used for genomic prediction of cross performance. We suggest using an across environment value."), call. = FALSE)}

  ############################
  # loading the dataset
  if(is.null(which(data()$predictions$effectType=="designationA")) && is.null(which(data()$predictions$effectType=="designationD")) && is.null(which(data()$predictions$effectType=="inbreeding"))){
  #if(is.null(phenoDTfile$GPCP)){
    stop("GPCP is only possible if the MTA analysis is done using model 'Main effects (A+D)' ")

  }else{
    if(verbose){
      message("GPCP slot detected in phenoDTfile.")
    }
  }

  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(!is.null(phenoDTfile$predictions)){
    if("effectType" %!in% colnames(phenoDTfile$predictions) ){
      phenoDTfile$predictions$effectType <- "general"
    }
  }

  mydata <- phenoDTfile$predictions #
  mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  if(!is.null(entryType)){
    mydata <- mydata[which(mydata$entryType %in% entryType),]
  }
  if(!is.null(effectType)){
    mydata <- mydata[which(mydata$effectType %in% effectType),]
  }

  if(var(mydata$predictedValue) == 0){
    stop("No variance in the trait to be used for merit. Please correct.", call. = FALSE)
  }

  if( phenoDTfile$status[phenoDTfile$status$analysisId == analysisId,"module"] == "indexD"){
    otherTraits <- setdiff( unique(phenoDTfile$modeling[phenoDTfile$modeling$analysisId == analysisId, "trait"]), c("inputObject","general"))
    analysisIdOtherTraits <- phenoDTfile$modeling[phenoDTfile$modeling$analysisId == analysisId & phenoDTfile$modeling$trait == "inputObject", "value"]
  }else{
    otherTraits <- setdiff(unique(mydata$trait), trait)
    analysisIdOtherTraits <- analysisId
  }

  if(! is.null(relDTfile)){ # we need to calculate backsolving matrices

    qas<-which( names(phenoDTfile$data$geno_imp)==analysisIdgeno )

    if (length(qas) == 0)
      stop("Genotype version ", analysisIdgeno, " not found in geno_imp.", call. = FALSE)

    gl <- phenoDTfile$data$geno_imp[[ qas[1] ]]
    ## MARKER KERNEL
    M  <- adegenet::tab(gl, NA.method = "mean")

    if(is.null(M)){stop("Markers are not available for this dataset. GPCP requires markers to work", call. = FALSE)}

    if(length(qas) > 0){
      if(length(which(is.na(M))) > 0){M <- apply(M,2,enhancer::imputev)}
    }else{
      missing <- apply(M,2,enhancer::propMissing)
      M <- apply(M[,which(missing < 0.9)],2,enhancer::imputev)
    }


    #Backsolving matrices for marker effects
    #Create helper function
    get_backsolving_mat = function(M){
      M.T <- M %*% t(M)
      M.T = M.T + 0.001*diag(nrow(M.T))
      M.Tinv <- solve(M.T) ## inverse
      M.TTinv <- t(M) %*% M.Tinv # M'%*% (M'M)-
      return(M.TTinv)
    }

    ploidyFactor = max(M)/2
    ploidy = ploidyFactor * 2

    Amat = M - ploidyFactor
    Amat = get_backsolving_mat(Amat)

    if(ploidyFactor == 1){
      Dmat = 1- abs(M)
      Dmat = get_backsolving_mat(Dmat)
    }else{
      MAF <- colMeans(M, na.rm = TRUE) / ploidy
      C_mat <- matrix(choose(ploidy, 2), nrow = nrow(t(M)), ncol = ncol(t(M)))
      Ploidy_mat <- matrix(ploidy, nrow = nrow(t(M)), ncol = ncol(t(M)))

      Q <- (MAF^2 * C_mat) -
        (Ploidy_mat - 1) * MAF * t(M) +
        0.5 * t(M) * (t(M)-1)

      Dmat = t(Q)
      Dmat = get_backsolving_mat(Dmat)
    }
  }

  ## Check traits
  utraits <- unique(mydata$trait) # traits available
  if (!trait %in% utraits){
    stop(paste0("'", trait, "' is not present in the given dataset or the entryType doesn't correspond."), call. = FALSE)
  }
  mydata <- mydata[which(mydata$trait %in% trait),]
  mydata <- mydata[which(mydata$environment %in% environment),] # make sure only across
  mydata <- mydata[with(mydata, order(-predictedValue)), ]
  mydata <- mydata[1:min(c(nrow(mydata), numberBest)),]

  if(nrow(mydata) == 0){stop("Please check the trait and environment selected since there's no phenotypic data for that combination",call. = "FALSE")}


  ##################################################
  #Get marker effects
  #GPCP_list = phenoDTfile$GPCP
  GPCP_list <- list()
  GPCP_list$BlupA<-phenoDTfile$predictions[which(phenoDTfile$predictions$effectType=="designationA"),]
  GPCP_list$BlupD<-phenoDTfile$predictions[which(phenoDTfile$predictions$effectType=="designationD"),]
  GPCP_list$f<-phenoDTfile$predictions[which(phenoDTfile$predictions$effectType=="inbreeding"),]
  uno=phenoDTfile$modeling[phenoDTfile$modeling$module=="indexD",]
  uno=uno[which(uno$environment=="(Intercept)"),]
  GPCP_list$index_weights=uno[which(uno$parameter=="weight"),]

  # make sure you have same blups and genotypes

  common <- intersect(rownames(M), GPCP_list$BlupA$designation)

  if(length(common) == 0){
    stop("There was no intersection between the IDs in the relationship matrix and the IDs in the blups provided by MTA analysis. Please check your input files.",call. = FALSE)
  }

  M <- M[common,]

  GPCP_list$BlupA <- GPCP_list$BlupA[which(GPCP_list$BlupA[,"designation"] %in% common),]
  GPCP_list$BlupD <- GPCP_list$BlupD[which(GPCP_list$BlupD[,"designation"] %in% common),]

  #Index check
  if(!is.null(GPCP_list$index_weights)){
    traits_in_index = unique(GPCP_list$index_weights$trait)

    if(any(grepl("scaled",traits_in_index))){
      traits_in_index = gsub("_scaled","",traits_in_index)
    }

    GPCP_list$BlupA = GPCP_list$BlupA[GPCP_list$BlupA$trait %in% traits_in_index,]
    GPCP_list$BlupD = GPCP_list$BlupD[GPCP_list$BlupD$trait %in% traits_in_index,]
    GPCP_list$f = GPCP_list$f[GPCP_list$f$trait %in% traits_in_index,]
  }else{
    warning("Index not present, all traits will have weight = 1")
    traits_in_mta = unique(GPCP_list$f$trait)
    GPCP_list$index_weights = data.frame(trait = traits_in_mta,
                                         value = rep(1, length(traits_in_mta)))
  }

  traits_in_mta = unique(GPCP_list$f$trait)

  add_eff = list()
  dom_eff = list()
  for(t in 1:length(traits_in_mta)){
    blupA = GPCP_list$BlupA[GPCP_list$BlupA$trait == traits_in_mta[t],]
    blupA = blupA[match(rownames(M),blupA$designation),]
    Amat = Amat[,blupA$designation]

    add_eff[[t]] = as.vector(Amat %*% matrix(blupA$predictedValue))

    blupD = GPCP_list$BlupD[GPCP_list$BlupD$trait == traits_in_mta[t],]
    blupD = blupD[match(rownames(M),blupD$designation),]
    Dmat = Dmat[,blupD$designation]

    fCoef = GPCP_list$f[GPCP_list$f$trait == traits_in_mta[t],]
    fCoef = fCoef$predictedValue / ncol(M)

    dom_eff[[t]] = as.vector(Dmat %*% matrix(blupD$predictedValue)) + fCoef

  }

  w = as.numeric(GPCP_list$index_weights$value)

  #Weighted additive and dominance effects:
  ai = Map('*',add_eff, w)
  di = Map('*',dom_eff, w)

  #Sum of weighted effects:
  ai = Reduce('+',ai)
  di = Reduce('+',di)

  #Free up memory
  Amat <- Dmat <- NULL
  if(ploidyFactor != 1){
    MAF <- C_mat <- Ploidy_mat <- Q <- NULL
  }

  mydata <- mydata[which(mydata[,"designation"] %in% common),]
  M <- M[match(mydata[,"designation"],rownames(M)),]

  #Compute relationship matrix
  K <- sommer::A.mat(M)
  K <- as.matrix(K)

  ############################
  ## gpcp + ocs analysis

  forLoop <- expand.grid(nCross, targetAngle)
  predictionsBindList <- list()
  meanCross <- meanFcross <- meanCrossSe <- meanFcrossSe <- numeric()

  for(iRow in 1:nrow(forLoop)){ # iRow=1

    tgv <- data.frame(mydata[,c("predictedValue")]); rownames(tgv) <- mydata[,"designation"]
    tgv <- data.frame(tgv[rownames(M),]); rownames(tgv) <- rownames(M)


    # GPCP+OCS: Determine a crossing plan
    plan = cgiarOcs:: selectCrossPlan(42, # cycle number, currently irrelevant.
                                      nCross = forLoop[iRow,1], # desired number of crosses-- to deal with incompatibility pull more than needed
                                      M = M, # SNP genotypes: ind in rows, snps in columns
                                      a = ai, # additive Effects
                                      d = di, # dominance Effecs
                                      targetAngle = ((forLoop[iRow,2])*pi)/180, # target angle in radians
                                      ploidy = ploidy, # Ploidy
                                      cores = 4)


    dim(plan$crossPlan)

    crossPlan <- as.data.frame(plan$crossPlan) # list of crosses to be made already sorted by best
    crossPlan[ ,1] <- rownames(K)[crossPlan[ ,1]]
    crossPlan[ ,2] <- rownames(K)[crossPlan[ ,2]]
    colnames(crossPlan) <- c("Parent1", "Parent2", "OCS.merit")
    eMPsel = (tgv[crossPlan[ ,1],] +     # expected TGVs of selected crosses based on
                tgv[crossPlan[ ,2],])/2  # mean parent TGVs

    inbreeding = diag(K)
    inbreedingSel = (inbreeding[crossPlan[ ,1]] + inbreeding[crossPlan[ ,2]])/2
    treatment <- paste(trait,"~", paste(forLoop[iRow,1],"crosses *",forLoop[iRow,2], "degrees"))
    predictionsBindList[[iRow]] <- data.frame(module="gpcp",analysisId=gpcpAnalysisId, pipeline= paste(sort(unique(mydata$pipeline)),collapse=", "),
                                              trait=trait, gid=1:nrow(crossPlan), designation=paste(crossPlan[,1],crossPlan[,2], sep=" x "),
                                              mother=crossPlan[,1],father=crossPlan[,2], entryType="predictedCross", effectType="designation",
                                              environment=treatment, predictedValue=eMPsel, stdError=inbreedingSel, reliability=crossPlan[,3]
    )
    # bind modeling for this treatment
    modeling <- data.frame(module="gpcp",  analysisId=gpcpAnalysisId, trait=trait, environment=treatment,
                           parameter= "gpcpFormula", value= treatment
    )
    phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
    # bind metric for this treatment
    metrics <- data.frame(module="gpcp",  analysisId=gpcpAnalysisId, trait=trait, environment=treatment,
                          parameter= c("meanValue","meanF"),method= "sum/n", value=c(mean(eMPsel),mean(inbreedingSel)),
                          stdError=c(sd(eMPsel)/sqrt(length(eMPsel)) ,  sd(inbreedingSel)/sqrt(length(inbreedingSel)) )  )
    phenoDTfile$metrics <- rbind(phenoDTfile$metrics, metrics[, colnames(phenoDTfile$metrics)])

    if(length(otherTraits) > 0){ # if there's more traits in the file, add the value of the crosses for those traits
      traitPredictions <- list()
      for(iTrait in otherTraits){ # iTrait <- otherTraits[1]

        if(grepl("scaled",iTrait)){
          iTrait = gsub("_scaled","",iTrait)
        }

        provPredictions <- phenoDTfile$predictions
        provPredictions <- provPredictions[which(provPredictions$analysisId %in% analysisIdOtherTraits),]
        if(!is.null(entryType)){
          provPredictions <- provPredictions[which(provPredictions$entryType %in% entryType),]
        }
        if(!is.null(effectType)){
          provPredictions <- provPredictions[which(provPredictions$effectType %in% effectType),]
        }
        provPredictions <- provPredictions[which(provPredictions$trait == iTrait), ]
        provPredictions <- provPredictions[which(provPredictions[,"designation"] %in% common),]
        ebv2 <- data.frame(provPredictions[,c("predictedValue")]); rownames(ebv2) <- provPredictions[,"designation"]
        eMPtrait = (ebv2[crossPlan[ ,1],] +  ebv2[crossPlan[ ,2],])/2  #

        traitPredictions[[iTrait]] <- data.frame(module="gpcp",  analysisId=gpcpAnalysisId, pipeline= paste(sort(unique(mydata$pipeline)),collapse=", "),
                                                 trait=iTrait, gid=1:nrow(crossPlan), designation=paste(crossPlan[,1],crossPlan[,2], sep=" x "),
                                                 mother=crossPlan[,1],father=crossPlan[,2], entryType="predictedCross", effectType="designation",
                                                 environment=treatment, predictedValue=eMPtrait, stdError=inbreedingSel, reliability=crossPlan[,3]
        )
        metrics <- data.frame(module="gpcp",  analysisId=gpcpAnalysisId, trait=iTrait, environment=treatment,
                              parameter= c("meanValue"),method= "sum/n", value=c(mean(eMPtrait)),
                              stdError=c(sd(eMPtrait)/sqrt(length(eMPtrait))   )  )
        phenoDTfile$metrics <- rbind(phenoDTfile$metrics, metrics[, colnames(phenoDTfile$metrics)])

      }
      predictionsBindList[[iRow]] <- rbind(predictionsBindList[[iRow]], do.call(rbind, traitPredictions))
    }
  }

  predictionsBind <- do.call(rbind, predictionsBindList)

  #########################################
  ## update structure
  # setdiff(colnames(predictionsBind), colnames(phenoDTfile$predictions))
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,  predictionsBind[, colnames(phenoDTfile$predictions)])
  newStatus <- data.frame(module="gpcp", analysisId=gpcpAnalysisId, analysisIdName=NA)
  phenoDTfile$status <- rbind(phenoDTfile$status, newStatus[,colnames(phenoDTfile$status)] )
  ## add which data was used as input
  modeling <- data.frame(module="gpcp",  analysisId=gpcpAnalysisId, trait="inputObject", environment="general",
                         parameter= "analysisId", value= analysisId)
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)
}
