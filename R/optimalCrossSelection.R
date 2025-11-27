ocs <- function(
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
  ocsAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the predictions", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the analysisId", call. = FALSE)}
  if(is.null(relDTfile)){stop("Please make sure that you have the marker or pedigree data to calculate one of the following relationship matrices: 'grm', 'nrm', or 'both' ", call. = FALSE)}
  if(nchar(relDTfile)==0){stop("Please make sure that you have the marker or pedigree data to calculate one of the following relationship matrices: 'grm', 'nrm', or 'both' ", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(length(trait) > 1){stop(paste0(" Only one trait can be used for optimal contribution. We suggest using an index."), call. = FALSE)}
  if(length(environment) > 1){stop(paste0(" Only one environment can be used for optimal contribution. We suggest using an across environment value."), call. = FALSE)}


  ############################
  # loading the dataset
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

  if(relDTfile %in% c("both","nrm")){ # we need to calculate NRM
    paramsPed <- phenoDTfile$metadata$pedigree
    if(length(intersect(paramsPed$value, colnames(phenoDTfile$data$pedigree)))  < 3){
      stop("Metadata for pedigree (mapping) and pedigree frame do not match. Please reupload and map your pedigree information.", call. = FALSE)
    }
    N <- cgiarBase::nrm2(pedData=phenoDTfile$data$pedigree,
                         indivCol = paramsPed[paramsPed$parameter=="designation","value"],
                         damCol = paramsPed[paramsPed$parameter=="mother","value"],
                         sireCol = paramsPed[paramsPed$parameter=="father","value"]
    )
  }

  if(relDTfile %in% c("grm","both")){ # we need to calculate GRM

    qas<-which( names(phenoDTfile$data$geno_imp)==analysisIdgeno )

    if (length(qas) == 0)
      stop("Genotype version ", analysisIdgeno, " not found in geno_imp.", call. = FALSE)

    gl <- phenoDTfile$data$geno_imp[[ qas[1] ]]
    ## MARKER KERNEL
    M  <- adegenet::tab(gl, NA.method = "mean")

    if(is.null(M)){stop("Markers are not available for this dataset. OCS requires pedigree or markers to work. Please upload any of these in the data retrieval tabs.", call. = FALSE)}

    if(length(qas) > 0){
      if(length(which(is.na(M))) > 0){M <- apply(M,2,enhancer::imputev)}
    }else{
      missing <- apply(M,2,enhancer::propMissing)
      M <- apply(M[,which(missing < 0.9)],2,enhancer::imputev)
    }

    if(ncol(M) > 5000){ # we remove that many markers if a big snp chip
      A <- sommer::A.mat(M[,sample(1:ncol(M), 5000)])
    }else{ A <- sommer::A.mat(M) };  M <- NULL
    if(relDTfile == "both"){ # only if ssgblup we merge
      A <- sommer::H.mat(N,A, tau=1,  omega=1, tolparinv=1e-6)
    }
  }
  if(relDTfile %in% c("nrm")){ A <- N  }
  myrel <- A

  utraits <- unique(mydata$trait) # traits available
  if (!trait %in% utraits){
    stop(paste0("'", trait, "' is not present in the given dataset or the entryType doesn't correspond."), call. = FALSE)
  }
  mydata <- mydata[which(mydata$trait %in% trait),]
  mydata <- mydata[which(mydata$environment %in% environment),] # make sure only across
  mydata <- mydata[with(mydata, order(-predictedValue)), ]
  mydata <- mydata[1:min(c(nrow(mydata), numberBest)),]

  if(nrow(mydata) == 0){stop("Please check the trait and environment selected since there's no phenotypic data for that combination",call. = "FALSE")}
  # make sure you have same phenotypes and genotypes

  common <- intersect(rownames(myrel), mydata[,"designation"])
  if(length(common) == 0){
    stop("There was no intersection between the IDs in the relationship matrix and the IDs in the phenotypes provided. Please check your input files.",call. = FALSE)
  }
  myrel <- myrel[common,common]
  mydata <- mydata[which(mydata[,"designation"] %in% common),]

  ############################
  ## ocs analysis

  forLoop <- expand.grid(nCross, targetAngle)
  predictionsBindList <- list()
  meanCross <- meanFcross <- meanCrossSe <- meanFcrossSe <- numeric()

  for(iRow in 1:nrow(forLoop)){ # iRow=1

    ebv <- data.frame(mydata[,c("predictedValue")]); rownames(ebv) <- mydata[,"designation"]
    ebv <- data.frame(ebv[rownames(myrel),]); rownames(ebv) <- rownames(myrel)
    crossComb = t(combn(1:nrow(myrel), 2)) # all possible cross combintations
    eMP = (ebv[crossComb[,1],] +  ebv[crossComb[,2],])/2  # expected EBVs of all crosses based on # mean parent EBVs
    K <- as.matrix(myrel)

    # OCS: Determine a crossing plan
    plan = cgiarOcs::selectCrosses(nCross=forLoop[iRow,1], # number of crossed to be identified using OCS
                                   targetAngle=((forLoop[iRow,2])*pi)/180, # 30 degrees in radians
                                   u=eMP, # expected cross mean EBVs
                                   maxRun=maxRun,
                                   G=K)   # GRM
    dim(plan$crossPlan)

    crossPlan <- as.data.frame(plan$crossPlan) # list of crosses to be made already sorted by best
    crossPlan[ ,1] <- rownames(K)[crossPlan[ ,1]]
    crossPlan[ ,2] <- rownames(K)[crossPlan[ ,2]]
    colnames(crossPlan) <- c("Parent1", "Parent2", "OCS.merit")
    eMPsel = (ebv[crossPlan[ ,1],] +     # expected EBVs of selected crosses based on
                ebv[crossPlan[ ,2],])/2  # mean parent EBVs
    inbreeding = diag(K)
    inbreedingSel = (inbreeding[crossPlan[ ,1]] + inbreeding[crossPlan[ ,2]])/2
    treatment <- paste(trait,"~", paste(forLoop[iRow,1],"crosses *",forLoop[iRow,2], "degrees"))
    predictionsBindList[[iRow]] <- data.frame(module="ocs",analysisId=ocsAnalysisId, pipeline= paste(sort(unique(mydata$pipeline)),collapse=", "),
                                              trait=trait, gid=1:nrow(crossPlan), designation=paste(crossPlan[,1],crossPlan[,2], sep=" x "),
                                              mother=crossPlan[,1],father=crossPlan[,2], entryType="predictedCross", effectType="designation",
                                              environment=treatment, predictedValue=eMPsel, stdError=inbreedingSel, reliability=crossPlan[,3]
    )
    # bind modeling for this treatment
    modeling <- data.frame(module="ocs",  analysisId=ocsAnalysisId, trait=trait, environment=treatment,
                           parameter= "ocsFormula", value= treatment
    )
    phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
    # bind metric for this treatment
    metrics <- data.frame(module="ocs",  analysisId=ocsAnalysisId, trait=trait, environment=treatment,
                          parameter= c("meanValue","meanF"),method= "sum/n", value=c(mean(eMPsel),mean(inbreedingSel)),
                          stdError=c(sd(eMPsel)/sqrt(length(eMPsel)) ,  sd(inbreedingSel)/sqrt(length(inbreedingSel)) )  )
    phenoDTfile$metrics <- rbind(phenoDTfile$metrics, metrics[, colnames(phenoDTfile$metrics)])

    if(length(otherTraits) > 0){ # if there's more traits in the file, add the value of the crosses for those traits
      traitPredictions <- list()
      for(iTrait in otherTraits){ # iTrait <- otherTraits[1]

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
        eMPtrait = (ebv2[crossPlan[ ,1],] +  ebv2[crossPlan[ ,2],])/2

        traitPredictions[[iTrait]] <- data.frame(module="ocs",  analysisId=ocsAnalysisId, pipeline= paste(sort(unique(mydata$pipeline)),collapse=", "),
                                                 trait=iTrait, gid=1:nrow(crossPlan), designation=paste(crossPlan[,1],crossPlan[,2], sep=" x "),
                                                 mother=crossPlan[,1],father=crossPlan[,2], entryType="predictedCross", effectType="designation",
                                                 environment=treatment, predictedValue=eMPtrait, stdError=inbreedingSel, reliability=crossPlan[,3]
        )
        metrics <- data.frame(module="ocs",  analysisId=ocsAnalysisId, trait=iTrait, environment=treatment,
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
  newStatus <- data.frame(module="ocs", analysisId=ocsAnalysisId, analysisIdName=NA)
  phenoDTfile$status <- rbind(phenoDTfile$status, newStatus[,colnames(phenoDTfile$status)] )
  ## add which data was used as input
  modeling <- data.frame(module="ocs",  analysisId=ocsAnalysisId, trait="inputObject", environment="general",
                         parameter= "analysisId", value= analysisId)
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)
}

