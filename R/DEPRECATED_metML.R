metML <- function(
    phenoDTfile= NULL, # analysis to be picked from predictions database
    analysisId=NULL,
    analysisIdForGenoModifications=NULL,
    fixedTerm= c("environment"),  randomTerm=c("designation"),  residualBy=NULL,
    interactionsWithGeno=NULL, envsToInclude=NULL,
    trait= NULL, traitFamily=NULL, useWeights=TRUE,
    heritLB= 0.15,  heritUB= 0.95,
    meanLB=0, meanUB=Inf,
    modelType="blup", # either "blup", "pblup", "gblup", "rrblup"
    nMarkersRRBLUP=1000,
    deregress=FALSE,  nPC=0,
    maxIters=50, batchSizeToPredict=500, tolParInv=1e-4,
    minimumNumberEnvsFW=6,
    verbose=TRUE,
    yourLayers = NULL, 
    epochs_M = 50, batchSizeToTrain=50
){
  ## THIS FUNCTION PERFORMS A MULT TRIAL ANALYSIS USING KERAS AND TENSORFLOW
  mtaAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the phenotype file", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the analysisId to be analyzed", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(traitFamily)){traitFamily <- rep("quasi(link = 'identity', variance = 'constant')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  if(length(fixedTerm) == 0 | is.null(fixedTerm)){fixedTerm <- "1"}
  if(modelType %in% c("gblup","ssgblup","rrblup") & is.null(phenoDTfile$data$geno)){
    stop("Please include marker information in your data structure to fit this model type.", call. = FALSE)
  }
  # if user will use markers read them and apply the modifications in the table
  if(modelType %in% c("gblup","ssgblup","rrblup") & !is.null(phenoDTfile$data$geno)){
    Markers <- phenoDTfile$data$geno
    if(is.null(analysisIdForGenoModifications)){ # user didn't provide a modifications id
      if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and you have not provided a modifications table to impute the genotype data. Please go to the 'Markers QA/QC' module prior to run a gBLUP or rrBLUP model.", call. = FALSE)}
    }else{ # user provided a modifications Id
      modificationsMarkers <- phenoDTfile$modifications$geno
      theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdForGenoModifications)
      if(length(theresMatch) > 0){ # there's a modification file after matching the Id
        modificationsMarkers <- modificationsMarkers[theresMatch,]
        Markers <- cgiarBase::applyGenoModifications(M=Markers, modifications=modificationsMarkers)
      }else{ # there's no match of the modification file
        if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and your Id didn't have a match in the modifications table to impute the genotype data.", call. = FALSE)}
      }
    }
    Markers <- Markers[,sample(1:min(c(ncol(Markers), nMarkersRRBLUP)))] # don't use all the markers if goes beyond 1K
  }
  if(is.null(phenoDTfile$metadata$weather)){ # avoid an error when there is no weather information
    provMet <- as.data.frame(matrix(nrow=0, ncol=4));  colnames(provMet) <- c("environment", "trait", "parameter" ,  "value")
    phenoDTfile$metadata$weather <- provMet
  }
  names(traitFamily) <- trait
  heritLB <- rep(heritLB,length(trait))
  heritUB <- rep(heritUB,length(trait))
  meanLB <- rep(meanLB,length(trait))
  meanUB <- rep(meanUB,length(trait))
  traitOrig <- trait
  common <- intersect(fixedTerm,randomTerm)
  fixedTerm <- setdiff(fixedTerm,common) # make sure fixed and random effects don't overlap
  if(length(fixedTerm) == 0 | is.null(fixedTerm)){fixedTerm <- "1"} # assign the intercept if there's no fixed effects
  # print(fixedTerm)
  if("designation" %in% randomTerm){returnFixedGeno=FALSE}else{interactionsWithGeno <- NULL; returnFixedGeno<- TRUE} # don't allow interactions if genotype is not random from start
  if(!"designation" %in% randomTerm){modelType <- "blue"} # assign blue method if the user provided designation in the fixed part
  ############################
  # loading the dataset
  mydata <- phenoDTfile$predictions #
  if (nrow(mydata) < 2) stop("Not enough data is available to perform a multi trial analysis. Please perform an STA before trying to do an MET.", call. = FALSE)
  mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  # add the other available columns to the dataset
  metaPheno <- phenoDTfile$metadata$pheno[which(phenoDTfile$metadata$pheno$parameter %in% c("environment","year","season","country","location","trial")),]
  otherMetaCols <- unique(phenoDTfile$data$pheno[,metaPheno$value,drop=FALSE])
  colnames(otherMetaCols) <- cgiarBase::replaceValues(Source = colnames(otherMetaCols), Search = metaPheno$value, Replace = metaPheno$parameter )
  otherMetaCols <- otherMetaCols[which(!duplicated(otherMetaCols[,"environment"])),,drop=FALSE] # we do this in case the users didn't define the environment properly
  mydata <- merge(mydata, otherMetaCols, by="environment", all.x = TRUE)
  # some checks after filtering
  if(nrow(mydata)==0){stop("No match for this analysisId. Please correct.", call. = FALSE)}
  if( length(setdiff(setdiff(fixedTerm,"1"),colnames(mydata))) > 0 ){stop(paste("column(s):", paste(setdiff(setdiff(fixedTerm,"1"),colnames(mydata)), collapse = ","),"couldn't be found."), call. = FALSE)}
  if( length(setdiff(setdiff(randomTerm,"1"),colnames(mydata))) > 0 ){stop(paste("column(s):", paste(setdiff(setdiff(randomTerm,"1"),colnames(mydata)), collapse = ","),"couldn't be found."), call. = FALSE)}
  # loading the metrics
  metrics <- phenoDTfile$metrics
  metrics <- metrics[which(metrics$analysisId %in% analysisId),]
  #############################
  # defining the name of the surrogate depending on the model
  if(modelType == "blup"){
    surrogate <- "TGV"
  }else if(modelType == "pblup"){
    surrogate <- "EBV"
  }else if(modelType %in% c("gblup","rrblup")){
    surrogate <- "GEBV"
  }else if(modelType == "ssgblup"){
    surrogate <- "ssGEBV"
  }else{surrogate <- "PV"}
  # if user didn't provide a table for which environments should be included, make it! include all environments as default
  if(is.null(envsToInclude)){
    envsToInclude=  as.data.frame( do.call( rbind, list (with(mydata, table(environment,trait)) ) ) )
    bad <- which(envsToInclude <= 1, arr.ind = TRUE)
    if(nrow(bad) > 0){envsToInclude[bad] = 0}
    envsToInclude[which(envsToInclude > 1, arr.ind = TRUE)] = 1
  }
  # check if traits specified are actually available and remove the ones that are not present
  utraits <- unique(mydata$trait)
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% utraits){
      if(verbose){ cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))}
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  if(length(trait)==0){stop("None of the traits specified are available. Please double check", call. = FALSE)}
  traitToRemove <- character()
  if(!is.null(envsToInclude)){
    for(k in 1:length(trait)){
      if( length(which(envsToInclude[,trait[k]] > 0)) == 0 ){
        if(verbose){ cat(paste0("'", trait[k], "' has not envs selected. It will be removed from trait list \n"))}
        traitToRemove <- c(traitToRemove,trait[k])
      }
    }
    trait <- setdiff(trait,traitToRemove)
  }
  if(length(trait)==0){stop("None of the traits specified are available. Please double check", call. = FALSE)}
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  meanLB <- meanLB[which(traitOrig %in% trait)]
  meanUB <- meanUB[which(traitOrig %in% trait)]
  ##############################
  ## met analysis
  allEnvironments <- na.omit(unique(mydata[,"environment"]))
  predictionsList <- list(); counter=counter2=1
  traitToRemove <- character()
  for(iTrait in trait){ # # iTrait = trait[1]  iTrait="value"
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    # subset data
    mydataSub <- droplevels(mydata[which(mydata$trait == iTrait),])
    # 
    if(!is.null(envsToInclude)){
      goodFieldsUser = rownames(envsToInclude)[which(envsToInclude[,iTrait] > 0)]
    }else{goodFieldsUser <- na.omit(unique(mydataSub[,"environment"]))}
    # remove bad environment based on h2 and r2
    pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("plotH2","H2","meanR2","r2")),]
    goodFields <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > heritLB[counter2]) & (pipeline_metricsSub$value < heritUB[counter2])),"environment"])
    goodFields <- intersect(goodFields, goodFieldsUser)
    mydataSub <- mydataSub[which(mydataSub$environment %in% goodFields),]
    # remove bad environment based on environment means
    pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("mean")),]
    if(nrow(pipeline_metricsSub) > 0){ # second reduction for good environments or a given threshold
      goodFieldsMean <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > meanLB[counter2]) & (pipeline_metricsSub$value < meanUB[counter2])),"environment"])
      goodFields <- intersect(goodFields, goodFieldsMean)
      mydataSub <- mydataSub[which(mydataSub$environment %in% goodFields),]
    }
    if(verbose){print(paste("Fields included:",paste(goodFields,collapse = ",")))}
    ## remove records without marker data if marker effects
    if(modelType == "rrblup"){ # if rrBLUP model we remove records without markers
      mydataSub <- mydataSub[which(mydataSub$designation %in% rownames(Markers)),]
      if(nrow(mydataSub) == 0 & verbose){
        cat("There was no match between markers and phenotypes. Maybe your threshold to discard individuals has lead to remove from marker information all the individuals that are present in the phenotypic dataset.")
      }
      Mtrait <- Markers[which(rownames(Markers) %in% unique(mydataSub$designation)),]
    }
    LGrp <- list();   groupTrait <- NULL # in case of GxE models to store the grouping
    ## next step
    if(nrow(mydataSub) == 0){
      print(paste("There was no predictions to work with in trait",iTrait,". Please look at your H2 boundaries. You may be discarding all envs."))
      traitToRemove <- c(traitToRemove,iTrait)
    }else{
      for(iTerm in setdiff(unique(c(fixedTerm, randomTerm)),"1")){
        mydataSub[,iTerm] <- as.factor(mydataSub[,iTerm]) # move to factor
      }
      mydataSub <- mydataSub[which(mydataSub$designation != ""),] # remove blank designations
      ## build the environmental index
      ei <- aggregate(predictedValue~environment, data=mydataSub,FUN=mean, na.rm=TRUE); colnames(ei)[2] <- "envIndex0"
      ei <- ei[with(ei, order(envIndex0)), ]
      ei$envIndex <- ei$envIndex0 - mean(ei$envIndex0)
      colnames(ei) <- cgiarBase::replaceValues(colnames(ei), Search = "envIndex0", Replace = "value")
      ei$parameter <- iTrait # paste0(iTrait,"-envIndex")
      ei$trait <- "envIndex" # paste0(iTrait,"-envIndex")
      # update the weather metadata
      phenoDTfile$metadata$weather <- rbind(phenoDTfile$metadata$weather,ei[,colnames(phenoDTfile$metadata$weather)])
      toKeep <- rownames(unique(phenoDTfile$metadata$weather[,c("environment","trait","parameter")])) # only keep unique records using rownames (alternatively we could use which(!duplicated(phenoDTfile$metadata$weather[,c("environment","parameter")])))
      phenoDTfile$metadata$weather <- phenoDTfile$metadata$weather[toKeep,]
      ## add metadata from environment(e.g., weather) as new columns of the phenotype dataset in case the user wants to model it
      if(!is.null(phenoDTfile$metadata$weather)){
        metas <- phenoDTfile$metadata$weather;
        # metas$feature <- paste(metas$parameter, metas$trait, sep="_")
        set1 <- which(metas$parameter == iTrait) # set of environmental means for iTrait
        set2 <- which(metas$parameter %in% c("mean","date","coordinate") ) # set of weather means
        metas <- metas[c(set1,set2),]
        metas$feature <- paste(metas$environment, metas$trait, metas$parameter)
        metas <- metas[!duplicated(metas$feature),]
        # metas <- metas[which(metas[,"trait"] == iTrait),]
        metas <- reshape(metas[,c("environment","trait","value")], direction = "wide",
                         idvar = "environment",
                         timevar = "trait", v.names = "value", sep= "_")
        colnames(metas) <- gsub("value_","", colnames(metas))
        metasClass <- unlist(lapply(metas,class))
        numericMetas <- names(metasClass)[which(metasClass %in% c("integer","numeric"))]
        # center variables
        for(iMeta in numericMetas){ # iMeta = numericMetas[9]
          if( ( var(as.vector(metas[,iMeta]), na.rm=TRUE) > 0 ) & (length(which(!is.na(metas[,iMeta]))) > minimumNumberEnvsFW) ) { # only add if there is variation in this environmental covariate and we have enough data points
            metas[,iMeta] <- metas[,iMeta] - mean(metas[,iMeta], na.rm=TRUE) # scale(metas[,iMeta]) #  
          }else{
            metas <- metas[,-which(colnames(metas) == iMeta), drop=FALSE]
            print(paste(iMeta, "removed from trait", iTrait, "because doesn't met conditions of variance or minimum number of environments."))
          }
        }
        metas <- metas[which(metas$environment %in% goodFields),, drop=FALSE]
        # colnames(metas) <- cgiarBase::replaceValues(colnames(metas), Search = paste0("envIndex_", iTrait), Replace = "envIndex")
        if(ncol(metas) > 1){ # if we kept some variables merge
          mydataSub <- merge(mydataSub,metas, by="environment", all.x = TRUE) # [,c("environment","envIndex")]
        }
      }
      ## define the interactions to be used
      if(!is.null(interactionsWithGeno)){
        interactionsWithGenoTrait <- interactionsWithGeno
        interactionsWithGenoToRemove <- character()
        for(iInter in 1:length(interactionsWithGenoTrait)){
          if(interactionsWithGenoTrait[iInter] %in% colnames(mydataSub)){ # if trait is even there in dataset
            checkInters <- length(unique(mydataSub[,interactionsWithGenoTrait[iInter]]))
            if (checkInters < minimumNumberEnvsFW){ # there needs to be at least more than one level
              interactionsWithGenoToRemove <- c(interactionsWithGenoToRemove,interactionsWithGenoTrait[iInter])
              print(paste(interactionsWithGenoTrait[iInter], "removed for trait", iTrait))
            }
          }else{ # remove straight away
            interactionsWithGenoToRemove <- c(interactionsWithGenoToRemove,interactionsWithGenoTrait[iInter])
          }
        }
        interactionsWithGenoTrait <- setdiff(interactionsWithGenoTrait,interactionsWithGenoToRemove)
      }else{
        interactionsWithGenoTrait <- interactionsWithGeno
      }
      ##############
      # do analysis
      if(!is.na(var(mydataSub[,"predictedValue"],na.rm=TRUE))){ # if there's variance
        if( var(mydataSub[,"predictedValue"], na.rm = TRUE) > 0 ){
          Ve <- var(mydataSub[,"predictedValue"], na.rm = TRUE)
          # make sure the terms to be fitted have more than one level
          if(deregress){
            mydataSub$predictedValue <- mydataSub$predictedValue/mydataSub$reliability
          }
          if(!is.null(randomTerm)){
            if(modelType == "rrblup"){
              reduced <- with(mydataSub,enhancer::redmm(x=designation,M=Mtrait, nPC=nPC, returnLam = TRUE))
              LGrp[["QTL"]] <- c((ncol(mydataSub)+1):(ncol(mydataSub)+ncol(reduced$Z)))
              mydataSub <- cbind(mydataSub,reduced$Z)
              rTermsTrait <- randomTerm[which(apply(data.frame(randomTerm),1,function(x){length(table(mydataSub[,x]))}) > 1)]
              rTermsTrait <- setdiff(rTermsTrait,"designation")
              rTermsTrait <- c("grp(QTL)",rTermsTrait)
            }else{
              rTermsTrait <- randomTerm[which(apply(data.frame(randomTerm),1,function(x){length(table(mydataSub[,x]))}) > 1)]
            }
          }else{rTermsTrait=NULL}
          
          if(!is.null(fixedTerm)){
            fixedTermTraitMinus <- setdiff(fixedTerm,"1") # remove intercept for a minute
            if(length(fixedTermTraitMinus) > 0){ # check that fixed effect terms have more than one level and can be fitted
              fixedTermTrait <- fixedTermTraitMinus[which(apply(data.frame(fixedTermTraitMinus),1,function(x){length(table(mydataSub[,x]))}) > 1)]
              fixedTermTrait <- c("1",fixedTermTrait)  # add it back
            }else{
              fixedTermTrait <- fixedTerm
            }
          }else{fixedTermTrait=NULL}
          
          if(!is.null(residualBy)){
            residualByTrait <- residualBy[which(apply(data.frame(residualBy),1,function(x){length(table(mydataSub[,x]))}) > 1)]
            if(length(residualByTrait) == 0){residualByTrait=NULL}
          }else{residualByTrait=NULL}
          #####
          if(length(interactionsWithGenoTrait) > 0){ # If there's interactions to be fitted build the formula terms
            if(modelType == "rrblup"){
              for(iInteraction in unique(interactionsWithGenoTrait)){ # iInteraction <- unique(interactionsWithGenoTrait)[1]
                LGrp[[paste0("QTL",iInteraction)]] <- c((ncol(mydataSub)+1):(ncol(mydataSub)+ncol(reduced$Z)))
                mydataSub <- cbind(mydataSub, apply(reduced$Z,2,function(z){z*mydataSub[,iInteraction]}) ) # reduced$Z*mydataSub[,iInteraction])
                rTermsTrait <- c(rTermsTrait,paste0("grp(QTL",iInteraction,")"))
              }
            }else{
              interacs <- expand.grid("designation",interactionsWithGenoTrait)
              interacs<- as.data.frame(interacs[which(as.character(interacs[,1]) != as.character(interacs[,2])),])
              interacsUnlist <- apply(interacs,1,function(x){paste(x,collapse = ":")})
              rTermsTrait <- c(rTermsTrait,interacsUnlist)
            }
          }else{
            rTermsTrait <- rTermsTrait
          }
          rTermsTrait <- setdiff(rTermsTrait, fixedTermTrait)
          
          if(length(rTermsTrait) == 0){ # if there's not a single random effect to be fitted
            ranran <- NULL
            myGinverse <- NULL # if no random effects we don't have a g inverse for LMMsolver
          }else{
            ranran <- paste("~",paste(rTermsTrait, collapse=" + ") )
          }
          fix <- paste("predictedValue ~",paste(fixedTermTrait, collapse=" + "))
          # final residual formula
          if(!is.null(residualByTrait)){ # if the user wants to fit residuals at the level of environment
            ranres <- as.formula(paste0("~",residualByTrait,""))
          }else{
            ranres <- NULL
          }
          
          mydataSub=mydataSub[with(mydataSub, order(environment)), ] # sort by environments
          mydataSub$w <- 1/(mydataSub$stdError^2) # add weights column
          if(verbose){
            cat(fix,"\n")
            cat(ranran,"\n")
          }
          
          if( modelType == "blup" ){ # if user doesn't have any special model
            # make sure the matrix only uses the leves for individuals with data
            designationFlevels <- unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"designation"])
            Ainv <- A <- diag(length(designationFlevels))
            colnames(Ainv) <- rownames(Ainv) <- designationFlevels
            colnames(A) <- rownames(A) <- designationFlevels
            inter <- character()
            onlyInA <- character() # genotypes only present in A and not in dataset
            differ <- character()
            myGinverse <- list(designation=Ainv)
            levelsInAinv <- colnames(Ainv)
            genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=designationFlevels, withMarkOnly=onlyInA)
          }else{
            if(modelType == "rrblup"){ # if user will use the rrblup
              designationFlevels <- unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"designation"])
              inter <- character()
              onlyInA <- character() # genotypes only present in A and not in dataset
              differ <- character()
              myGinverse <- NULL
              groupTrait <- LGrp
              genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=designationFlevels, withMarkOnly=onlyInA)
            }else if(modelType %in% "gblup"){ # if user wants to do a gblup, pblup or ssgblup model
              designationFlevels <- as.character(unique(mydataSub[which(!is.na(mydataSub[,"predictedValue"])),"designation"]))
              
              if(modelType %in% c("pblup","ssgblup")){ # we need to calculate NRM
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
              if(modelType %in% c("gblup","ssgblup")){ # we need to calculate GRM
                commonBetweenMandP <- intersect(rownames(Markers),designationFlevels)
                if(length(commonBetweenMandP) < 2){ stop("Markers could not be matched with phenotypes. Please ensure that you have used the right marker file or check the rownames of your marker matrix and ids of your phenotypes.", call. = FALSE)                }
                M <- Markers[commonBetweenMandP,]
                if(ncol(M) > 5000){ # we remove that many markers if a big snp chip
                  A <- sommer::A.mat(M[,sample(1:ncol(M), 5000)])
                }else{ A <- sommer::A.mat(M) };  M <- NULL
                if(modelType == "ssgblup"){ # only if ssgblup we merge
                  A <- sommer::H.mat(N,A, tau=1,  omega=1, tolparinv=1e-6)
                }
              }
              if(modelType %in% c("pblup")){ A <- N  }
              
              badGeno <- which(rownames(A) == "") # should have no covariance with anyone
              if(length(badGeno) > 0){A[badGeno,2:ncol(A)]=0; A[2:nrow(A),badGeno]=0} # make zero covariance with this genotype
              badBlankGenotype <- which(colnames(A)=="")
              if(length(badBlankGenotype) > 0){A <- A[-badBlankGenotype,-badBlankGenotype]}
              inter <- intersect(designationFlevels,colnames(A)) # go for sure
              if(modelType %in% c("gblup","ssgblup")){
                onlyInA <- setdiff(rownames(Markers),designationFlevels)
              }else{
                onlyInA <- setdiff(rownames(A),designationFlevels)
              }
              differ <- setdiff(designationFlevels,inter) # are missing in A matrix
              genoMetaData <- list(withMarkandPheno=inter, withPhenoOnly=differ, withMarkOnly=onlyInA)
              # get inverse matrix
              if(length(inter) > 0){ #
                A1inv <- solve(A[inter,inter] + diag(tolParInv,length(inter), length(inter)))
              }else{A1inv <- matrix(0,0,0)}
              if(length(differ) > 0){ # we have to add individuals without markers or not being part of the GRM?
                if(length(differ) > 1){ # there's at least 2 inds to be added
                  A2inv <- diag(x=rep(mean(diag(A1inv)),length(differ)))
                }else{ A2inv <- diag(1)*mean(diag(A1inv)) } # there's only one individual to be added
                colnames(A2inv) <- rownames(A2inv) <- differ
              }else{A2inv <- matrix(0,0,0)}
              Ainv <- enhancer::adiag1(A1inv,A2inv)
              Ainv[lower.tri(Ainv)] <- t(Ainv)[lower.tri(Ainv)] # fill the lower triangular
              colnames(Ainv) <- rownames(Ainv) <- c(colnames(A1inv), colnames(A2inv))
              A1inv <- NULL; A2inv <- NULL;
              levelsInAinv <- colnames(Ainv)
              myGinverse <- list(designation=Ainv)
            }else{ # blue model
              inter <- character()
              onlyInA <- character() # genotypes only present in A and not in dataset
              differ <- character()
            }
          }
          if(returnFixedGeno){myGinverse <- NULL}
          if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
          if(useWeights){
            weightsFormulation="w"
            if(verbose){print("Using weights in the analysis. Residual variance will be fixed to 1.")  }
          }else{
            weightsFormulation=NULL
            if(verbose){print("Ignoring weights in the analysis. Residual variance will be estimated.")  }
          }
          # options(spam.cholsymmetrycheck=FALSE)
          mix <- try(
            LMMsolver::LMMsolve(fixed =as.formula(fix),
                                random = ranFormulation,
                                residual=ranres,
                                weights = weightsFormulation,
                                ginverse = myGinverse,
                                group = groupTrait,
                                family = eval(parse(text = traitFamily[iTrait])),
                                data = mydataSub, maxit = maxIters),
            silent = TRUE
          )
          myGinverse <- NULL      # mix
          
          # print(mix$VarDf)
          if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
            ## save the modeling used
            currentModeling <- data.frame(module="mta", analysisId=mtaAnalysisId,trait=iTrait, environment="across",
                                          parameter=c("fixedFormula","randomFormula","family","designationEffectType"), 
                                          value=c(fix,ifelse(returnFixedGeno,NA,ranran),traitFamily[iTrait], toupper(modelType) ))
            phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
            ## save the environments used goodFields
            currentModeling <- data.frame(module="mta", analysisId=mtaAnalysisId,trait=iTrait, environment=allEnvironments,
                                          parameter="includedInMta", 
                                          value=ifelse(allEnvironments%in%goodFields, TRUE, FALSE))
            phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
            ##
            if(is.null(phenoDTfile$metadata$weather)){numericMetas <- character()}
            for(iIndex in c(numericMetas,"envIndex")){ # iIndex="latitude"
              if( (iIndex %in% interactionsWithGenoTrait) ){names(mix$ndxCoefficients[[paste0("designation:",iIndex)]]) <- names(mix$ndxCoefficients$designation) } # copy same names than main designation effect
            }
            
            iGenoUnit <- "designation" # in MET iGenoUnit is always "designation" only in STA we allow for different
            
            ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(modelType == "rrblup"){ # user wants to predict using marker effects
              ##
              Mfull <- Markers
              if((ncol(Mfull) < nrow(Mfull)) | nPC==0){M2 <- Mfull}else{  M2 <- tcrossprod(Mfull)}
              xx2 = with(mydataSub, enhancer::redmm(x=designation, M=M2, nPC=nPC, returnLam = TRUE)) # we need the new lambda for the fullt marker matrix
              ss <- mix$VarDf;  rownames(ss) <- ss$VarComp
              pp <- list()
              for(iGroup in names(LGrp)){ # iGroup <- names(LGrp)[1]  for each rrBLUP effect
                shouldBeOne <- which(mix$ndxCoefficients[[iGroup]] == 0)
                if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGroup]][shouldBeOne] = 1}
                blup <- mix$coefMME[mix$ndxCoefficients[[iGroup]]]
                names(blup) <- names(mix$ndxCoefficients[[iGroup]])
                names(blup) <- gsub(paste0(iGroup,"_"),"",names(blup))
                predictedValue <- (xx2$Lam %*% blup[colnames(xx2$Lam)]) + mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
                if(length(shouldBeOne) > 0){predictedValue[1] = mix$coefMME[mix$ndxCoefficients$`(Intercept)`]}
                # names(predictedValue) <- names(mix$ndxCoefficients[[iGroup]])
                dims <- mix$EDdf
                start <- sum(dims[1:(which(dims$Term == iGroup) - 1),"Model"]) # we don't add a one because we need the intercept
                nEffects <- length(mix$coefMME[mix$ndxCoefficients[[iGroup]]])
                pev <- as.matrix(solve(mix$C))[start:(start+nEffects-1),start:(start+nEffects-1)]
                pev <- xx2$Lam %*% pev %*% t(xx2$Lam)
                stdError <- (sqrt(Matrix::diag(pev)))
                Vg <- (ss[iGroup,"Variance"]*ncol(Markers))/2
                reliability <- abs((Vg - Matrix::diag(pev))/Vg)
                badRels <- which(reliability > 1); if(length(badRels) > 0){reliability[badRels] <- 0.9999}
                badRels2 <- which(reliability < 0); if(length(badRels2) > 0){reliability[badRels2] <- 0}
                pp[[iGroup]] <- data.frame(designation=rownames(predictedValue), predictedValue=predictedValue, stdError=stdError, reliability=reliability,
                                           trait=paste0(iTrait,"_",iGroup) )
                cv <- (sd(predictedValue,na.rm=TRUE)/mean(predictedValue,na.rm=TRUE))*100
                Ve <- Ve - Vg
                ## save metrics
                phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                             data.frame(module="mta",analysisId=mtaAnalysisId, trait= paste0(iTrait,"_",iGroup), environment="across",
                                                        parameter=c("mean","CV", "r2","Vg"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML"),
                                                        value=c(mean(predictedValue, na.rm=TRUE), cv, median(reliability), var(predictedValue, na.rm=TRUE)  ),
                                                        stdError=c(0,0,sd(reliability, na.rm = TRUE)/sqrt(length(reliability)),0)
                                             )
                )
                counter <- counter+1
              }
              phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                           data.frame(module="mta",analysisId=mtaAnalysisId, trait= iTrait, environment="across",
                                                      parameter=c("Vr"), method=c("REML"), value=c( Ve ),stdError=c(0) )
              )
              pp <- do.call(rbind,pp)
              ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            }else{ # user wants to do relationship-based genetic evaluation
              
              
              if(length(grep("designation", ranran)) > 0){ # geno was random
                ######################################################
                ## do genetic evaluation only if there was genotype as random effect
                ######################################################
                if(modelType %in% c("gblup","ssgblup")){ # if user provided markers
                  allIndsNames <- rownames(Markers) # individuals in the marker matrix
                }else{# if user provided relationship matrix
                  allIndsNames <- colnames(A) # individuals in the relationship matrix
                }
                indNamesInAinv <- levelsInAinv # reduced set of individuals in the met analysis
                indNamesInAinvNoDiffer <- setdiff(indNamesInAinv,differ) # individuals with no markers data that were just expanded in the A matrix
                indsToBePredicted <- c(indNamesInAinv, setdiff(allIndsNames, indNamesInAinv) )# additional individuals to be predicted in the A matrix
                # Ainv <- NULL
                if(length(indsToBePredicted) > 0 ){ # if there's additional inds to be predicted
                  groups0 <- sort(rep(paste0("g",1:1000),batchSizeToPredict)) # groups of 500
                  grouping <- data.frame(id=indsToBePredicted,grp=groups0[1:length(indsToBePredicted)])
                  grouping$grp[1:length(indNamesInAinv)] <- "g1" ## all the ones in the model go to group 1
                  groups <- unique(grouping$grp)
                }else{ # no additional inds
                  grouping <- data.frame(id=indNamesInAinv,grp="g1")
                  groups <- unique(grouping$grp)
                }
                nEffectsToStore <- length(interactionsWithGenoTrait)+1
                genEva <- lapply(vector(mode="list",nEffectsToStore), function(x,nn){vector(mode="list",nn)}, nn=4) #list of lists
                names(genEva) <- names(mix$ndxCoefficients)[grep("designation",names(mix$ndxCoefficients))]
                genEva <- lapply(genEva, function(x){names(x) <- c("predictedValue","stdError","pev","R2"); return(x)})
                ttGrouping <- table(grouping$grp)
                for(iGroup in 1:length(groups)){ # iGroup=1  # for each group do the genetic evaluation
                  if(verbose){print(paste("Predicting batch",iGroup))}
                  iGrupNames <- grouping[which(grouping$grp == groups[iGroup]),"id"]
                  # make a new A for the ith group
                  present <- unique(c( intersect(colnames(A), indNamesInAinv), iGrupNames ))
                  mydataSub2 <- mydataSub[which(!is.na(mydataSub$predictedValue)),]
                  mydataSub2 <- droplevels(mydataSub2[which(mydataSub2$designation %in% c(present)),])
                  if(modelType %in% c("gblup","ssgblup")){ # if user provided markers
                    AGE <- sommer::A.mat(Markers[intersect(present,rownames(Markers)),])
                  }else{
                    present2 <- intersect(colnames(A),present)
                    AGE <- A[present2,present2]
                  }
                  ######################################
                  ## do the genetic evaluation
                  genEvaProv <- cgiarBase::geneticEvaluation3(
                    fixed =as.formula(fix),
                    random = as.formula(ranran),
                    rcov=ranres,
                    weights = "w",
                    data = mydataSub2,
                    varComp=mix$VarDf,
                    keep=ttGrouping[iGroup],
                    AFGE=AGE
                  )
                  AGE <- NULL
                  # print("finished")
                  # fill the list
                  for(uu in 1:length(genEvaProv)){ # for each genetic effect
                    genEva[[uu]]$predictedValue <- c(genEva[[uu]]$predictedValue,genEvaProv[[uu]]$predictedValue )
                    genEva[[uu]]$stdError <- c(genEva[[uu]]$stdError,genEvaProv[[uu]]$stdError )
                    if(iGroup == 1){
                      genEva[[uu]]$pev <- genEvaProv[[uu]]$pev
                      genEva[[uu]]$R2 <- genEvaProv[[uu]]$R2
                    }else{
                      genEva[[uu]]$pev <- Matrix::bdiag(genEva[[uu]]$pev,genEvaProv[[uu]]$pev )
                      genEva[[uu]]$R2 <- c(genEva[[uu]]$R2,genEvaProv[[uu]]$R2 ) #Matrix::bdiag(genEva[[uu]]$R2,genEvaProv[[uu]]$R2 )
                    }
                  }
                  genEvaProv <- NULL
                }
                ######################################################
                ## end of genetic evaluation function
                ######################################################
              }
              ## now extract the needed values
              if(iGenoUnit %in% fixedTermTrait){ # user wants fixed effect predictions for genotype
                shouldBeOne <- which(mix$ndxCoefficients[[iGenoUnit]] == 0)
                if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGenoUnit]][shouldBeOne] = 1}
                predictedValue <- mix$coefMME[mix$ndxCoefficients[[iGenoUnit]]] + mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
                if(length(shouldBeOne) > 0){predictedValue[1] = mix$coefMME[mix$ndxCoefficients$`(Intercept)`]}
                names(predictedValue) <- names(mix$ndxCoefficients[[iGenoUnit]])
                dims <- mix$EDdf
                start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                pev <- as.matrix(solve(mix$C))[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                stdError <- (sqrt(diag(pev)))
                designation <- gsub("designation_","", names(predictedValue))
              }else{ # user wants random effect predictions for genotype (main effect)
                predictedValue <- genEva$designation$predictedValue
                stdError <- genEva$designation$stdError
                pev <- genEva$designation$pev
                designation <- gsub("designation","", names(predictedValue))
              }
              pp <- data.frame(designation,predictedValue,stdError)
              ss = mix$VarDf; rownames(ss) <- ss$VarComp
              if(iGenoUnit %in% fixedTermTrait){Vg=0}else{Vg <- ss["designation",2]; }
              Ve <- Ve - Vg
              if(iGenoUnit %in% fixedTermTrait){ # add reliabilities to the data frame
                pp$reliability <- NA
              }else{ # if random, reliability can be calculated for main effect
                # print(names(genEva))
                pp$reliability <- genEva$designation$R2
                badRels <- which(pp$reliability > 1); if(length(badRels) > 0){pp$reliability[badRels] <- 0.9999}
                badRels2 <- which(pp$reliability < 0); if(length(badRels2) > 0){pp$reliability[badRels2] <- 0}
              }
              pp$trait <- iTrait # add trait
              cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
              ## save metrics
              phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                           data.frame(module="mta",analysisId=mtaAnalysisId, trait=iTrait, environment="across",
                                                      parameter=c("mean","CV", "r2","Vg"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML"),
                                                      value=c(mean(pp$predictedValue, na.rm=TRUE), cv, mean(pp$reliability), Vg),
                                                      stdError=c(0,0,sd(pp$reliability, na.rm = TRUE)/sqrt(length(pp$reliability)),0)
                                           )
              )
              ## genetic variances
              lpv <- sum(mix$EDdf$Model[1:which(mix$EDdf$Term == "designation")])+1 # to be used as a starting point if random regression is requested
              # extract sensitivities if interaction is requested
              if(length(interactionsWithGenoTrait) > 0){ # if there's interactions
                if( length(intersect(interactionsWithGenoTrait, c("envIndex","timePoint","latitude","longitude","altitude","weather1","weather2"))) > 0 ){
                  
                  for(iInteractionTrait in c("envIndex","timePoint","latitude","longitude","altitude","weather1","weather2")){ # iInteractionTrait="envIndex"
                    
                    if( (iInteractionTrait %in% interactionsWithGenoTrait) ){ # iInteractionTrait="envIndex"
                      
                      counter <- counter+1
                      iGenoUnit <- paste0("designation:",iInteractionTrait)#"designation:envIndex"
                      if(iGenoUnit %in% fixedTermTrait){ # user wants fixed effect predictions for genotype:envIndex
                        shouldBeOne <- which(mix$ndxCoefficients[[iGenoUnit]] == 0)
                        if(length(shouldBeOne) > 0){mix$ndxCoefficients[[iGenoUnit]][shouldBeOne] = 1}
                        predictedValue <- mix$coefMME[mix$ndxCoefficients[[iGenoUnit]]] #+ mix$coefficients$`(Intercept)`
                        dims <- mix$EDdf
                        start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                        stdError <- (sqrt(diag(as.matrix(solve(mix$C)))))[start:(start+length(predictedValue)-1)]
                        designation <- gsub("designation_","", names(predictedValue))
                      }else{ # user wants random effect predictions for genotype:envIndex
                        predictedValue <- genEva[[iGenoUnit]]$predictedValue #+ mix$coefficients$`(Intercept)`
                        stdError <- genEva[[iGenoUnit]]$stdError
                        designation <- gsub("designation","", names(predictedValue))
                      }
                      pp2 <- data.frame(designation,predictedValue,stdError)
                      Vg <- ss[iGenoUnit,2];
                      Ve <- Ve - Vg
                      if(iGenoUnit %in% fixedTermTrait){pp$reliability <- 1e-6}else{pp2$reliability <- genEva[[iGenoUnit]]$R2}
                      pp2$trait <- paste(iTrait,iInteractionTrait,sep="-")
                      cv <- (sd(pp2$predictedValue,na.rm=TRUE)/mean(pp2$predictedValue,na.rm=TRUE))*100
                      ## save metrics
                      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                                   data.frame(module="mta",analysisId=mtaAnalysisId, trait=paste(iTrait,iInteractionTrait,sep="-"),
                                                              environment="across",
                                                              parameter=c("mean","CV", "r2","Vg"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML"),
                                                              value=c(mean(pp2$predictedValue, na.rm=TRUE), cv, mean(pp2$reliability), Vg),
                                                              stdError=c(0,0,sd(pp2$reliability, na.rm = TRUE)/sqrt(length(pp2$reliability)),0)
                                                   )
                      )
                      pp <- rbind(pp,pp2) # bind predictions
                      
                    }
                    
                  }
                  
                }
              } # if there's even interactions
              # add the error variance
              phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                           data.frame(module="mta",analysisId=mtaAnalysisId, trait= iTrait, environment="across",
                                                      parameter=c("Vr"), method=c("REML"), value=c( Ve ),stdError=c(0) )
              )
              
            }
            ###%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
          }else{ # if model failed
            if(modelType == "rrblup"){
              cat(paste("Model failed for trait:", iTrait))
            }else{
              if(verbose){ cat(paste("Mixed model failed for this combination. Aggregating and assuming h2 = 0 \n"))}
              pp <- aggregate(predictedValue ~ designation, FUN=mean, data=mydataSub)
              pp$stdError <- sd(pp$predictedValue)  # 1
              # pp$status <- "Aggregated"
              pp$reliability <- 1e-6
              pp$trait <- iTrait
              cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
              ## save metrics
              phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                           data.frame(module="mta",analysisId=mtaAnalysisId, trait=iTrait,
                                                      environment="across",
                                                      parameter=c("mean","CV", "r2","Vg","Vr"), method=c("sum(x)/n","sd/mu","(G-PEV)/G","REML","REML"),
                                                      value=c(mean(pp$predictedValue, na.rm=TRUE), cv, 0, 0, Ve),
                                                      stdError=c(0,0,0,0, 0)
                                           )
              )
              currentModeling <- data.frame(module="mta", analysisId=mtaAnalysisId,trait=iTrait, environment="across",
                                            parameter=c("fixedFormula","randomFormula","family","designationEffectType"), 
                                            value=c("None","None","None","mean"))
              phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
              
            }
          }
          #######################################
          #######################################
          #######################################
          # Trait run is finished. If model worked well save the results
          if(!inherits(mix,"try-error") ){ 
            pp$environment <- "across"
            mydataForEntryType <- droplevels(mydata[which(mydata$trait == iTrait),])
            pp$entryType <- apply(data.frame(pp$designation),1,function(x){
              found <- which(mydataForEntryType$designation %in% x)
              if(length(found) > 0){
                x2 <- paste(sort(unique(toupper(trimws(mydataForEntryType[found,"entryType"])))), collapse = "#");
              }else{x2 <- "unknown"}
              return(x2)
            })
            mydataForEntryType <- NULL
            if(modelType == "rrblup"){
              pp$entryType <- ifelse(as.character(pp$designation) %in% rownames(Mtrait) ,"GEBV_tested", "GEBV_predicted")
            }else{
              pp$entryType <- paste(ifelse(as.character(pp$designation) %in% differ, "TGV", surrogate),
                                    pp$entryType,
                                    ifelse(as.character(pp$designation) %in% onlyInA, "predicted", "tested"),
                                    sep="_")
            }
            ###
            predictionsList[[counter2]] <- pp;
            counter=counter+1
          }
          
        }
      }
    }
    counter2 = counter2+1
  }
  #
  trait <- setdiff(trait,traitToRemove)
  heritLB <- heritLB[which(traitOrig %in% trait)]
  heritUB <- heritUB[which(traitOrig %in% trait)]
  if(length(predictionsList) == 0){stop("There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all envs.",call. = FALSE)}
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- mtaAnalysisId
  ##########################################
  ## add timePoint of origin, stage and designation code
  entries <- unique(mydata[,"designation"])
  baseOrigin <- do.call(rbind, apply(data.frame(entries),1,function(x){
    out1 <- (sort(mydata[which(mydata$designation %in% x),"gid"], decreasing = FALSE))[1]
    out2 <- (sort(mydata[which(mydata$designation %in% x),"mother"], decreasing = FALSE))[1]
    out3 <- (sort(mydata[which(mydata$designation %in% x),"father"], decreasing = FALSE))[1]
    out4 <- paste(unique(sort(mydata[which(mydata$designation %in% x),"pipeline"], decreasing = FALSE)),collapse=", ")
    y <- data.frame(designation=x,gid=out1,mother=out2,father=out3,pipeline=out4)
    return(y)
  }))
  predictionsBind <- merge(predictionsBind,baseOrigin, by="designation", all.x=TRUE)
  predictionsBind$module <- "mta"
  #########################################
  ## update databases
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,
                                   predictionsBind[,colnames(phenoDTfile$predictions)])
  phenoDTfile$status <- rbind( phenoDTfile$status, data.frame(module="mta", analysisId=mtaAnalysisId))
  ## add which data was used as input
  modeling <- data.frame(module="mta",  analysisId=mtaAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId ))
  # if(!is.null(analysisIdForGenoModifications)){
  #   modeling <- rbind(modeling, data.frame(module="mta",  analysisId=mtaAnalysisId, trait=c("inputObject"), environment="general",
  #                                          parameter= c("analysisId"), value= c(analysisIdForGenoModifications ))
  #   )
  # }
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)
}
