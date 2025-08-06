metLMMsolver <- function(
    phenoDTfile= NULL, analysisId=NULL, analysisIdGeno = NULL,
    fixedTerm= list("1"),  randomTerm=NULL, expCovariates=NULL,
    envsToInclude=NULL, trait= NULL, traitFamily=NULL, useWeights=TRUE,
    calculateSE=TRUE, heritLB= 0.15,  heritUB= 0.95,
    meanLB=0, meanUB=Inf, nPC=NULL,   # subsetVariable=NULL, subsetVariableLevels=NULL,
    maxIters=50,  verbose=TRUE
){

  #print(nPC)
  ## THIS FUNCTION PERFORMS A MULT TRIAL ANALYSIS USING LMM SOLVER
  mtaAnalysisId <- as.numeric(Sys.time())
  namesSeq <- function(x){
    nCharX <- nchar(x)
    maxZeros <- max(nCharX)
    nZeros <- abs(nCharX - maxZeros)
    zeros <- apply(data.frame(nZeros),1,function(x){paste(rep("0",x), collapse = "")})
    res <- paste0(zeros, as.character(x))
    return(res)
  }
  sommerVersion <- as.numeric(paste(strsplit(as.character(packageVersion("sommer")),"[.]")[[1]][1:2], collapse = ""))
  ##########################################
  ##########################################
  ## CONTROLS FOR MISSPECIFICATION (6 lines)
  if(is.null(phenoDTfile)){stop("Please provide the phenotype file", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the STA analysisId to be analyzed", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}else{
    baseData <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId ),]
    if(length(intersect(trait, unique(baseData[,"trait"]))) == 0){stop("The traits you have specified are not present in the analysisId provided.", call. = FALSE)}
  }
  if(is.null(traitFamily)){traitFamily <- rep("quasi(link = 'identity', variance = 'constant')", length(trait))}
  if(!is.null(randomTerm)){
    if(length(randomTerm) == 0){randomTerm <- expCovariates <- NULL}else{
      flatCovars <- unlist(expCovariates)
      if (!"genoD" %in% flatCovars) {
        randomTerm <- unique(randomTerm)
      }
      if(is.null(expCovariates)){expCovariates <- randomTerm; expCovariates <- lapply(expCovariates, function(x){rep("none",length(x))})}else{
        if(length(expCovariates) != length(randomTerm)){
          stop("Please ensure that expCovariates and randomTerm arguments have the same length.", call. = FALSE)
        }else{
          if( sum(unlist(mapply('-', lapply(randomTerm,length), lapply(expCovariates,length), SIMPLIFY = FALSE))) != 0){
            stop("Please ensure that expCovariates and randomTerm arguments have the same length.", call. = FALSE)
          }
        }
      }
    }
  }
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  if(length(fixedTerm) == 0 | is.null(fixedTerm)){fixedTerm <- "1"}else{fixedTerm <- unique(fixedTerm)}
  traitsForExpCovariates <- unique(phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId),"trait"])
  if(!is.null(nPC)){if(length(intersect(names(nPC), c("geno","weather","pedigree", traitsForExpCovariates))) == 0){stop("The nPC argument needs to be a named numeric vector with names 'geno', 'weather' and 'pedigree' or traits available.", call. = FALSE)}}

  ##########################################
  ##########################################
  ## EXTRACT POSSIBLE EXPLANATORY COVARIATES AND FORM KERNELS (30 lines)
  Weather <- cgiarPipeline::summaryWeather(phenoDTfile, wide=TRUE) # in form of covariates
  if(nrow(Weather) > 1){
    Weather <- apply(Weather,2,sommer::imputev)
    colnames(Weather) <- gsub(" ","",colnames(Weather))
  }
  covars <- unique(unlist(expCovariates))
  randomTermForCovars <- unique(unlist(randomTerm))
  if(!is.null(randomTermForCovars)){
    if( any( covars %in% c("genoA","genoAD", "weather","pedigree", traitsForExpCovariates ) ) ){
      if(verbose){message("Checking and calculating kernels requested")}
      ## MARKER KERNEL
      Markers <- as.matrix(phenoDTfile$data$geno) # in form of covariates
      if(any(c("genoA","genoAD") %in% covars) & !is.null(Markers)){
        classify <- unique(unlist(randomTerm)[which(unlist(expCovariates) %in% c("genoA","genoAD","genoD") )])
        # eventually we may have to do a for loop
        if(verbose){message(paste("   Marker kernel for",paste(classify,collapse = " and "), "requested"))}
        if(analysisIdGeno == '' | is.null(analysisIdGeno)){ # user didn't provide a modifications id
          if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and you have not provided a modifications table to impute the genotype data. Please go to the 'Markers QA/QC' module prior to run a model with genoA or genoAD covariate.", call. = FALSE)}
        }else{ # user provided a modifications Id
          if(class(phenoDTfile$data$geno)[1] == "genlight"){
            theresMatch <- which(as.character(analysisIdGeno) %in% names(phenoDTfile$data$geno_imp))
          } else{
            modificationsMarkers <- phenoDTfile$modifications$geno
            theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdGeno)
          }

          if(length(theresMatch) > 0){ # there's a modification file after matching the Id
            if(class(phenoDTfile$data$geno)[1] == "genlight"){
              Markers <- as.matrix(phenoDTfile$data$geno_imp[[as.character(analysisIdGeno)]])
            } else{
              modificationsMarkers <- modificationsMarkers[theresMatch,]
              Markers <- cgiarBase::applyGenoModifications(M=Markers, modifications=modificationsMarkers)
            }
          }else{ # there's no match of the modification file
            if(length(which(is.na(Markers))) > 0){stop("Markers have missing data and your Id didn't have a match in the modifications table to impute the genotype data.", call. = FALSE)}
          }
        }
        # qas <- which( phenoDTfile$status$module == "qaGeno" ); qas <- qas[length(qas)]
        # if(length(qas) > 0){
        #   modificationsMarkers <- phenoDTfile$modifications$geno[which(phenoDTfile$modifications$geno$analysisId %in% qas ),]
        #   Markers <- cgiarBase::applyGenoModifications(M=Markers, modifications=modificationsMarkers)
        #   # if(length(which(is.na(Markers))) > 0){Markers <- apply(Markers,2,sommer::imputev)}
        # }else{
        #   missing <- apply(Markers,2,sommer::propMissing)
        #   Markers <- apply(Markers[,which(missing < 0.9)],2,sommer::imputev)
        # }
        if(nPC["geno"] < 0){ # do not include extra individuals
          mydataX <-  phenoDTfile$predictions[which( phenoDTfile$predictions$analysisId %in% analysisId),]
          Markers <- Markers[which(rownames(Markers) %in% unique(mydataX$designation) ), ]
          if(verbose){message(paste("Subsetting marker to",nrow(Markers),"individuals present"))}
        }
        ploidyFactor <- max(Markers)/2
        if("genoA" %in% covars){G <- sommer::A.mat(Markers-ploidyFactor);} # additive model
        if("genoD" %in% covars){ #Dominance kernel
          if(ploidyFactor == 1){

            D <- 1-abs(Markers)

            f <- rowSums(D) / ncol(D) #inbreeding fixed eff
            names(f) <- rownames(Markers)

            D <- sommer::D.mat(D)
            D <- D + diag(1e-5, ncol(D), ncol(D))
          }else{ #autopolyploid formula for digenic dominance (Batista et al. 2022)

            ploidy <- ploidyFactor * 2

            dom_matrix = Markers/ploidy
            dom_matrix = 4*dom_matrix - 4*(dom_matrix*dom_matrix)

            f <- rowSums(dom_matrix) / ncol(dom_matrix) #inbreeding fixed eff
            names(f) <- rownames(Markers)

            MAF <- colMeans(Markers, na.rm = TRUE) / ploidy
            tMarkers <- t(Markers)

            C_mat <- matrix(choose(ploidy, 2), nrow = nrow(tMarkers), ncol = ncol(tMarkers))
            Ploidy_mat <- matrix(ploidy, nrow = nrow(tMarkers), ncol = ncol(tMarkers))

            Q <- (MAF^2 * C_mat) -
              (Ploidy_mat - 1) * MAF * tMarkers +
              0.5 * tMarkers * (tMarkers-1)

            D <- crossprod(Q)
            denomDom <- sum(C_mat[,1]*MAF^2*(1-MAF)^2)
            D <- D/denomDom
            D <- D + diag(1e-5, ncol(D), ncol(D))
          }
          Dchol <- t(chol(D))
        }
        if("genoAD" %in% covars){ # if genetic value is desired let's do a log marker model
          Markers <- apply(Markers+1,2,log)
          G <- sommer::A.mat(Markers)
        } # additive + dominance model
        G <- G + diag(1e-5, ncol(G), ncol(G))
        Gchol <- t(chol(G))
        if(nPC["geno"] > 0){
            if(verbose){message("   Eigen decomposition of marker kernel requested")}
            decomp <- RSpectra::svds(Gchol, k = min(c(nPC["geno"], ncol(Gchol))), which = "LM")
            rownames(decomp$u) <- rownames(G); colnames(decomp$u) <- paste0("PC",namesSeq(1:ncol(decomp$u)))
            Gchol <- decomp$u
        }
      }
      ## WEATHER KERNEL
      if("weather" %in% covars & !is.null(Weather)){
        classify <- unique(unlist(randomTerm)[which(unlist(expCovariates) %in% "weather")])
        if(verbose){message(paste("   Weather kernel for",paste(classify,collapse = " and "), "requested"))}
        WeatherK <- Weather
        rownamesWeather <- rownames(WeatherK)
        WeatherK <- apply(WeatherK, 2, scale)
        WeatherK <- WeatherK[,which( !is.na(apply(WeatherK,2,var)) ), drop=FALSE]
        rownames(WeatherK) <- rownamesWeather
        if(nPC["weather"] < 0){ # do not include extra individuals
          mydataX <-  phenoDTfile$predictions[which( phenoDTfile$predictions$analysisId %in% analysisId),]
          WeatherK <- WeatherK[which(rownames(WeatherK) %in% unique(mydataX$environment) ), ]
          if(verbose){message(paste("Subsetting weather matrix to",nrow(WeatherK),"environments present"))}
        }
        W <- sommer::A.mat(WeatherK)
        W <- W + diag(1e-5, ncol(W), ncol(W))
        Wchol <- t(chol(W))
        if(nPC["weather"] > 0){
          if(verbose){message("   Eigen decomposition of weather kernel requested")}
          decomp <- RSpectra::svds(Wchol, k = min(c(nPC["weather"], ncol(Wchol))), which = "LM")
          rownames(decomp$u) <- rownames(W); colnames(decomp$u) <- paste0("PC",namesSeq(1:ncol(decomp$u)))
          Wchol <- decomp$u
        }
      }
      ## PEDIGREE KERNEL
      Pedigree <- phenoDTfile$data$pedigree
      if("pedigree" %in% covars  &  !is.null(Pedigree)){
        classify <- unique(unlist(randomTerm)[which(unlist(expCovariates) %in% "pedigree")])
        if(verbose){message(paste("   Pedigree kernel for",paste(classify,collapse = " and "), "requested"))}
        paramsPed <- phenoDTfile$metadata$pedigree
        N <- cgiarBase::nrm2(pedData=phenoDTfile$data$pedigree,
                             indivCol = paramsPed[paramsPed$parameter=="designation","value"],
                             damCol = paramsPed[paramsPed$parameter=="mother","value"],
                             sireCol = paramsPed[paramsPed$parameter=="father","value"]
        )
        if(nPC["pedigree"] < 0){ # do not include extra individuals
          mydataX <-  phenoDTfile$predictions[which( phenoDTfile$predictions$analysisId %in% analysisId),]
          N <- N[which(rownames(N) %in% unique(mydataX$designation) ), which(rownames(N) %in% unique(mydataX$designation) ) ]
          if(verbose){message(paste("Subsetting pedigree to",nrow(N),"individuals present"))}
        }
        Nchol <- t(chol(N))
        if(nPC["pedigree"] > 0){
          if(verbose){message("   Eigen decomposition of pedigree kernel requested")}
          decomp <- RSpectra::svds(Nchol, k = min(c(nPC["pedigree"], ncol(Nchol))), which = "LM")
          rownames(decomp$u) <- rownames(N); colnames(decomp$u) <- paste0("PC",namesSeq(1:ncol(decomp$u)))
          Nchol <- decomp$u
        }
      } # now is in the form of covariates
      # TRAIT-BASED KERNEL (ALWAYS ROW-GROUPED BY DESIGNATION)
      if( any(covars %in% traitsForExpCovariates) ){
        TraitKernels <- list()
        covarsTraits <- intersect(covars,traitsForExpCovariates)
        # covarsTraits <- unlist(expCovariates)[which(unlist(expCovariates) %in% comset)]
        Schol <- list()
        for(iCovar in covarsTraits){ # iCovar = covarsTraits[1] # for each trait specified in covar
          classify <- unlist(randomTerm)[which(unlist(expCovariates) %in% iCovar)] # identify at what levels should the trait be classified
          for(iClassify in classify){
            if(verbose){message(paste("  ",iCovar,"kernel for",iClassify, "requested"))}
            if(iClassify == "designation"){ # not allowed
              Schol <- diag(1); rownames(Schol) <- colnames(Schol) <- "A"
            }else{
              baseData <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId & (phenoDTfile$predictions$trait == iCovar) ),]
              ww <- as.data.frame(Weather); ww$environment <- rownames(ww)
              baseData <- merge(baseData, ww, by="environment", all.x = TRUE)
              if( unlist(lapply(baseData,class))[iClassify] %in% c("numeric","integer") ){
                Schol <- diag(1); rownames(Schol) <- colnames(Schol) <- "A"
              }else{ # iClassify is a character or a factor
                wideTrait <- reshape(baseData[,c(iClassify,"designation","predictedValue")], direction = "wide",
                                     idvar = "designation", timevar = iClassify, v.names = "predictedValue", sep= "_")
                wideTrait <- apply(wideTrait[,-1],2,sommer::imputev)
                colnames(wideTrait) <- gsub("predictedValue_","",colnames(wideTrait))
                S <- cov(wideTrait)
                S <- as.matrix(Matrix::nearPD(x = S, corr = FALSE,
                                              keepDiag = FALSE, base.matrix = FALSE, do2eigen = TRUE,
                                              doSym = FALSE, doDykstra = TRUE, only.values = FALSE,
                                              ensureSymmetry = !isSymmetric(S), eig.tol = 1e-06,
                                              conv.tol = 1e-07, posd.tol = 1e-08, maxit = 100, conv.norm.type = "I",
                                              trace = FALSE)$mat)
                Schol <- t(chol(S))
                if(nPC[iCovar] > 0){
                  if(verbose){message(paste("   Eigen decomposition of",iCovar," classified by",classify, "kernel requested"))}
                  decomp <- RSpectra::svds(Schol, k = min(c(nPC[iCovar], ncol(Schol))), which = "LM")
                  rownames(decomp$u) <- rownames(S); colnames(decomp$u) <- paste0("PC",namesSeq(1:ncol(decomp$u)))
                  Schol <- decomp$u
                }
              }
            }
            TraitKernels[[iCovar]][[iClassify]] <- Schol
          } # end of for each classify
        } # end of for each iCovar or trait
      }# end of if statement for a trait-based kernel
    }# end of if statement for any kernel
  }
  # print(dim(Nchol))
  ##########################################
  ##########################################
  ## COMPLETE THE CLEANING PARAMETERS (7 lines)
  names(traitFamily) <- trait
  heritLB <- rep(heritLB,length(trait)); heritLB <- heritLB[1:length(trait)]; names(heritLB) <- trait
  heritUB <- rep(heritUB,length(trait)); heritUB <- heritUB[1:length(trait)]; names(heritUB) <- trait
  meanLB <- rep(meanLB,length(trait)); meanLB <- meanLB[1:length(trait)]; names(meanLB) <- trait
  meanUB <- rep(meanUB,length(trait)); meanUB <- meanUB[1:length(trait)]; names(meanUB) <- trait
  traitOrig <- trait ## ?????????
  if(length(fixedTerm) == 0 | is.null(fixedTerm)){fixedTerm <- list("1")} # assign the intercept if there's no fixed effects
  ##########################################
  ##########################################
  # LOAD THE DATASET AND EXTEND IT TO INCLUDE METADATA (16 lines)
  if(verbose){message("Loading the dataset and adding metadata.")}
  mydata <- phenoDTfile$predictions #
  mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  if (nrow(mydata) < 2) stop("Not enough data is available to perform a multi trial analysis. Please perform an STA before trying to do an MET.", call. = FALSE)
  metaPheno <- phenoDTfile$metadata$pheno[which(phenoDTfile$metadata$pheno$parameter %in% c("pipeline","stage","environment","year","season","timepoint","country","location","trial","study","management")),]
  otherMetaCols <- unique(phenoDTfile$data$pheno[,metaPheno$value,drop=FALSE])
  colnames(otherMetaCols) <- cgiarBase::replaceValues(Source = colnames(otherMetaCols), Search = metaPheno$value, Replace = metaPheno$parameter )
  otherMetaCols <- otherMetaCols[which(!duplicated(otherMetaCols[,"environment"])),,drop=FALSE] # we do this in case the users didn't define the environment properly
  mydata <- merge(mydata, otherMetaCols, by="environment", all.x = TRUE)
  WeatherRow <- as.data.frame(Weather); WeatherRow$environment <- rownames(WeatherRow)
  mydata <- merge(mydata, WeatherRow, by="environment", all.x = TRUE)
  # if(!is.null(subsetVariable) & !is.null(subsetVariableLevels)){
  #   if(subsetVariable %in% colnames(otherMetaCols)){
  #     forSubset <- which(mydata[,subsetVariable] %in% subsetVariableLevels)
  #     if(length(forSubset) > 0){mydata <- mydata[forSubset,]}
  #   }
  # }
  ##########################################
  ##########################################
  # CHECK THE ENVS TO INCLUDE PER TRAIT (6 lines)
  if(is.null(envsToInclude)){
    envsToInclude=  as.data.frame( do.call( rbind, list (with(mydata, table(environment,trait)) ) ) )
    bad <- which(envsToInclude <= 1, arr.ind = TRUE)
    if(nrow(bad) > 0){envsToInclude[bad] = 0}
    envsToInclude[which(envsToInclude > 1, arr.ind = TRUE)] = 1
  }; allEnvironments <- rownames(envsToInclude)
  ##########################################
  ##########################################
  # BUILD THE DATASETS FOR MODEL FITTING (100 lines)
  if(verbose){message("Building trait datasets.")}
  metrics <- phenoDTfile$metrics
  metrics <- metrics[which(metrics$analysisId %in% analysisId),]
  myDataTraits <- fixedTermTrait <- randomTermTrait <- groupingTermTrait <- Mtrait <- envsTrait <- entryTypesTrait <- envCount <- list()
  for(iTrait in trait){ # iTrait = trait[1]
    # filter for records available
    vt <- which(mydata[,"trait"] == iTrait)
    if(length(vt) > 0){ # we have data for the trait
      prov <- mydata[vt,]
      # filter by the environments to include
      vte <- which(prov[,"environment"] %in% rownames(envsToInclude)[as.logical(envsToInclude[,iTrait])])
      prov <- prov[vte,]
      # remove bad environment based on h2 and r2
      pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("plotH2","H2","meanR2","r2", apply(expand.grid( c("plotH2","H2","meanR2","r2"), c("designation","mother","father")),1,function(f){paste(f,collapse = "_")}) )),]
      goodFields <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value >= heritLB[iTrait]) & (pipeline_metricsSub$value <= heritUB[iTrait])),"environment"])
      prov <- prov[which(prov$environment %in% goodFields),]
      # remove bad environment based on environment means
      pipeline_metricsSub <- metrics[which(metrics$trait == iTrait & metrics$parameter %in% c("plotH2","H2","meanR2","r2", apply(expand.grid( c("mean"), c("designation","mother","father")),1,function(f){paste(f,collapse = "_")}) ) ),]
      goodFieldsMean <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > meanLB[iTrait]) & (pipeline_metricsSub$value < meanUB[iTrait])),"environment"])
      prov <- prov[which(prov$environment %in% goodFieldsMean),]
      envCount[[iTrait]] <- unique(prov$environment)

      #Add inbreeding coefficient to prov
      if ("genoD" %in% covars & (!"inbreeding" %in% colnames(prov))) {
        prov$inbreeding <- f[match(prov$designation, names(f))]
      }

      fixedTermProv <- fixedTerm
          for(iFixed in 1:length(fixedTermProv)){ # for each element in the list # iFixed=1
            fixedTermProv2 <- fixedTermProv[[iFixed]]
            for(iFixed2 in  fixedTermProv2){ # for each factor in the interactions # iFixed2 = fixedTermProv2[1]
              if(iFixed2 != "1"){
                if( length( table(prov[,iFixed2]) ) == 1 ){ fixedTermProv[[iFixed]] <- setdiff( fixedTermProv[[iFixed]], iFixed2 )}
              }
            }
          }
          fixedTermTrait[[iTrait]] <- unique(fixedTermProv[which(unlist(lapply(fixedTermProv,length)) > 0)])

      if(nrow(prov) > 0){ # if after filters there's still data for this trait we can continue and save the data
        if( var(prov[,"predictedValue"], na.rm = TRUE) > 0 ){ # check that there is variance
          # make new formula for this specific trait if it passed all the filters
          fixedTermProv <- fixedTerm
          for(iFixed in 1:length(fixedTermProv)){ # for each element in the list # iFixed=1
            fixedTermProv2 <- fixedTermProv[[iFixed]]
            for(iFixed2 in  fixedTermProv2){ # for each factor in the interactions # iFixed2 = fixedTermProv2[1]
              if(iFixed2 != "1"){
                if( length( table(prov[,iFixed2]) ) == 1 ){ fixedTermProv[[iFixed]] <- setdiff( fixedTermProv[[iFixed]], iFixed2 )}
              }
            }
          }
          fixedTermTrait[[iTrait]] <- unique(fixedTermProv[which(unlist(lapply(fixedTermProv,length)) > 0)])
          # random formula per trait
          randomTermProv <- randomTerm
          if(!is.null(randomTermProv)){
            for(irandom in 1:length(randomTermProv)){ # for each element in the list # irandom=2
              randomTermProv2 <- randomTermProv[[irandom]]
              for(irandom2 in  1:length(randomTermProv2)){ # for each factor in the interactions # irandom2 = 2
                if( length( table(prov[,randomTermProv2[irandom2]]) ) <= 1 ){ randomTermProv[[irandom]] <- setdiff( randomTermProv[[irandom]], randomTermProv2[irandom2] )}
              }
            }
          }
          # any term that is modified from what user specified we remove it totally, is better than fitting something undesired
          goodTerms <- which( ( unlist(lapply(randomTerm,length)) - unlist(lapply(randomTermProv,length)) ) == 0 )
          randomTermProv <- randomTerm[goodTerms]
          expCovariatesProv <- expCovariates[goodTerms]
          if (!"genoD" %in% unlist(expCovariatesProv)) {
            randomTermProv <- unique(randomTermProv)
          }

          # if reduced models reduce the datasets to the needed explanatory covariates
          if(!is.null(randomTermProv)){
            # reduce datasets
            for(irandom in 1:length(randomTermProv)){ # for each element in the list # irandom=1
              randomTermProv2 <- randomTermProv[[irandom]]
              for(irandom2 in  1:length(randomTermProv2)){ # for each factor in the interactions # irandom2 = 1
                # print(expCovariatesProv[[irandom]][irandom2])
                if( expCovariatesProv[[irandom]][irandom2] == "weather"){
                  M <- Wchol
                }else if(expCovariatesProv[[irandom]][irandom2] %in% c("geno","genoA","genoAD") ){
                  M <- Gchol
                }else if(expCovariatesProv[[irandom]][irandom2] == "genoD"){
                  M <- Dchol
                }else if(expCovariatesProv[[irandom]][irandom2] == "pedigree"){
                  M <- Nchol
                }else if(expCovariatesProv[[irandom]][irandom2] %in% traitsForExpCovariates){ # Trait kernel
                  classify <- randomTermForCovars[which(covars == expCovariatesProv[[irandom]][irandom2])]
                  M <- TraitKernels[[expCovariatesProv[[irandom]][irandom2]]][[classify]] # Schol equivalent
                }else{ # No kernel
                  namesZ <- unique(prov[,randomTermProv2[irandom2]])
                  M <- Matrix::Diagonal(n=length(namesZ)); rownames(M) <- colnames(M) <- namesZ
                }
                goodLevels <- intersect(unique(prov[,randomTermProv2[irandom2]]), rownames(M) )
                if(length(goodLevels) > 0){ # only if we make a match we reduce the dataset
                  prov <- prov[which(prov[,randomTermProv2[irandom2]] %in% goodLevels),]
                }else{expCovariatesProv[[irandom]][irandom2]="none"}# else we don't and change to "none" the kernel for that effect
              }
            }
          }
          ## build and add the incidence matrices
          groupingTermProv <- Mprov <- envsProv <- entryTypeProv <- list()
          if(!is.null(randomTermProv)){
            for(irandom in 1:length(randomTermProv)){ # for each element in the list # irandom=1
              randomTermProv2 <- randomTermProv[[irandom]]
              expCovariatesProv2 <- expCovariatesProv[[irandom]]
              # nExp <- numeric() # to save the number of effects
              xxList <- Mlist <- list()
              for(irandom2 in  1:length(randomTermProv2)){ # for each factor in the interactions # irandom2 = 2
                # get kernel
                if( expCovariatesProv[[irandom]][irandom2] == "weather"){
                  M <- Wchol # Weather
                }else if(expCovariatesProv[[irandom]][irandom2] %in% c("geno","genoA","genoAD") ){
                  M = Gchol # Markers
                }else if(expCovariatesProv[[irandom]][irandom2] == "genoD"){
                  M <- Dchol
                }else if(expCovariatesProv[[irandom]][irandom2] == "pedigree"){
                  M <- Nchol # Pedigree
                }else if(expCovariatesProv[[irandom]][irandom2] %in% traitsForExpCovariates){ # Trait kernel
                  classify <- randomTermForCovars[which(covars == expCovariatesProv[[irandom]][irandom2])]
                  M <- TraitKernels[[expCovariatesProv[[irandom]][irandom2]]][[classify]] # Schol equivalent
                }else{ # No kernel
                  if( unlist(lapply(prov, class))[randomTermProv2[irandom2]] %in% c("factor","character") ){
                    namesZ <- unique(prov[,randomTermProv2[irandom2]])
                    M <- Matrix::Diagonal(n=length(namesZ)); rownames(M) <- colnames(M) <- namesZ
                  } else{ # numeric or integer
                    M <- Matrix::Diagonal(n=1); rownames(M) <- colnames(M) <- randomTermProv2[irandom2]
                  }
                }
                # build incidence matrix
                if( unlist(lapply(prov, class))[randomTermProv2[irandom2]] %in% c("factor","character") ){
                  goodLevels <- intersect(unique(prov[,randomTermProv2[irandom2]]), rownames(M) )
                  if(length(goodLevels) == 0){ # if no match then use the regular model matrix
                    namesZ <- unique(prov[,randomTermProv2[irandom2]])
                    M <- Matrix::Diagonal(n=length(namesZ)); rownames(M) <- colnames(M) <- namesZ
                  }
                  xx = lme4breeding::redmm(x=prov[,randomTermProv2[irandom2]], M=M, nPC=0)
                }else{
                  if(sommerVersion < 44){
                    xx <- sommer::isc(prov[,randomTermProv2[irandom2]])$Z
                  }else{
                    xx <- sommer::ism(prov[,randomTermProv2[irandom2]])$Z
                  }
                }
                xxList[[irandom2]] = xx # model matrix for ith effect saved
                Mlist[[irandom2]] = M # single factor kernel M saved
                # if irandom > 1 expand M
                if(irandom2 > 1){
                  if(ncol(xxList[[irandom2-1]]) > 1){
                    if(sommerVersion < 44){
                      m1 <- sommer::dsc(xxList[[irandom2-1]])
                    }else{
                      m1 <- sommer::dsm(xxList[[irandom2-1]])
                    }
                  }else{
                    if(sommerVersion < 44){
                      m1 <- sommer::isc(xxList[[irandom2-1]][,1])
                    }else{
                      m1 <- sommer::ism(xxList[[irandom2-1]][,1])
                    }

                  }
                  if(ncol(xxList[[irandom2]]) > 1){
                    if(sommerVersion < 44){
                      m2 <- sommer::isc(xx)
                    }else{
                      m2 <- sommer::ism(xx)
                    }

                  }else{
                    if(sommerVersion < 44){
                      m2 <- sommer::isc(xx[,1])
                    }else{
                      m2 <- sommer::ism(xx[,1])
                    }

                  }
                  if(sommerVersion < 44){
                    m3 <- sommer::vsc( m1  , m2  )
                  }else{
                    m3 <- sommer::vsm( m1  , m2  )
                  }

                  environmentCol <- list()
                  for(o in 1:length(m3$Z)){environmentCol[[o]] <- rep(colnames(m3$theta)[o],nrow(M))}
                  ff <- do.call( "cbind", m3$Z )
                  M <- kronecker(Mlist[[irandom2-1]] , M, make.dimnames = TRUE) # update M
                }else{
                  ff <- xxList[[irandom2]]
                  # M <- M
                }
              }
              # compute environment column for later
              namesForEnvs <- lapply(Mlist,function(x){rownames(x)})
              namesForEnvs=do.call(expand.grid, rev(namesForEnvs))
              if(ncol(namesForEnvs)==1){ # if there's no interactions
                envs <- rep("(Intercept)",nrow(M))
              }else{ # if there's interactions
                nLevsInEnvs <- apply(namesForEnvs,2, function(x){length(unique(x))})
                # remove the one with the biggest number of levels
                namesForEnvs <- namesForEnvs[,-c(which(nLevsInEnvs == max(nLevsInEnvs))),drop=FALSE]
                envs <- apply(namesForEnvs,1,function(x){paste(x,collapse = ":")})
              }
              xxList=NULL;Mlist=NULL
              groupingTermProv[[irandom]] <- c( (ncol(prov)+1) : ( ncol(prov)+ncol(ff) ) ) # build grouping term
              prov <- cbind(prov, as.matrix(ff)) # bind matrix to dataset
              Mprov[[irandom]] <- M # save M matrix that combines previous effects
              envsProv[[irandom]] <- envs # save levels for environment
              # compute entry type column for later
              entryTypeProv[[irandom]] <- paste(expCovariatesProv2,collapse = ":") # save info for kernels used in the different effects
            }
          }

          if ("genoD" %in% unlist(expCovariatesProv)) { #rename designation effects
            for (i in seq_along(randomTermProv)) {
              for (j in seq_along(randomTermProv[[i]])) {
                if (randomTermProv[[i]][j] == "designation") {
                  if (expCovariatesProv[[i]][j] == "genoA") {
                    randomTermProv[[i]][j] <- "designationA"
                  } else if (expCovariatesProv[[i]][j] == "genoD") {
                    randomTermProv[[i]][j] <- "designationD"
                  }
                }
              }
            }
          }

          randomTermTrait[[iTrait]] <- unique(randomTermProv) # random formula for the trait
          myDataTraits[[iTrait]] <- prov # dataset for this trait
          names(groupingTermProv) <- names(envsProv) <- names(entryTypeProv) <- names(randomTermTrait[[iTrait]]) <- names(Mprov) <- unlist(lapply(randomTermProv, function(x){paste(x,collapse = "_")}))
          groupingTermTrait[[iTrait]] <- groupingTermProv # grouping for this trait
          Mtrait[[iTrait]] <- Mprov # save the M matrix that combines all single M kernel matrices to later recover the BLUPs
          envsTrait[[iTrait]] <- envsProv # save the values for environment column
          entryTypesTrait[[iTrait]] <- entryTypeProv
          # end of formula formation
        }
      }
    }
  }
  # print(groupingTermTrait)
  ##########################################
  ##########################################
  ## MODEL FITTING
  if(verbose){message("Fitting a model.")}
  predictionsList <- list();
  for(iTrait in names(myDataTraits)){ # # iTrait = trait[1]  iTrait="value"
    if(verbose){message(paste("Analyzing trait", iTrait))}
    mydataSub <- myDataTraits[[iTrait]] # extract dataset
    groupingSub <- groupingTermTrait[[iTrait]] # extract grouping indices
    Msub <- Mtrait[[iTrait]] # extract the M kernel matrices
    envsSub <- envsTrait[[iTrait]] # extract the values for the environment column
    entryTypesSub <- entryTypesTrait[[iTrait]] # extract values for the entryType column (kernels used)
    fixedTermSub <- fixedTermTrait[[iTrait]] # extract fixed formula
    # names(fixedTermSub) <- unlist(lapply(fixedTermSub,function(x){paste(x,collapse = ":")}))
    randomTermSub <- randomTermTrait[[iTrait]] # extract random formula
    ## deregress if needed
    VarFull <- var(mydataSub[,"predictedValue"], na.rm = TRUE) # total variance
    effectTypeTrait <- phenoDTfile$modeling[which(phenoDTfile$modeling$analysisId == analysisId & phenoDTfile$modeling$trait == iTrait & phenoDTfile$modeling$parameter == "designationEffectType"),"value"]
    if(names(sort(table(effectTypeTrait), decreasing = TRUE))[1] == "BLUP"){ # if STA was BLUPs deregress
      mydataSub$predictedValue <- mydataSub$predictedValue/mydataSub$reliability
    }
    ## calculate weights
    mydataSub=mydataSub[with(mydataSub, order(environment)), ] # sort by environments
    mydataSub$w <- 1/(mydataSub$stdError^2) # add weights column
    ## get formula
    fix <- paste( unlist(lapply(fixedTermSub, function(x){paste(x, collapse = ":")})), collapse = " + ")
    fix <- paste("predictedValue ~", fix)
    if(length(groupingSub) > 0){
      ranran <- paste("~", paste(paste0("grp(",names(groupingSub),")"), collapse=" + "))
    }else{ranran <- character()}
    if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
    # warnin messages in weights use
    if(useWeights){
      weightsFormulation="w"
      if(verbose){message("   Using weights in the analysis. Residual variance will be fixed to 1.")  }
    }else{
      weightsFormulation=NULL
      if(verbose){message("   Ignoring weights in the analysis. Residual variance will be estimated.")  }
    }
    if(is.null(randomTermSub)){groupingSub=NULL}
    # print(groupingSub)
    # print(ranFormulation)

    ## model fit
    mix <- try(
      LMMsolver::LMMsolve(fixed =as.formula(fix),
                          random = ranFormulation,
                          # residual=ranres,
                          weights = weightsFormulation,
                          # ginverse = myGinverse,
                          group = groupingSub,
                          family = eval(parse(text = traitFamily[iTrait])),
                          data = mydataSub, maxit = maxIters),
      silent = TRUE
    )
    # print(mix$VarDf)

    pp <- list()
    if(!inherits(mix,"try-error") ){ # if random model runs well try the fixed model
      ## save the modeling used
      currentModeling <- data.frame(module="mtaLmms", analysisId=mtaAnalysisId,trait=iTrait, environment=c(rep("across",3), names(unlist(entryTypesSub))),
                                    parameter=c("fixedFormula","randomFormula","family",rep("kernels",length(unlist(entryTypesSub)))),
                                    value=c(fix,ifelse(length(ranran)>0,ranran,NA),traitFamily[iTrait],unlist(entryTypesSub) ))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
      ## save the environments used goodFields
      currentModeling <- data.frame(module="mtaLmms", analysisId=mtaAnalysisId,trait=iTrait, environment=allEnvironments,
                                    parameter="includedInMta",
                                    value=ifelse(allEnvironments%in%unique(mydataSub$environment), TRUE, FALSE))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
      # get variance components
      ss <- mix$VarDf;  rownames(ss) <- ss$VarComp
      Ve <- ss["residual","Variance"]
      mu <- mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
      Ci <- solve(mix$C)
      if(length(mu) > 0){
        pp[["(Intercept)"]] <- data.frame(designation="(Intercept)", predictedValue=mu, stdError=sqrt(Ci[1,1]), reliability=NA,
                                          trait=iTrait, effectType="(Intercept)", entryType="(Intercept)", environment="(Intercept)" )
      }
      fixedEffects <- setdiff(mix$EDdf$Term, mix$VarDf$VarComp)
      fixedEffects <- setdiff(fixedEffects, "(Intercept)")

      for(iGroupFixed in fixedEffects){ # iGroupFixed = fixedEffects[1]

        pick <- mix$ndxCoefficients[[iGroupFixed]]
        pick <- pick[which(pick!=0)]

        # shouldBeOne <- which(pick == 0)
        # if(length(shouldBeOne) > 0){pick[shouldBeOne] = 1}
        blue <- mix$coefMME[pick] + mu; names(blue) <- names(pick); #blue[1] <- blue[1]-mu
        start <- sum(mix$EDdf[1:(which(mix$EDdf$Term == iGroupFixed) - 1),"Model"]) # we don't add a one because we need the intercept
        nEffects <- mix$EDdf[which(mix$EDdf$Term == iGroupFixed),"Effective"]#length(blue)
        pev <- Ci[start:(start+nEffects-1),start:(start+nEffects-1)]
        if(is.matrix(pev)){ stdError <- (sqrt(Matrix::diag(pev)))}else{stdError <- pev}

        prov <- data.frame(designation=names(blue), predictedValue=blue, stdError=stdError, reliability=NA,
                           trait=iTrait, effectType=iGroupFixed, environment="(Intercept)" )

        for(iLabel in unique(unlist(fixedTermSub))){
          prov$designation <- gsub(paste0(iLabel,"_"),"",prov$designation)
        }

        # add additional entry type labels
        mydataSub[,"designationXXX"] <- apply(mydataSub[,unlist(strsplit(iGroupFixed,":")),drop=FALSE],1,function(x){paste(x,collapse = ":")})
        prov$entryType <- apply(data.frame(prov$designation),1,function(x){
          found <- which(mydataSub[,"designationXXX"] %in% x)
          if(length(found) > 0){
            x2 <- paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#");
          }else{x2 <- "unknown"}
          return(x2)
        })
        prov$entryType <- cgiarBase::replaceValues(prov$entryType, Search = "", Replace = "unknown")

        # save
        pp[[iGroupFixed]] <- prov
      };
      if(!is.null(randomTermSub)){
        for( iGroup in names(groupingSub)){ # iGroup=names(groupingSub)[2]
          pick <- mix$ndxCoefficients[[iGroup]]
          shouldBeOne <- which(pick == 0)
          if(length(shouldBeOne) > 0){pick[shouldBeOne] = 1}
          blup <- (Msub[[iGroup]] %*% mix$coefMME[pick]); blup <- as.vector(blup)
          names(blup) <- rownames(Msub[[iGroup]]) # Msub will always be a matrix wither a diagonal or different than but do it across for consistency
          start <- sum(mix$EDdf[1:(which(mix$EDdf$Term == iGroup) - 1),"Model"]) # we don't add a one because we need the intercept
          nEffects <- ncol(Msub[[iGroup]])
          Vg <- ss[iGroup,"Variance"]
          if(calculateSE){
            if(verbose){message(paste("   Calculating standar errors for",iTrait, iGroup,"predictions"))}
            # Ci <- solve(mix$C)
            Cinv <- Ci[start:(start+nEffects-1),start:(start+nEffects-1)]
            if(is.matrix(Cinv)){ # when there is more than one effect
              Cinv <- as(as(as( Cinv,  "dMatrix"), "generalMatrix"), "CsparseMatrix") # as(Cinv, Class = "dgCMatrix")
              startPev <- seq(1, length(blup), 500)
              endPev  <- c(startPev - 1, length(blup)); endPev <- endPev[-1]
              stdError <- list()
              for(s in 1:length(startPev)){
                use <- (startPev[s]:endPev[s])
                stdError[[s]] <-  sqrt(Matrix::diag( Msub[[iGroup]][use,] %*% Matrix::tcrossprod( Cinv, Msub[[iGroup]][use,]) ) )
              }
              stdError <- unlist(stdError)
            }else{stdError <- Cinv} # random effect was just one column
            reliability <- abs((Vg - (stdError^2)) /Vg) # reliability <- abs((Vg - Matrix::diag(pev))/Vg)
          }else{stdError <- reliability <- rep(NA,length(blup))}
          badRels <- which(reliability > 1); if(length(badRels) > 0){reliability[badRels] <- 0.9999}
          badRels2 <- which(reliability < 0); if(length(badRels2) > 0){reliability[badRels2] <- 0}
          prov <- data.frame(designation=names(blup), predictedValue=blup, stdError=stdError, reliability=reliability,
                             trait=iTrait, effectType=iGroup , environment=envsSub[[iGroup]] )
          # add fixed effects if present in the random term
          feToAdd <- intersect( randomTermSub[[iGroup]], fixedEffects ) # unlist(fixedTermSub)
          if(length(feToAdd) > 0){
            varInppGroup <- strsplit( prov[,"designation"], ":")
            for(iFe in feToAdd){ # iFe = feToAdd[1]
              provFe <- pp[[iFe]]
              rownames(provFe) <- gsub(paste0(iFe,"_"),"", provFe[,"designation"])
              pickVarInppGroup <- which(randomTermSub[[iGroup]] == feToAdd)
              feUsed <- unlist(lapply(varInppGroup, function(x){x[pickVarInppGroup]}))
              mu0 <- provFe[feUsed,"predictedValue"]; mu0[which(is.na(mu0))]=mu
              prov[,"predictedValue"] <-  prov[,"predictedValue"] + mu0
            }
          }else{
            prov[,"predictedValue"] <-  prov[,"predictedValue"] + mu
          }
          # end of adding fixed effects
          sdP <- sd(prov[,"predictedValue"],na.rm=TRUE)
          cv <- (sd(prov[,"predictedValue"],na.rm=TRUE)/mean(prov[,"predictedValue"],na.rm=TRUE))*100
          # add additional entry type labels
          colsToUse <- unlist(randomTermSub[[iGroup]])
          colsToUse[colsToUse %in% c("designationA", "designationD")] <- "designation"
          mydataSub[,"designationXXX"] <- apply(mydataSub[,colsToUse,drop=FALSE],1,function(x){paste(x,collapse = ":")})
          prov$entryType <- apply(data.frame(prov$designation),1,function(x){
            found <- which(mydataSub[,"designationXXX"] %in% x)
            if(length(found) > 0){
              x2 <- paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#");
            }else{x2 <- "unknown"}
            return(x2)
          })
          prov$entryType <- cgiarBase::replaceValues(prov$entryType, Search = "", Replace = "unknown")
          # save
          pp[[iGroup]] <- prov
          phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                       data.frame(module="mtaLmms",analysisId=mtaAnalysisId, trait= iTrait,
                                                  environment=paste(unique(envsSub[[iGroup]]), collapse = "_"),
                                                  parameter=c( paste(c("mean","sd", "r2","Var"),iGroup,sep="_") ),
                                                  method=c("sum(x)/n","sd","(G-PEV)/G","REML"),
                                                  value=c(mean(prov[,"predictedValue"], na.rm=TRUE), sdP, median(reliability), var(prov[,"predictedValue"], na.rm=TRUE) ),
                                                  stdError=c(NA,NA,sd(reliability, na.rm = TRUE)/sqrt(length(reliability)),NA )
                                       )
          )
        }
      }
    }else{ # if model failed
      if(verbose){ cat(paste("Mixed model failed for trait",iTrait,". Aggregating and assuming h2 = 0 \n"))}
      means <- aggregate(predictedValue ~ designation, FUN=mean, data=mydataSub)
      Ve <- var(mydataSub[,"predictedValue"], na.rm=TRUE)
      means$environment <- "(Intercept)"
      means$stdError <- sd(means$predictedValue)
      means$reliability <- 1e-6
      means$trait <- iTrait
      means$effectType <- "designation"
      means$entryType <- "unknown"
      sdP <- sd(means$predictedValue,na.rm=TRUE)
      cv <- (sd(means$predictedValue,na.rm=TRUE)/mean(means$predictedValue,na.rm=TRUE))*100
      ## save metrics
      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                   data.frame(module="mtaLmms",analysisId=mtaAnalysisId, trait=iTrait,
                                              environment="across",
                                              parameter=c("mean","sd", "r2","Var_designation","Var_residual"),
                                              method=c("sum(x)/n","sd","(G-PEV)/G","REML","REML"),
                                              value=c(mean(means$predictedValue, na.rm=TRUE), sdP, 0, 0, Ve ),
                                              stdError=NA
                                   )
      )
      currentModeling <- data.frame(module="mtaLmms", analysisId=mtaAnalysisId,trait=iTrait, environment="across",
                                    parameter=c("fixedFormula","randomFormula","family","designationEffectType"),
                                    value=c("None","None","None","mean"))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
      pp[["designation"]] <- means
    }

    predictionsTrait <- do.call(rbind,pp)
    #############################################################
    ## add across env estimate for DESIGNATION effect type fitted
    #############################################################
    match1 <- unlist(lapply(fixedTermSub,function(x){sum(as.numeric(x=="designation"))}))
    names(match1) <- unlist(lapply(fixedTermSub,function(x){paste(x,collapse = "_")}))
    match2 <- unlist(lapply(randomTermSub,function(x){sum(as.numeric(x=="designation"))}))
    match3 <- c(match1,match2)
    useForPreds <- names(match3)[which(match3 > 0)]
    doublematch <- table(predictionsTrait$effectType, predictionsTrait$environment)
    rownames(doublematch) <- gsub(":", "_", rownames(doublematch) )
    interceptCheck <- sum(apply(data.frame(useForPreds),1,function(x){
      if(x %in% rownames(doublematch)){
        return(
          sum(as.numeric("(Intercept)" %in% colnames(doublematch)[which(doublematch[x,]>0)] ))
        )
      }else{ return(0) }
    }))
    # interceptCheck <- sum(apply(data.frame(useForPreds),1,function(x){sum(as.numeric("(Intercept)" %in% colnames(doublematch)[which(doublematch[x,]>0)] ))}))
    '%!in%' <- function(x,y)!('%in%'(x,y))
    if( length(useForPreds) > 0 & interceptCheck==0 ){ # only if there was designation and no main effect exist then we aggregate
      provx <- predictionsTrait
      provx <- provx[which(provx$effectType %in% useForPreds),]
      provx$designation <- apply(provx[,c("environment","designation")],1,function(x){gsub(paste0(x[1],":"),"",x[2])})
      provx <- aggregate(cbind(predictedValue,stdError,reliability)~designation+trait, FUN=mean, data=provx)
      provx$environment <- "(Intercept)"
      provx$effectType <- "designation"
      provx$entryType <- apply(data.frame(provx$designation),1,function(x){
        found <- which(mydataSub[,"designationXXX"] %in% x)
        if(length(found) > 0){
          x2 <- paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#");
        }else{x2 <- "unknown"}
        return(x2)
      })
      provx$entryType <- cgiarBase::replaceValues(provx$entryType, Search = "", Replace = "unknown")
      predictionsTrait <- rbind(predictionsTrait, provx[,colnames(predictionsTrait)])
    }
    #############################################################
    ## add across env estimate for GID effect type fitted
    #############################################################
    match1 <- unlist(lapply(fixedTermSub,function(x){sum(as.numeric(x=="gid"))}))
    names(match1) <- unlist(lapply(fixedTermSub,function(x){paste(x,collapse = "_")}))
    match2 <- unlist(lapply(randomTermSub,function(x){sum(as.numeric(x=="gid"))}))
    match3 <- c(match1,match2)
    useForPreds <- names(match3)[which(match3 > 0)]
    doublematch <- table(predictionsTrait$effectType, predictionsTrait$environment)
    interceptCheck <- sum(apply(data.frame(useForPreds),1,function(x){sum(as.numeric("(Intercept)" %in% colnames(doublematch)[which(doublematch[x,]>0)]))}))
    '%!in%' <- function(x,y)!('%in%'(x,y))
    if( length(useForPreds) > 0 & interceptCheck==0 ){ # only if there was designation and no main effect exist then we aggregate
      provx <- predictionsTrait
      provx <- provx[which(provx$effectType %in% useForPreds),]
      provx$designation <- apply(provx[,c("environment","designation")],1,function(x){gsub(paste0(x[1],":"),"",x[2])})
      provx <- aggregate(cbind(predictedValue,stdError,reliability)~designation+trait, FUN=mean, data=provx)
      provx$environment <- "(Intercept)"
      provx$effectType <- "gid"
      provx$entryType <- apply(data.frame(provx$designation),1,function(x){
        found <- which(mydataSub[,"designationXXX"] %in% x)
        if(length(found) > 0){
          x2 <- paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#");
        }else{x2 <- "unknown"}
        return(x2)
      })
      provx$entryType <- cgiarBase::replaceValues(provx$entryType, Search = "", Replace = "unknown")
      predictionsTrait <- rbind(predictionsTrait, provx[,colnames(predictionsTrait)])
    }
    #
    predSta <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId &
                                               phenoDTfile$predictions$trait == iTrait &
                                               phenoDTfile$predictions$environment %in% envCount[[iTrait]]),]
    phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                 data.frame(module="mtaLmms",analysisId=mtaAnalysisId, trait= iTrait, environment="across",
                                            parameter=c("Var_residual","nEnv","nEntries"),
                                            method=c("REML","n","n"),
                                            value=c( Ve, length(envCount[[iTrait]]), length(unique(predSta$designation)) ),
                                            stdError=c(NA,NA,NA) )
    )
    predictionsList[[iTrait]] <- predictionsTrait
  }
  ## enf of model fitting
  if(length(predictionsList) == 0){stop("There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all envs.",call. = FALSE)}
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- mtaAnalysisId

  ##########################################
  ## add timePoint of origin, stage and designation code
  if(verbose){message("Wrapping the results.")}
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
  predictionsBind$module <- "mtaLmms"; rownames(predictionsBind) <- NULL

  #print(head(predictionsBind))
  #########################################
  ## update databases
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(!is.null(phenoDTfile$predictions)){
    if("effectType" %!in% colnames(phenoDTfile$predictions) ){
      phenoDTfile$predictions$effectType <- NA
    }
  }

  ##Adapt exports to GPCP model

  if(all(c("designationA", "designationD") %in% predictionsBind$effectType)){
    desA <- predictionsBind[predictionsBind$effectType == "designationA", ]
    desD <- predictionsBind[predictionsBind$effectType == "designationD", ]

    # Make sure rows align by designation + trait
    keyCols <- c("designation", "trait")
    desA <- desA[order(desA[[keyCols[1]]], desA[[keyCols[2]]]), ]
    desD <- desD[order(desD[[keyCols[1]]], desD[[keyCols[2]]]), ]

    #Get averages
    avgDes <- desA
    avgDes$predictedValue <- rowMeans(cbind(desA$predictedValue, desD$predictedValue), na.rm = TRUE)
    avgDes$stdError <- rowMeans(cbind(desA$stdError, desD$stdError), na.rm = TRUE)
    avgDes$effectType <- "designation"

    # Add the averaged designation rows
    predictionsBind <- rbind(predictionsBind, avgDes)
  }


  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,
                                   predictionsBind[,colnames(phenoDTfile$predictions)])



  newStatus <- data.frame(module="mtaLmms", analysisId=mtaAnalysisId, analysisIdName=NA)
  phenoDTfile$status <- rbind( phenoDTfile$status, newStatus[,colnames(phenoDTfile$status)] )
  ## add which data was used as input
  modeling <- data.frame(module="mtaLmms",  analysisId=mtaAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId ))
  if(!is.null(nPC)){
    modeling <- rbind(modeling,
                      data.frame(module="mtaLmms",  analysisId=mtaAnalysisId, trait=names(nPC), environment="general",
                                 parameter= c("nPC"), value= nPC )
    )
  }
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)
}
