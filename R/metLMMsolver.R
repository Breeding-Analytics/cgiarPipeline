metLMMsolver <- function(
    phenoDTfile= NULL, analysisId=NULL, analysisIdGeno = NULL,
    fixedTerm= list("1"),  randomTerm=NULL, expCovariates=NULL,
    envsToInclude=NULL, trait= NULL, traitFamily=NULL,
    useWeights=TRUE,estHybrids = TRUE,
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
  
  .tokens_from_terms <- function(termList) {
    if (is.null(termList) || length(termList) == 0) return(character())
    out <- unlist(lapply(termList, function(v) unlist(strsplit(v, ":", fixed = TRUE))), use.names = FALSE)
    unique(out[nzchar(out)])
  }
  
  ## configurable safety limits 
  MAX_DENSE_BYTES <- getOption("bioflow.max_dense_bytes", 10 * 1024^3)  # ~6 GB
  MAX_PROJ_COLS   <- getOption("bioflow.max_proj_cols",   1e5)    # 100k columns
  
  ## helper: check projected size and stop 
  .check_cross_guard <- function(n_rows, n_cols, tag) {
    est_bytes <- as.double(n_rows) * as.double(n_cols) * 8  # double matrix
    if (!is.finite(est_bytes) || n_cols > MAX_PROJ_COLS || est_bytes > MAX_DENSE_BYTES) {
      stop(sprintf(
        paste0("Model too large (", tag, "). ",
               "Projected crossed design has ~%s rows x %s cols (≈ %.2f GB). ",
               "Reduce factors in interaction terms or try the model without covariates"),
        format(n_rows, big.mark=","), format(n_cols, big.mark=","), est_bytes/1073741824
      ), call. = FALSE)
    }
  }
  
  #Helpers to get the levels of a fitted factor
  .canon_term <- function(s, sep=":") paste(sort(strsplit(s, sep, fixed=TRUE)[[1]]), collapse=":")
  
  .find_ndx_key <- function(keys, query) {
    want <- .canon_term(query, ":")
    hit <- vapply(keys, function(k) .canon_term(k, ":") == want, logical(1))
    if (!any(hit)) return(NA_character_)
    keys[which(hit)[1]]
  }
  
  .levels_for_term <- function(dat, term_label) {
    # Build design matrix for this term
    ff <- as.formula(paste0("~ 0 + ", term_label))
    mm <- model.matrix(ff, data = dat)
    cn <- colnames(mm)
    
    # Strip variable-name prefixes from the column names
    # e.g. "designation2369:environmentNYH3y2021" -> "2369:NYH3y2021"
    vars <- strsplit(term_label, ":", fixed = TRUE)[[1]]
    for (v in vars) {
      cn <- sub(paste0("^", v), "", cn)   # at the start
      cn <- sub(paste0(":", v), ":", cn)  # after a colon
    }
    
    cn
  }
  
  #helper for models without covariates
  all_none_covariates <- function(expCovariates) {
    x <- expCovariates
    # normalize: lower, trim, strip trailing dots
    x <- tolower(trimws(x))
    x <- sub("\\.*$", "", x)   # remove any number of trailing dots
    x <- x[nzchar(x)]          # drop empty after normalization
    ok <- length(x) == 0 || all(x == "none")
    ok
  }
  
  #Flags for SCA and GCA models
  is_SCA_GCA = FALSE
  is_GCA = FALSE
  
  find_sca_gca = intersect(unlist(randomTerm),c("mother","father","designation"))
  if (setequal(find_sca_gca, c("mother","father","designation"))){
    is_SCA_GCA = TRUE
  }else if(setequal(find_sca_gca, c("mother","father"))){
    is_GCA = TRUE
  }else if(setequal(find_sca_gca, c("mother"))){
    is_GCA = TRUE
  }else if(setequal(find_sca_gca, c("father"))){
    is_GCA = TRUE
  }
  
  
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
    Weather <- apply(Weather,2,enhancer::imputev)
    colnames(Weather) <- gsub(" ","",colnames(Weather))
  }
  covars <- unique(unlist(expCovariates))
  randomTermForCovars <- unique(unlist(randomTerm))
  if(!is.null(randomTermForCovars)){
    if( any( covars %in% c("genoA","genoAD", "weather","pedigree", traitsForExpCovariates ) ) ){
      if(verbose){message("Checking and calculating kernels requested")}
      ## MARKER KERNEL
      Markers <- NULL
      if(!is.null(phenoDTfile$data$geno)){
        Markers <- as.matrix(phenoDTfile$data$geno) # in form of covariates
      }
       
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
        #   # if(length(which(is.na(Markers))) > 0){Markers <- apply(Markers,2,enhancer::imputev)}
        # }else{
        #   missing <- apply(Markers,2,enhancer::propMissing)
        #   Markers <- apply(Markers[,which(missing < 0.9)],2,enhancer::imputev)
        # }
        if(nPC["geno"] < 0){ # do not include extra individuals
          mydataX <-  phenoDTfile$predictions[which( phenoDTfile$predictions$analysisId %in% analysisId),]
          if(!is_SCA_GCA & !is_GCA){
            Markers <- Markers[which(rownames(Markers) %in% unique(mydataX$designation) ), ]
          }
          if(verbose){message(paste("Subsetting marker to",nrow(Markers),"individuals present"))}
        }
        
        ploidyFactor <- max(Markers)/2
        
        if (is_SCA_GCA) {
          #This is the SCA/GCA model
          if(any(unlist(expCovariates) %in% c("genoA","genoD"))){
            #Only records with parental information will be kept
            myData = phenoDTfile$predictions
            myData <- myData[which(myData$analysisId %in% analysisId),]
            myData = myData[!is.na(myData$mother),]
            myData = myData[!is.na(myData$father),]
            
            has_genoMother = myData$mother %in% rownames(Markers)
            has_genoFather = myData$father %in% rownames(Markers)
            
            to_model = myData[has_genoMother & has_genoFather,]
            if(nrow(to_model)==0){
              stop("Parental genotypes are required to run this model")
            }
            hybrids_for_sca = unique(to_model$designation)
            mothers_for_sca = unique(to_model$mother)
            fathers_for_sca = unique(to_model$father)
            
            if(ploidyFactor == 1 & estHybrids){ #we can estimate hybrid diploid genotypes from parental lines
              desig_noGeno = myData[!myData$designation %in% rownames(Markers),]
              
              if(nrow(desig_noGeno)>0){
                has_genoMother = desig_noGeno$mother %in% rownames(Markers)
                has_genoFather = desig_noGeno$father %in% rownames(Markers)
                
                desig_noGeno = desig_noGeno[has_genoMother & has_genoFather,]
                desig_noGeno = unique(desig_noGeno[,c("designation","mother","father")])
                if(nrow(desig_noGeno)>0){
                  n_hybrids = nrow(desig_noGeno)
                  message("Parental genotypes found for: ",n_hybrids," hybrids \n")
                  
                  M_mother = Markers[desig_noGeno$mother,]
                  M_father = Markers[desig_noGeno$father,]
                  M_hybrid = (M_mother + M_father)/2 
                  rownames(M_hybrid) = desig_noGeno$designation
                  
                  Markers = rbind(Markers,M_hybrid)
                  
                }else{
                  message("No parental genotypes found for missing hybrids.")
                }
              }
            }
          }
        }else if(is_GCA){
          
          if(any(unlist(expCovariates) == "genoA")){
          #Only records with parental information will be kept
          myData = phenoDTfile$predictions
          myData = myData[!is.na(myData$mother),]
          myData = myData[!is.na(myData$father),]
          
          has_genoMother = myData$mother %in% rownames(Markers)
          has_genoFather = myData$father %in% rownames(Markers)
          
          to_model = myData[has_genoMother & has_genoFather,]
          mothers_for_gca = unique(to_model$mother)
          fathers_for_gca = unique(to_model$father)
          }
        }else if(setequal(find_sca_gca, c("mother"))){
          myData = phenoDTfile$predictions
          myData = myData[!is.na(myData$mother),]
          
          has_genoMother = myData$mother %in% rownames(Markers)
          to_model = myData[has_genoMother,]
          mothers_for_gca = unique(to_model$mother)
          
        }else if(setequal(find_sca_gca, c("father"))){
          myData = phenoDTfile$predictions
          myData = myData[!is.na(myData$father),]
          
          has_genoFather = myData$father %in% rownames(Markers)
          to_model = myData[has_genoFather,]
          fathers_for_gca = unique(to_model$father)
        }  
        
        
        if("genoA" %in% covars){G <- sommer::A.mat(Markers-ploidyFactor);} # additive model
        if("genoD" %in% covars){ #Dominance kernel
          if(ploidyFactor == 1){

            D <- 1-abs(Markers)

            f <- rowSums(D) / ncol(D) #inbreeding fixed eff
            names(f) <- rownames(Markers)

            D <- sommer::D.mat(Markers-ploidyFactor)
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
                wideTrait <- apply(wideTrait[,-1],2,enhancer::imputev)
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
      
      ## Enforce factor vs covariate semantics  ---
      ## Anything selected in Fixed/Random menus => factor,
      ## EXCEPT explicit covariates (right-side boxes) and a tiny numeric whitelist.
      
      # What the user asked for on the left (Fixed/Random)
      left_fixed_vars  <- .tokens_from_terms(fixedTerm)         # before per-trait pruning
      left_random_vars <- .tokens_from_terms(randomTerm)
      
      # What the user explicitly set as covariates on the right
      right_covars <- unique(unlist(expCovariates))  # e.g., "genoA","weather","pedigree","envIndex...", trait kernels, etc.
      
      # Numeric fixed terms you intend to stay numeric (add/remove as needed)
      numeric_fixed_whitelist <- c("inbreeding", grep("^envIndex", colnames(prov), value = TRUE))
      
      # Columns that MUST be treated as categorical (factor) if they appear in prov
      must_be_factor <- setdiff(unique(c(left_fixed_vars, left_random_vars)),
                                unique(c(right_covars, numeric_fixed_whitelist)))
      
      # Coerce present columns in 'prov' to factor
      to_factor <- intersect(must_be_factor, colnames(prov))
      for (vv in to_factor) {
        if (!is.factor(prov[[vv]])) {
          prov[[vv]] <- factor(trimws(as.character(prov[[vv]])))
        } else {
          prov[[vv]] <- droplevels(prov[[vv]])
        }
      }

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
          
          use_formula = all_none_covariates(unlist(expCovariatesProv))
          
          #################################################################################
          #When no covariates are used, we can just pass the formulas directly to LMMsolver
          
          if(use_formula){
            # build envsProv for every random term in randomTermProv
            envsProv <- lapply(randomTermProv, function(term) {
              term <- term[term %in% colnames(prov)]  # drop any names not in prov
              if (length(term) <= 1L) {
                # single (or empty) factor: repeat "(Intercept)" once per level of that factor
                if (length(term) == 1L) {
                  levs1 <- sort(unique(as.character(prov[[term[1]]])))
                  rep("(Intercept)", length(levs1))
                } else {
                  "(Intercept)"
                }
              } else {
                # interaction: build level lists, drop the dimension with the most levels, expand/paste the rest
                levlist <- lapply(term, function(v) sort(unique(as.character(prov[[v]]))))
                nlev    <- vapply(levlist, length, integer(1))
                drop    <- which.max(nlev)
                keep    <- setdiff(seq_along(levlist), drop)
                if (length(keep) == 0L) {
                  rep("(Intercept)", nlev[drop])
                } else {
                  grid <- do.call(expand.grid, c(rev(levlist[keep]), stringsAsFactors = FALSE))
                  apply(grid, 1L, function(x) paste(x, collapse = ":"))
                }
              }
            })
            
           
            randomTermTrait[[iTrait]] <- unique(randomTermProv[which(unlist(lapply(randomTermProv, length)) > 0)])
            entryTypeProv <- list()
            
            for(irandom in 1:length(randomTermProv)){
              expCovariatesProv2 <- expCovariatesProv[[irandom]]
              entryTypeProv[[irandom]] <- paste(expCovariatesProv2,collapse = ":")
            }
            
            names(envsProv) <- names(entryTypeProv) <- names(randomTermTrait[[iTrait]]) <- vapply(randomTermProv, function(x) paste(x, collapse = "_"), character(1))
            entryTypesTrait[[iTrait]] <- entryTypeProv
          }else{## We need to pass the full incidence matrices to LMMsolver
          ##################################################################
            if(!is.null(randomTermProv)){
              # reduce datasets
              for(irandom in 1:length(randomTermProv)){ # for each element in the list # irandom=1
                randomTermProv2 <- randomTermProv[[irandom]]
                for(irandom2 in  1:length(randomTermProv2)){ # for each factor in the interactions # irandom2 = 1
                  # print(expCovariatesProv[[irandom]][irandom2])
                  if( expCovariatesProv[[irandom]][irandom2] == "weather"){
                    M <- Wchol
                  }else if(expCovariatesProv[[irandom]][irandom2] %in% c("geno","genoA","genoAD") ){
                    if(is_SCA_GCA){
                      if(randomTermProv[[irandom]][irandom2]=="mother"){
                        M <- Gchol[mothers_for_sca,]
                      }else if(randomTermProv[[irandom]][irandom2]=="father"){
                        M <- Gchol[fathers_for_sca,]
                      }
                    }else if(is_GCA){
                      if(randomTermProv[[irandom]][irandom2]=="mother"){
                        M <- Gchol[mothers_for_gca,]
                      }else if(randomTermProv[[irandom]][irandom2]=="father"){
                        M <- Gchol[fathers_for_gca,]
                      }
                    }else{
                      M <- Gchol
                    }
                  }else if(expCovariatesProv[[irandom]][irandom2] == "genoD"){
                    if(is_SCA_GCA){
                      if(randomTermProv[[irandom]][irandom2]=="designation"){
                        M <- Dchol[hybrids_for_sca,]
                      }else{
                        stop("SCA requires genoD to be fitted for designation")
                      }
                    }else{
                      M <- Dchol
                    }
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
                # DEBUG: sizes and names
                for(irandom2 in  1:length(randomTermProv2)){ # for each factor in the interactions # irandom2 = 2
                  # get kernel
                  if( expCovariatesProv[[irandom]][irandom2] == "weather"){
                    M <- Wchol # Weather
                  }else if(expCovariatesProv[[irandom]][irandom2] %in% c("geno","genoA","genoAD") ){
                    if(is_SCA_GCA){
                      if(randomTermProv[[irandom]][irandom2]=="mother"){
                        M <- Gchol[mothers_for_sca,]
                      }else if(randomTermProv[[irandom]][irandom2]=="father"){
                        M <- Gchol[fathers_for_sca,]
                      }
                    }else if(is_GCA){
                      if(randomTermProv[[irandom]][irandom2]=="mother"){
                        M <- Gchol[mothers_for_gca,]
                      }else if(randomTermProv[[irandom]][irandom2]=="father"){
                        M <- Gchol[fathers_for_gca,]
                      }
                    }else{
                      M <- Gchol
                    }
                  }else if(expCovariatesProv[[irandom]][irandom2] == "genoD"){
                    if(is_SCA_GCA){
                      if(randomTermProv[[irandom]][irandom2]=="designation"){
                        M <- Dchol[hybrids_for_sca,]
                      }
                    }else{
                      M <- Dchol
                    }
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
                    term_name <- randomTermProv2[irandom2]
                    
                    # 1) Make the grouping a clean factor with only the levels that are present
                    xf <- droplevels(factor(as.character(prov[[term_name]])))
                    
                    # 2) Align the kernel rows to those (dropped) levels
                    M  <- M[levels(xf), , drop = FALSE]
                    
                    # 3) Call redmm
                    xx <- enhancer::redmm(x = xf, M = M, nPC = 0)
                  }else{
                    if(sommerVersion < 44){
                      xx <- sommer::isc(prov[,randomTermProv2[irandom2]])$Z
                    }else{
                      xx <- sommer::ism(prov[,randomTermProv2[irandom2]])$Z
                    }
                  }
                  
                  xxList[[irandom2]] = xx # model matrix for ith effect saved
                  Mlist[[irandom2]] = M # single factor kernel M saved
                  
                }
                
                
                make_design <- function(xx) {
                  if (ncol(xx) > 1) {
                    if (sommerVersion < 44) sommer::dsc(xx) else sommer::dsm(xx)
                  } else {
                    if (sommerVersion < 44) sommer::isc(xx[,1]) else sommer::ism(xx[,1])
                  }
                }
                
                if (length(xxList) >= 2) {
                  
                  n_rows <- nrow(xxList[[1]])
                  ncols_each <- vapply(xxList, ncol, integer(1))
                  proj_total  <- prod(ncols_each)
                  .check_cross_guard(n_rows, proj_total, "initial block")
                  
                  m_prev <- make_design(xxList[[1]])
                  for (k in 2:length(xxList)) {
                    m_cur <- make_design(xxList[[k]])
                    m_prev <- if (sommerVersion < 44) sommer::vsc(m_prev, m_cur) else sommer::vsm(m_prev, m_cur)
                    if(k < length(xxList)){
                      if (!is.list(m_prev$Z) || any(vapply(m_prev$Z, is.null, TRUE))) {
                        stop("m_prev$Z malformed (NULL or not a list) before cbind", call. = FALSE)
                      }
                      rows <- vapply(m_prev$Z, function(Z) nrow(Z), integer(1))
                      if (length(unique(rows)) != 1L) {
                        stop("Incompatible row counts in crossed design components before cbind", call. = FALSE)
                      }
                      m_prev <- tryCatch({
                        make_design(as.matrix(do.call("cbind", m_prev$Z)))
                      }, error = function(e) {
                        stop(sprintf(
                          paste0("Insufficient memory while materializing crossed design (after %d factors). ",
                                 "Projected size ~%s rows x %s cols. Original error: %s"),
                          k, format(n_rows, big.mark=","), format(proj_cols_k, big.mark=","), conditionMessage(e)
                        ), call. = FALSE)
                      })
                    }
                  }
                  ff <- do.call("cbind", m_prev$Z)  
                }else{
                  ff <- xxList[[1]]
                }
                
                Mlist_used <- vector("list", length(xxList))
                for (jj in seq_along(xxList)) {
                  xxj  <- xxList[[jj]]
                  Mj   <- Mlist[[jj]]
                  pick <- match(colnames(xxj), colnames(Mj))
                  if (anyNA(pick)) {
                    bad <- which(is.na(pick))
                    stop(sprintf("Column mismatch: %d/%d columns in model matrix not found in kernel for factor %s",
                                 length(bad), ncol(xxj), randomTermProv2[jj]),
                         call. = FALSE)
                  }
                  Mlist_used[[jj]] <- Mj[, pick, drop = FALSE]
                }
                
                Macc <- Reduce(function(A,B) kronecker(A, B, make.dimnames = TRUE), Mlist_used)
                
                # safety check: columns must match ff so Macc %*% coef works
                if (ncol(Macc) != ncol(ff)) {
                  stop(sprintf("Internal mismatch: ncol(Macc)=%d, ncol(ff)=%d for random term %s",
                               ncol(Macc), ncol(ff), paste(randomTermProv2, collapse=":")))
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
                Mprov[[irandom]] <- Macc # save M matrix that combines previous effects
                envsProv[[irandom]] <- envs # save levels for environment
                # compute entry type column for later
                entryTypeProv[[irandom]] <- paste(expCovariatesProv2,collapse = ":") # save info for kernels used in the different effects
              }
            }
            
            if (("genoD" %in% unlist(expCovariatesProv))&(!is_SCA_GCA)) { #rename designation effects
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
            names(groupingTermProv) <- names(envsProv) <- names(entryTypeProv) <- names(randomTermTrait[[iTrait]]) <- names(Mprov) <- unlist(lapply(randomTermProv, function(x){paste(x,collapse = "_")}))
            groupingTermTrait[[iTrait]] <- groupingTermProv # grouping for this trait
            Mtrait[[iTrait]] <- Mprov # save the M matrix that combines all single M kernel matrices to later recover the BLUPs
            entryTypesTrait[[iTrait]] <- entryTypeProv
          }
          
          myDataTraits[[iTrait]] <- prov # dataset for this trait
          envsTrait[[iTrait]] <- envsProv # save the values for environment column
          
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
    entryTypesSub <- entryTypesTrait[[iTrait]] # extract values for the entryType column (kernels used)
    
    if(!use_formula){
      groupingSub <- groupingTermTrait[[iTrait]] # extract grouping indices
      Msub <- Mtrait[[iTrait]] # extract the M kernel matrices
    }
    
    envsSub <- envsTrait[[iTrait]] # extract the values for the environment column
    fixedTermSub <- fixedTermTrait[[iTrait]] # extract fixed formula
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
    
    if(use_formula){
      
      if (length(randomTermSub)) {
        ranran <- paste(vapply(randomTermSub, function(x) paste(x, collapse=":"), ""),
                         collapse = " + ")
        ranran <- paste0("~",ranran)
        ranFormulation <- as.formula(ranran)
      } else {
        ranFormulation <- NULL
      }
      
    }else{
      
      if(length(groupingSub) > 0){
        ranran <- paste("~", paste(paste0("grp(",names(groupingSub),")"), collapse=" + "))
      }else{ranran <- character()}
      if(length(ranran) == 0){ranFormulation=NULL}else{ranFormulation=as.formula(ranran)}
      if(is.null(randomTermSub)){groupingSub=NULL}
      
    }
    
    
    # warnin messages in weights use
    if(useWeights){
      weightsFormulation="w"
      if(verbose){message("   Using weights in the analysis. Residual variance will be fixed to 1.")  }
    }else{
      weightsFormulation=NULL
      if(verbose){message("   Ignoring weights in the analysis. Residual variance will be estimated.")  }
    }
    
    # print(groupingSub)
    # print(ranFormulation)

    ## model fit

    if(use_formula){
      mix <- try(
        LMMsolver::LMMsolve(fixed =as.formula(fix),
                            random = ranFormulation,
                            weights = weightsFormulation,
                            family = eval(parse(text = traitFamily[iTrait])),
                            data = mydataSub, maxit = maxIters),
        silent = TRUE
      )
    }else{
      mix <- try(
        LMMsolver::LMMsolve(fixed =as.formula(fix),
                            random = ranFormulation,
                            weights = weightsFormulation,
                            group = groupingSub,
                            family = eval(parse(text = traitFamily[iTrait])),
                            data = mydataSub, maxit = maxIters),
        silent = TRUE
      )
    }
    
    # print(mix$VarDf)

    pp <- list()
    if(!inherits(mix,"try-error") ){ 
      
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
      
      
      
      ss <- mix$VarDf;  rownames(ss) <- ss$VarComp
      Ve <- ss["residual","Variance"]
      mu <- mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
      #Get coefficient matrix
      C_sp  <- mix$C 
      
      #helper function for inverse trick
      rhs_eye_spam <- function(n, cols) {
        if (length(cols) == 0L)
          return(spam::spam(0, nrow = n, ncol = 0))
        # i = row indices, j = col indices (1..k), entries = 1
        spam::spam(list(i = cols, j = seq_along(cols), entries = rep(1, length(cols))),
                   nrow = n, ncol = length(cols))
      }
      
      if (length(mu) > 0) {
        # 1) find intercept coefficient position in the mixed-model equations
        idx0 <- as.integer(mix$ndxCoefficients$`(Intercept)`)
        if (length(idx0) == 1L && !is.na(idx0) && idx0 > 0) {
          C_sp <- mix$C 
          nC   <- nrow(C_sp)
          E    <- rhs_eye_spam(nC, idx0)   
          
          #solve C %*% X = E  => X = C^{-1} E
          X    <- spam::solve(C_sp, E)
          
          var_mu <- if (is.null(dim(X))) {
            # vector
            as.numeric(X[idx0])
          } else {
            # 2D (base or spam)
            as.numeric(X[idx0, 1])
          }
          
          se_mu <- sqrt(max(var_mu, 0))
          }
        
        pp[["(Intercept)"]] <- data.frame(
          designation   = "(Intercept)",
          predictedValue= mu,
          stdError      = se_mu,
          reliability   = NA,
          trait         = iTrait,
          effectType    = "(Intercept)",
          entryType     = "(Intercept)",
          environment   = "(Intercept)"
        )
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
        nEffects <- length(pick)
        # indices of this fixed term within C
        idx   <-  start:(start + nEffects - 1L) 
        nC    <- nrow(mix$C)
        
        chunk <- 400L
        dvals <- numeric(length(idx))
        pos   <- 1L
        while (pos <= length(idx)) {
          cols <- idx[pos:min(pos + chunk - 1L, length(idx))]
          E    <- rhs_eye_spam(nC, cols)                      # nC x k (very skinny)
          X    <- spam::solve(C_sp, E)                        # solves C %*% X = E
          # each needed diagonal element is X[cols[j], j]
          k    <- length(cols)
          dvals[pos:(pos + k - 1L)] <- X[cbind(cols, seq_len(k))]
          pos  <- pos + k
        }
        stdError    <- sqrt(pmax(dvals, 0))

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
        if(use_formula){
          rnd_labels <- vapply(randomTermSub, function(x) paste(x, collapse="_"), "")
        }else{
          rnd_labels <- names(groupingSub)
        }
        
        for( iGroup in rnd_labels){ # iGroup=names(groupingSub)[2]
          
          iGroup_colon <- gsub("_", ":", iGroup, fixed = TRUE)
          
          if(use_formula){
            ndx_keys <- names(mix$ndxCoefficients)
            key <- .find_ndx_key(ndx_keys, iGroup_colon)
            pick <- mix$ndxCoefficients[[key]]
          }else{
            pick <- mix$ndxCoefficients[[iGroup]]
          }
          
          shouldBeOne <- which(pick == 0)
          if(length(shouldBeOne) > 0){pick[shouldBeOne] = 1}
          
          if(use_formula){
            blup <- as.vector(mix$coefMME[pick])
            names(blup) <- .levels_for_term(mydataSub, iGroup_colon)
            nEffects <- length(pick)
          }else{
            blup <- (Msub[[iGroup]] %*% mix$coefMME[pick]); blup <- as.vector(blup)
            names(blup) <- rownames(Msub[[iGroup]]) 
            nEffects <- ncol(Msub[[iGroup]])
          }
          
          if(use_formula){
            ed_keys <- mix$EDdf$Term
            key_ed  <- .find_ndx_key(ed_keys, iGroup_colon)
            start <- sum(mix$EDdf[1:(which(ed_keys == key_ed) - 1), "Model"])
            Vg    <- ss[key_ed, "Variance"]
          }else{
            start <- sum(mix$EDdf[1:(which(mix$EDdf$Term == iGroup) - 1),"Model"]) # we don't add a one because we need the intercept
            Vg <- ss[iGroup,"Variance"]
          }
          
          C_sp  <- mix$C 
          
          if(calculateSE){
            if(verbose){message(paste("   Calculating standard errors for",iTrait, iGroup,"predictions"))}
              if(use_formula){
                # indices of this random term within C
                idx   <-  start:(start + nEffects - 1L)             # coefficient positions for this random term
                nC    <- nrow(mix$C)
                
                chunk <- 400L
                dvals <- numeric(length(idx))
                pos   <- 1L
                while (pos <= length(idx)) {
                  cols <- idx[pos:min(pos + chunk - 1L, length(idx))]
                  E    <- rhs_eye_spam(nC, cols)                      # nC x k (very skinny)
                  X    <- spam::solve(C_sp, E)                        # solves C %*% X = E
                  # each needed diagonal element is X[cols[j], j]
                  k    <- length(cols)
                  dvals[pos:(pos + k - 1L)] <- X[cbind(cols, seq_len(k))]
                  pos  <- pos + k
                }
                stdError    <- sqrt(pmax(dvals, 0))
              
              }else{
                #Covariates were used (Msub)
                stop <- start + nEffects - 1L
                idx_block <- start:stop
                
                startPev <- seq(1L, length(blup), by = 500L)
                endPev   <- c(startPev - 1L, length(blup)); endPev <- endPev[-1L]
                stdError <- vector("list", length(startPev))
                
                # Precompute C^{-1} columns for this block once
                nC  <- nrow(mix$C)
                Ebk <- rhs_eye_spam(nC, idx_block)
                # nC x k, k = nEffects
                
                Xbk <- spam::solve(C_sp, Ebk)               # nC x k  (C^{-1}[:, idx_block])
                Cinv_block <- Xbk[idx_block, , drop = FALSE]# k x k   (C^{-1}[idx_block, idx_block])
                
                for (b in seq_along(startPev)) {
                  use <- startPev[b]:endPev[b]
                  Muse <- Msub[[iGroup]][use, , drop = FALSE]
                  Mnum <- as.matrix(Muse)
                  
                  # diag( Muse %*% Cinv_block %*% t(Muse) ) = rowSums( (Muse %*% Cinv_block) * Muse )
                  Tmat <- as.matrix(Muse %*% Cinv_block)        # r x k
                  stdError[[b]] <- sqrt(pmax(rowSums(Tmat * Mnum), 0))
                }
                stdError <- unlist(stdError, use.names = FALSE)
              }
            reliability <- abs((Vg - (stdError^2)) /Vg) # reliability <- abs((Vg - Matrix::diag(pev))/Vg)
          }else{stdError <- reliability <- rep(NA,length(blup))}
          
          badRels <- which(reliability > 1); if(length(badRels) > 0){reliability[badRels] <- 0.9999}
          badRels2 <- which(reliability < 0); if(length(badRels2) > 0){reliability[badRels2] <- 0}
          
          envCol <- envsSub[[iGroup]]
          if (length(envCol) == 1L) envCol <- rep(envCol, length(blup))
          envCol <- unname(envCol)  # drop any names attribute on the short vector
          
          if (use_formula && grepl("_", iGroup)) {
            parts     <- strsplit(names(blup), ":", fixed = TRUE)
            term_vars <- unlist(randomTermSub[[iGroup]])
            
            # which slots are designation
            geno_vars <- c("designation", "designationA", "designationD",
                           "gid", "mother", "father")
            idx_env <- which(!(term_vars %in% geno_vars))
            
            if (length(idx_env) > 0L) {
              envCol <- vapply(parts, function(p) {
                if (length(p) < max(idx_env)) {
                  NA_character_
                } else {
                  paste(p[idx_env], collapse = ":")
                }
              }, character(1L))
            } else {
              envLevels <- envsSub[[iGroup]]
              envCol <- vapply(parts, function(p) {
                hit <- p[p %in% envLevels]
                if (length(hit) == 1L)      hit
                else if (length(hit) == 0L) NA_character_
                else                        hit[1L]
              }, character(1L))
            }
          }
          
          prov <- data.frame(designation=names(blup), predictedValue=blup, stdError=stdError, reliability=reliability,
                             trait=iTrait, effectType=iGroup , environment=envCol)
          
          # add fixed effects if present in the random term
          feToAdd <- intersect(unlist(randomTermSub[[iGroup]]), fixedEffects)
          if (length(feToAdd) > 0) {
            # split the compound label (designation column here encodes the crossed levels)
            varInppGroup <- strsplit(prov[,"designation"], ":", fixed = TRUE)
            
            for (iFe in feToAdd) {
              # which position of the crossed term corresponds to this fixed effect?
              pickVarInppGroup <- which(unlist(randomTermSub[[iGroup]]) == iFe)
              if (length(pickVarInppGroup) != 1L) next  # skip if ambiguous or absent
              
              # fixed-effect BLUEs table for this term
              provFe <- pp[[iFe]]
              if (is.null(provFe) || !NROW(provFe)) next
              
              # rownames = pure level name (strip iFe_ prefix)
              rn <- gsub(paste0("^", iFe, "_"), "", provFe[,"designation"])
              rownames(provFe) <- rn
              
              # extract the iFe level used in each crossed label of this random effect
              feUsed <- vapply(varInppGroup,
                               function(x) if (length(x) >= pickVarInppGroup) x[[pickVarInppGroup]] else NA_character_,
                               FUN.VALUE = character(1))
              
              # map to fixed BLUE; fallback to mu when not found
              mu0 <- provFe[feUsed, "predictedValue"]
              mu0[is.na(mu0)] <- mu
              
              prov[,"predictedValue"] <- prov[,"predictedValue"] + mu0
            }
          } else {
            prov[,"predictedValue"] <- prov[,"predictedValue"] + mu
          }
          
          # end of adding fixed effects
          sdP <- sd(prov[,"predictedValue"],na.rm=TRUE)
          cv <- (sd(prov[,"predictedValue"],na.rm=TRUE)/mean(prov[,"predictedValue"],na.rm=TRUE))*100
          
          ## PEV-corrected variance of this random effect (Fernández-González & Isidro y Sánchez, 2025, Eq. 4)
          var_PEV <- NA_real_
          
          if (isTRUE(calculateSE) && length(blup) > 0L) {
            ok <- which(!is.na(blup) & !is.na(stdError))
            if (length(ok) > 0L) {
              n_eff     <- length(ok)
              var_hat   <- sum(blup[ok]^2) / n_eff
              mean_pev  <- sum(stdError[ok]^2) / n_eff   # tr(PEV)/n using only diagonal
              var_PEV   <- var_hat + mean_pev            # σ^2 = var(â) + tr(PEV)/n
            }
          }
          
          
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
                                                  parameter=c( paste(c("mean","sd", "r2","Var","Var_PEVcorr"),iGroup,sep="_") ),
                                                  method=c("sum(x)/n","sd","(G-PEV)/G","REML","REML"),
                                                  value=c(mean(prov[,"predictedValue"], na.rm=TRUE), sdP, median(reliability), var(prov[,"predictedValue"], na.rm=TRUE), var_PEV),
                                                  stdError=c(NA,NA,sd(reliability, na.rm = TRUE)/sqrt(length(reliability)),NA,NA)
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
    names(match1) <- unlist(lapply(fixedTermSub,function(x){paste(x,collapse = ":")}))
    match2 <- unlist(lapply(randomTermSub,function(x){sum(as.numeric(x=="designation"))}))

    match3 <- c(match1,match2)
    useForPreds <- names(match3)[which(match3 > 0)]
    useForPreds <- gsub(":", "_", useForPreds)
    
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
  ## end of model fitting
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

  ##Adapt exports to SCA and GCA models
  if(is_SCA_GCA | is_GCA){
    predictionsBind[predictionsBind$effectType == "designation", "effectType"] = "designation_SCA"
    predictionsBind[predictionsBind$effectType == "mother", "effectType"] = "mother_GCA"
    predictionsBind[predictionsBind$effectType == "father", "effectType"] = "father_GCA"
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
