metASREML <- function(phenoDTfile = NULL,
                      analysisId = NULL,
                      analysisIdgeno = NULL,
                      #gsca = FALSE,
                      fixedTerm = NULL,
                      randomTerm = NULL,
                      covMod = NULL,
                      addG = NULL,
                      nFA = NULL,
                      envsToInclude = NULL,
                      trait = NULL,
                      traitFamily = NULL,
                      useWeights = TRUE,
                      calculateSE = TRUE,
                      heritLB = 0.15,
                      heritUB = 0.95,
                      meanLB = 0,
                      meanUB = Inf,
                      maxIters = 50,
                      verbose = TRUE)
{
  #save(phenoDTfile,analysisId,analysisIdgeno,fixedTerm,randomTerm,covMod, addG, nFA,envsToInclude, trait, traitFamily, useWeights,calculateSE, heritLB,  heritUB, meanLB, meanUB, maxIters,  verbose, file="NewAsr.RData")
  #library(asreml)
  '%!in%' <- function(x, y){! ('%in%'(x, y))}
  covMod <- lapply(covMod, gsub, pattern = "\\.", replacement = "")
  addG <- covMod
  ## THIS FUNCTION PERFORMS A MULT TRIAL ANALYSIS USING asreml
  mtaAnalysisId <- as.numeric(Sys.time())
  namesSeq <- function(x) {
    nCharX <- nchar(x)
    maxZeros <- max(nCharX)
    nZeros <- abs(nCharX - maxZeros)
    zeros <- apply(data.frame(nZeros), 1, function(x) {
      paste(rep("0", x), collapse = "")
    })
    res <- paste0(zeros, as.character(x))
    return(res)
  }
  ##########################################
  ##########################################
  ## CONTROLS FOR MISSPECIFICATION (6 lines)
  if (is.null(phenoDTfile)) {
    stop("Please provide the phenotype file", call. = FALSE)
  }
  if (is.null(analysisId)) {
    stop("Please provide the STA analysisId to be analyzed", call. = FALSE)
  }
  if (is.null(trait)) {
    stop("Please provide traits to be analyzed", call. = FALSE)
  } else{
    baseData <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId), ]
    if (length(intersect(trait, unique(baseData[, "trait"]))) == 0) {
      stop("The traits you have specified are not present in the analysisId provided.",
           call. = FALSE)
    }
  }
  if (is.null(traitFamily)) {
    traitFamily <- rep("asr_gaussian(link = 'identity', dispersion = 1.0)", length(trait))
  }
  if (length(traitFamily) != length(trait)) {
    stop("Trait distributions should have the same length than traits to be analyzed.",
         call. = FALSE)
  }
  traitsForExpCovariates <- unique(phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId), "trait"])
  ##########################################
  ##########################################
  ## EXTRACT POSSIBLE EXPLANATORY COVARIATES AND FORM KERNELS (30 lines)
  Weather <- cgiarPipeline::summaryWeather(phenoDTfile, wide = TRUE) # in form of covariates
  Weather <- apply(Weather, 2, sommer::imputev)
  colnames(Weather) <- gsub(" ", "", colnames(Weather))
  
  covkernel = c(
    "Relationship structure_Pedigree",
    "Relationship structure_GenoA",
    "Relationship structure_GenoD",
    "Relationship structure_GenoAD"
  )
  covars <- unique(unlist(addG))
  #We can assign random term to relationship structure
  asignarRandom <- function(covkernel, covMod, randomTerm, fill = "none") {
    # ============================================
    #  VALIDACIÓN DE conflictOS ENTRE ESTRUCTURAS
    # ============================================
    # Crear data.frame combinando estructuras y random terms
    df <- data.frame(
      cov = unlist(covMod),
      rnd = unlist(randomTerm),
      stringsAsFactors = FALSE
    )
    covModtmp<-covMod
    # Filtrar estructuras válidas (no "none")
    df_struct <- df[df$cov != "none", ]
    # Identificar estructuras repetidas
    duplicated_structs <- names(which(table(df_struct$cov) > 1))
    # Procesar cada estructura repetida
    for (s in duplicated_structs) {
      subset_df <- df_struct[df_struct$cov == s, ]
      terms <- unique(subset_df$rnd)
      # Caso permitido solo si aparecen mother y father
      allowed <- (s == "Relationship structure_GenoA") && all(c("mother", "father") %in% terms) && length(terms) <= 2
      if (!allowed) {
        message(
          "The structure '", s,
          "' appears in different random terms: ",
          paste(terms, collapse = ", "),
          ". is only available when you have 'mother' and 'father' in GenoA. The structure was deleted."
        )
        # Eliminar esta estructura en todas las entradas de covMod
        covModtmp <- lapply(covMod, function(x) x[x != s])
      }else{
        message("Structure '",s,"' kept renaming by associated randomTerm: ",paste(terms, collapse=", "))
        covModtmp <- lapply(seq_along(covMod), function(i) {
          x <- covMod[[i]]
          r <- randomTerm[[i]]
          # buscar coincidencias exactas donde x == s
          idx <- which(x == s)
          if (length(idx) == 0) return(x)
          # por cada posición encontrada, renombrarla con su randomTerm asociado
          for (k in seq_along(idx)) {
            x[idx[k]] <- paste0(s, "_", r[idx[k]])
          }
          return(x)
        })
      }
    }
    # ============================================
    #  ASIGNACIÓN ORIGINAL
    # ============================================
    match_idx <- sapply(covMod, function(x) match(x, covkernel))
    asignacion <- vector("list", length(covkernel))
    for (i in seq_along(covMod)) {
      rT <- which(!is.na(match_idx[[i]]))
      if (length(rT)!=0) {
        pos <- match_idx[[i]][rT]
        if(is.null(asignacion[[pos]])){
          asignacion[[pos]] <- randomTerm[[i]][rT]
        }else{
          asignacion[[pos]] <- c(asignacion[[pos]],randomTerm[[i]][rT])
        }
      }
    }
    if (!is.null(fill)) {
      asignacion[sapply(asignacion, is.null)] <- fill
    }
    return(list(asignacion,covModtmp))
  }
  
  randomTermForCovars <- asignarRandom(covkernel,covMod,randomTerm)
  covMod<-randomTermForCovars[[2]]
  covars<-unlist(covMod)
  randomTermForCovars<-randomTermForCovars[[1]]
  fixedTermForCovars <- setdiff(unique(unlist(fixedTerm)), c("environment", "designation"))
  G=D=N=Gad=NULL
  
  if (any(covkernel %in% covars)) {
    #New structure Geno info
    if (class(phenoDTfile$data$geno)[1] == "genlight") {
      qas <- which(names(phenoDTfile$data$geno_imp) == analysisIdgeno)
      Markers <- as.data.frame(phenoDTfile$data$geno_imp[qas])
    } else{#Old structure Geno info
      qas <- which(phenoDTfile$status$module == "qaGeno")
      qas <- qas[length(qas)]
      if (length(qas) > 0) {
        modificationsMarkers <- phenoDTfile$modifications$geno[which(phenoDTfile$modifications$geno$analysisId %in% qas), ]
        Markers <- cgiarBase::applyGenoModifications(M = Markers, modifications =
                                                       modificationsMarkers)
        if (length(which(is.na(Markers))) > 0) {
          Markers <- apply(Markers, 2, sommer::imputev)
        }
      } else{
        missing <- apply(Markers, 2, sommer::propMissing)
        Markers <- apply(Markers[, which(missing < 0.9)], 2, sommer::imputev)
      }
    }
    ploidyFactor <- max(Markers) / 2
    mydata <- phenoDTfile$predictions
    mydata <- mydata[which(mydata$analysisId %in% analysisId), ]
    if ("Relationship structure_GenoA_mother" %in% covars) {
      # additive model
      message(paste(" Marker A mother kernel is requested"))
      #males <- as.character(unique(mydata$designation))
      males <- as.character(unique(mydata[,c("mother")]))
      check1<-length(intersect(rownames(Markers), c(males)))
      if(check1<length(males)/2){
        stop("Individuals associated with marker matrix should have the same length. Please check both sources",call. = FALSE)
      }
      MarkersAm <- Markers[intersect(rownames(Markers), c(males)), ]
      G <- sommer::A.mat(as.matrix(MarkersAm - ploidyFactor))
      missing <- setdiff(c(males), rownames(G))
      A1m <- diag(mean(diag(G)),
                  nrow = length(missing),
                  ncol = length(missing))
      rownames(A1m) <- colnames(A1m) <- missing
      G <- sommer::adiag1(G, A1m)
      Gm <- G + diag(1e-5, ncol(G), ncol(G))
    } # additive model mother
    if ("Relationship structure_GenoA_father" %in% covars) {
      # additive model
      message(paste(" Marker A father kernel is requested"))
      #males <- as.character(unique(mydata$designation))
      males <- as.character(unique(mydata[,c("father")]))
      check1<-length(intersect(rownames(Markers), c(males)))
      if(check1<length(males)/2){
        stop("Individuals associated with marker matrix should have the same length. Please check both sources",call. = FALSE)
      }
      MarkersAf <- Markers[intersect(rownames(Markers), c(males)), ]
      G <- sommer::A.mat(as.matrix(MarkersAf - ploidyFactor))
      missing <- setdiff(c(males), rownames(G))
      A1m <- diag(mean(diag(G)),
                  nrow = length(missing),
                  ncol = length(missing))
      rownames(A1m) <- colnames(A1m) <- missing
      G <- sommer::adiag1(G, A1m)
      Gf <- G + diag(1e-5, ncol(G), ncol(G))
    } # additive model father
    if ("Relationship structure_GenoA" %in% covars) {
      # additive model
      message(paste(" Marker A kernel is requested"))
      #males <- as.character(unique(mydata$designation))
        males <- as.character(unique(mydata[,randomTermForCovars[[2]]]))
        check1<-length(intersect(rownames(Markers), c(males)))
        if(check1<length(males)/2){
          stop("Individuals associated with marker matrix should have the same length. Please check both sources",call. = FALSE)
        }
        MarkersA <- Markers[intersect(rownames(Markers), c(males)), ]
        G <- sommer::A.mat(as.matrix(MarkersA - ploidyFactor))
        missing <- setdiff(c(males), rownames(G))
        A1m <- diag(mean(diag(G)),
                  nrow = length(missing),
                  ncol = length(missing))
        rownames(A1m) <- colnames(A1m) <- missing
        G <- sommer::adiag1(G, A1m)
        G <- G + diag(1e-5, ncol(G), ncol(G))
    } # additive model
    if ("Relationship structure_GenoD" %in% covars) {
      #Dominance kernel
      if (ploidyFactor == 1) {
        message(paste(" Marker D kernel is requested"))
        males <- as.character(unique(mydata[,randomTermForCovars[[3]]]))
        check1<-length(intersect(rownames(Markers), c(males)))
        if(check1<length(males)/2){
          stop("Individuals associated with marker matrix should have the same length. Please check both sources",call. = FALSE)
        }
        #males <- as.character(unique(mydata$designation))
        MarkersD <- Markers[intersect(rownames(Markers), c(males)), ]
        D <- 1 - abs(MarkersD)
        f <- rowSums(D) / ncol(D) #inbreeding fixed eff
        names(f) <- rownames(MarkersD)
        D <- sommer::D.mat(as.matrix(D))
        missing <- setdiff(c(males), rownames(D))
        A1m <- diag(mean(diag(D)),
                    nrow = length(missing),
                    ncol = length(missing))
        rownames(A1m) <- colnames(A1m) <- missing
        D <- sommer::adiag1(D, A1m)
        Gd <- D + diag(1e-5, ncol(D), ncol(D))
      } else{
        #autopolyploid formula for digenic dominance (Batista et al. 2022)
        message(paste(" Marker D kernel is requested"))
        ploidy <- ploidyFactor * 2
        males <- as.character(unique(mydata[,randomTermForCovars[[3]]]))
        check1<-length(intersect(rownames(Markers), c(males)))
        if(check1<length(males)/2){
          stop("Individuals associated with marker matrix should have the same length. Please check both sources",call. = FALSE)
        }
        MarkersD <- Markers[intersect(rownames(Markers), c(males)), ]
        dom_matrix = as.matrix(MarkersD / ploidy)
        dom_matrix = 4 * dom_matrix - 4 * (dom_matrix * dom_matrix)
        f <- rowSums(dom_matrix) / ncol(dom_matrix) #inbreeding fixed eff
        names(f) <- rownames(MarkersD)
        MAF <- colMeans(MarkersD, na.rm = TRUE) / ploidy
        tMarkers <- t(MarkersD)
        C_mat <- matrix(choose(ploidy, 2),
                        nrow = nrow(tMarkers),
                        ncol = ncol(tMarkers))
        Ploidy_mat <- matrix(ploidy,
                             nrow = nrow(tMarkers),
                             ncol = ncol(tMarkers))
        Q <- (MAF^2 * C_mat) - (Ploidy_mat - 1) * MAF * tMarkers + 0.5 * tMarkers * (tMarkers -
                                                                                       1)
        D <- crossprod(Q)
        denomDom <- sum(C_mat[, 1] * MAF^2 * (1 - MAF)^2)
        D <- D / denomDom
        Gd <- D + diag(1e-5, ncol(D), ncol(D))
      }
    }#Dominance kernel
    if ("Relationship structure_GenoAD" %in% covars) {
      # additive + dominance model
      message(paste(" Marker AD kernel is requested"))
      males <- as.character(unique(mydata[,randomTermForCovars[[4]]]))
      check1<-length(intersect(rownames(Markers), c(males)))
      if(check1<length(males)/2){
        stop("Individuals associated with marker matrix should have the same length. Please check both sources",call. = FALSE)
      }
      MarkersAD <- Markers[intersect(rownames(Markers), c(males)), ]
      MarkersAD <- apply(MarkersAD + 1, 2, log)
      Gad <- sommer::A.mat(as.matrix(MarkersAD))
      Gad <- Gad + diag(1e-5, ncol(Gad), ncol(Gad))
    } # additive + dominance model
    ## PEDIGREE KERNEL
    if ("Relationship structure_Pedigree" %in% covars) {
      message(paste(" Pedigree kernel is requested"))
      Pedigree <- phenoDTfile$data$pedigree[, 2:4]
      mydataX <-  phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId), ]
      Pedigree <- Pedigree[!duplicated(Pedigree[, 1]), ]
      Pedigree <- Pedigree[which(Pedigree[, 1] %in% unique(mydataX$designation)), ]
      N <- asreml::ainverse(Pedigree)
    }
  }
  ## COMPLETE THE CLEANING PARAMETERS (7 lines)
  names(traitFamily) <- trait
  heritLB <- rep(heritLB, length(trait))
  heritLB <- heritLB[1:length(trait)]
  names(heritLB) <- trait
  heritUB <- rep(heritUB, length(trait))
  heritUB <- heritUB[1:length(trait)]
  names(heritUB) <- trait
  meanLB <- rep(meanLB, length(trait))
  meanLB <- meanLB[1:length(trait)]
  names(meanLB) <- trait
  meanUB <- rep(meanUB, length(trait))
  meanUB <- meanUB[1:length(trait)]
  names(meanUB) <- trait
  traitOrig <- trait 
  # LOAD THE DATASET AND EXTEND IT TO INCLUDE METADATA (16 lines)
  message("Loading the dataset and adding metadata.")
  mydata <- phenoDTfile$predictions #
  mydata <- mydata[which(mydata$analysisId %in% analysisId), ]
  if (nrow(mydata) < 2){
    stop(
      "Not enough data is available to perform a multi trial analysis. Please perform an STA before trying to do an MET.",
      call. = FALSE
    )}
  metaPheno <- phenoDTfile$metadata$pheno[which(
    phenoDTfile$metadata$pheno$parameter %in% c(
      "pipeline",
      "stage",
      "environment",
      "year",
      "season",
      "timepoint",
      "country",
      "location",
      "trial",
      "study",
      "management"
    )
  ), ]
  otherMetaCols <- unique(phenoDTfile$data$pheno[, metaPheno$value, drop =
                                                   FALSE])
  colnames(otherMetaCols) <- cgiarBase::replaceValues(
    Source = colnames(otherMetaCols),
    Search = metaPheno$value,
    Replace = metaPheno$parameter
  )
  otherMetaCols <- otherMetaCols[which(!duplicated(otherMetaCols[, "environment"])), , drop =
                                   FALSE] # we do this in case the users didn't define the environment properly
  mydata <- merge(mydata, otherMetaCols, by = "environment", all.x = TRUE)
  WeatherRow <- as.data.frame(Weather)
  WeatherRow$environment <- rownames(WeatherRow)
  mydata <- merge(mydata, WeatherRow, by = "environment", all.x = TRUE)
  ##########################################
  ##########################################
  # CHECK THE ENVS TO INCLUDE PER TRAIT (6 lines)
  if (is.null(envsToInclude)) {
    envsToInclude =  as.data.frame(do.call(rbind, list (with(
      mydata, table(environment, trait)
    ))))
    bad <- which(envsToInclude <= 1, arr.ind = TRUE)
    if (nrow(bad) > 0) {
      envsToInclude[bad] = 0
    }
    envsToInclude[which(envsToInclude > 1, arr.ind = TRUE)] = 1
  }
  allEnvironments <- rownames(envsToInclude)
  ##########################################
  ##########################################
  # BUILD THE DATASETS FOR MODEL FITTING (100 lines)
  message("Building trait datasets.")
  metrics <- phenoDTfile$metrics
  metrics <- metrics[which(metrics$analysisId %in% analysisId), ]
  myDataTraits <- fixedTermTrait <- randomTermTrait <- groupingTermTrait <- Mtrait <- envsTrait <- entryTypesTrait <- list()
  randomTermModel <- fixedTermModel <- list()
  x_option <- function(x, y, z) {
    switch(
      x,
      "none" = paste0("idv(", y, ")"),
      "Structure model_fa" = paste0("fa(", y, ",", z, ")"),
      "Structure model_diag" = paste0("diag(", y, ")"),
      "Structure model_us" = paste0("us(", y, ")"),
      "Relationship structure_GenoA" = paste0("vm(", y, ",source = G)"),
      "Relationship structure_GenoA_mother" = paste0("vm(", y, ",source = Gm)"),
      "Relationship structure_GenoA_father" = paste0("vm(", y, ",source = Gf)"),
      "Relationship structure_Pedigree" = paste0("vm(", y, ",source = N)"),
      "Relationship structure_GenoAD" = paste0("vm(", y, ",source = Gad)"),
      "Relationship structure_GenoD" = paste0("vm(", y, ",source = Gd)"),
      stop("Invalid `x` value")
    )
  }
  for (iTrait in trait) {
    # iTrait = trait[1]
    # filter for records available
    vt <- which(mydata[, "trait"] == iTrait)
    if (length(vt) > 0) {
      # we have data for the trait
      prov <- mydata[vt, ]
      # filter by the environments to include
      vte <- which(prov[, "environment"] %in% rownames(envsToInclude)[as.logical(envsToInclude[, iTrait])])
      prov <- prov[vte, ]
      # remove bad environment based on h2 and r2
      pipeline_metricsSub <- metrics[which(
        metrics$trait == iTrait &
          metrics$parameter %in% c(
            "plotH2",
            "H2",
            "meanR2",
            "r2",
            apply(expand.grid(
              c("plotH2", "H2", "meanR2", "r2"),
              c("designation", "mother", "father")
            ), 1, function(f) {
              paste(f, collapse = "_")
            })
          )
      ), ]
      goodFields <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value >= heritLB[iTrait]) &
                                                       (pipeline_metricsSub$value <= heritUB[iTrait])
      ), "environment"])
      prov <- prov[which(prov$environment %in% goodFields), ]
      # remove bad environment based on environment means
      pipeline_metricsSub <- metrics[which(
        metrics$trait == iTrait &
          metrics$parameter %in% c(
            "plotH2",
            "H2",
            "meanR2",
            "r2",
            apply(expand.grid(
              c("mean"), c("designation", "mother", "father")
            ), 1, function(f) {
              paste(f, collapse = "_")
            })
          )
      ), ]
      goodFieldsMean <- unique(pipeline_metricsSub[which((pipeline_metricsSub$value > meanLB[iTrait]) &
                                                           (pipeline_metricsSub$value < meanUB[iTrait])), "environment"])
      prov <- prov[which(prov$environment %in% goodFieldsMean), ]
      
      #Add inbreeding coefficient to prov
      if ("Relationship structure_GenoD" %in% covars &
          (!"f" %in% colnames(prov))) {
        prov$inbreeding <- f[match(prov$designation, names(f))]
      }
      
      if (nrow(prov) > 0) {
        # if after filters there's still data for this trait we can continue and save the data
        if (var(prov[, "predictedValue"], na.rm = TRUE) > 0) {
          # check that there is variance
          # make new formula for this specific trait if it passed all the filters
          
          #fixed formula per trait
          #fixedTermTrait[[iTrait]] <- paste("asreml::asreml(fixed = predictedValue ~ 1")
          fixedTermTrait[[iTrait]] <- paste("predictedValue ~ 1")
          if (length(fixedTerm) == 1 &
              "none" %in% unique(unlist(fixedTerm))) {
            fixedTermprov = "none"
          } else{
            fixedTermprov = unique(unlist(lapply(fixedTerm, paste, collapse = ":")))
            if (length(grep("none", fixedTermprov)) != 0) {
              fixedTermprov = fixedTermprov[-grep("none", fixedTermprov)]
            }
          }
          fixedTermModel[[iTrait]] = fixedTermprov
          if (length(fixedTermprov) != 0 | !is.null(fixedTermprov)) {
            for (iFixed in 1:length(fixedTermprov)) {
              # for each element in the list # iFixed=1
              fixedTermTrait[[iTrait]] <- paste(fixedTermTrait[[iTrait]], fixedTermprov[iFixed], sep =
                                                  " + ")
            }
          }
          #fixedTermTrait[[iTrait]] <- paste0(fixedTermTrait[[iTrait]], ",")
          
          # random formula per trait
          #randomTermTrait[[iTrait]] <- paste("random = ~")
          randomTermTrait[[iTrait]] <- "~"
          randomTermprov = list()
          for (w in 1:length(randomTerm)) {
            if (length(covMod[[w]]) == 1) {
              randomTermprov[[w]] = x_option(covMod[[w]], randomTerm[[w]], nFA[[w]])
            } else{
              tmprand = list()
                for(rt in 1:length(randomTerm[[w]])){
                  tmprand[[rt]]=x_option(covMod[[w]][rt], randomTerm[[w]][rt], nFA[[w]])
                }
              randomTermprov[[w]] = unlist(tmprand)
            }
          }
          randomTermprov = unique(unlist(lapply(
            randomTermprov, paste, collapse = ":"
          )))
          randomTermModel[[iTrait]] = randomTermprov
          if (length(randomTermprov) != 0 | !is.null(randomTermprov)) {
            randomTermTrait[[iTrait]] <- paste(randomTermTrait[[iTrait]], randomTermprov[1], sep =" ")
            if (length(randomTermprov)>1){
              for (iRand in 2:length(randomTermprov)) {
                # for each element in the list # iFixed=1
                randomTermTrait[[iTrait]] <- paste(randomTermTrait[[iTrait]], randomTermprov[iRand], sep =" + ")
              }
            }
          }
          #randomTermTrait[[iTrait]] <- paste0(randomTermTrait[[iTrait]], ",")
          
          myDataTraits[[iTrait]] <- prov # dataset for this trait
          # end of formula formation
        }
      }
    }
  }
  ## MODEL FITTING
  message("Fitting a model.")
  
  predictionsList <- list()
  
  for (iTrait in names(myDataTraits)) {
    # # iTrait = trait[1]  iTrait="value"
    message(paste("Analyzing trait", iTrait))
    
    mydataSub <- myDataTraits[[iTrait]] # extract dataset
    fixedTermSub <- fixedTermTrait[[iTrait]]
    randomTermSub <- randomTermTrait[[iTrait]]
    
    ## deregress if needed
    VarFull <- var(mydataSub[, "predictedValue"], na.rm = TRUE) # total variance
    effectTypeTrait <- phenoDTfile$modeling[which(
      phenoDTfile$modeling$analysisId == analysisId &
        phenoDTfile$modeling$trait == iTrait &
        phenoDTfile$modeling$parameter == "designationEffectType"
    ), "value"]
    if (names(sort(table(effectTypeTrait), decreasing = TRUE))[1] == "BLUP") {
      # if STA was BLUPs deregress
      mydataSub$predictedValue <- mydataSub$predictedValue / mydataSub$reliability
    }
    ## calculate weights
    mydataSub <- mydataSub[with(mydataSub, order(environment)), ] # sort by environments
    #mydataSub$w <- 1 / (mydataSub$stdError^2) # add weights column
    w <- 1 / (mydataSub$stdError^2) # add weights column
    # warnin messages in weights use
    message(" Using weights in the analysis. Residual variance will be fixed to 1. ")
    
    #Generate envIndex for finlay and wilkinson model
    if(grepl("envIndex",randomTermSub)){
      envIndex <- ave(mydataSub$predictedValue, mydataSub$environment, FUN = function(x) mean(x) - mean(mydataSub$predictedValue))
      mydataSub$envIndex <- envIndex
    }
    #trasform to factor variables
    labfactor <- c(
      "environment",
      "module",
      "analysisId",
      "pipeline",
      "trait",
      "gid",
      "designation",
      "mother",
      "father",
      "entryType",
      "effectType",
      "year",
      "location",
      "management",
      "trial"
    )
    tranfactfixed <- which(unique(unlist(fixedTerm)) %in% labfactor == T)    
    if (length(tranfactfixed) != 0) {
      mydataSub[, unique(unlist(fixedTerm))[tranfactfixed]] = lapply(mydataSub[unique(unlist(fixedTerm))[tranfactfixed]], as.factor)
    }
    tranfactrandom <- which(unique(unlist(randomTerm)) %in% labfactor == T)
    if (length(tranfactrandom) != 0) {
      mydataSub[, unique(unlist(randomTerm))[tranfactrandom]] = lapply(mydataSub[unique(unlist(randomTerm))[tranfactrandom]], as.factor)
    }
    
    asreml::asreml.options(workspace = 300e7,pworkspace = 200e7,trace = T,ai.sing = T)
    
    fixed_formula <- tryCatch(as.formula(fixedTermSub), error = function(e) return(NULL))
    random_formula <- tryCatch(as.formula(randomTermSub), error = function(e) return(NULL))
    
    if (is.null(fixed_formula) | is.null(random_formula)) {
      return("⚠️ Error in the formula: check the sintaxis.")
    }
    
    family_arg <- eval(parse(text = traitFamily[iTrait]))
    # Adjust model with asreml
    tryCatch({
      mix <<- asreml::asreml(
        fixed = fixed_formula,
        random = random_formula,
        data = mydataSub,
        na.action = na.method(x='include', y='include'),
        maxit = maxIters,
        weights = w,
        family = family_arg,	
        envir = .GlobalEnv
      )      
      #summary(mix)
    }, error = function(e) {
      paste("❌ Error to adjust model:", e$message)
    })
  
  

    update_until_converged <- function(model, max_updates = 10, verbose = TRUE) {
      count <- 0
      while (!model$converge && count < max_updates) {
        count <- count + 1
        if (verbose) message("🔄 Update #", count, "...")
        model <- tryCatch(asreml::update.asreml(model), error = function(e) {
          warning("⚠️ Error while update(): ", e$message)
          return(model)
        })
      }
      
      if (model$converge) {
        if (verbose) message("✅ Convergence in " , count, " updates")
        } else {
        warning("❌ No convergence in ", max_updates, " updates")
      }
      
      return(model)
    }

    if (!inherits(mix, "try-error")) {
      mix <<- update_until_converged(mix, max_updates = 10)
    }
    
    pp <- list()
    if (!inherits(mix, "try-error")) {
    # get variance components
    ss <- summary(mix)$varcomp
    Ve <- ss["units!R", "component"]
    mu <- coef(mix)[["fixed"]]["(Intercept)", ]
    Ci <- mix[["vcoeff"]][["fixed"]][1] * mix$sigma2
    if (length(mu) > 0) {
      pp[["(Intercept)"]] <- data.frame(
        designation = "(Intercept)",
        predictedValue = mu,
        stdError = sqrt(Ci),
        reliability = NA,
        trait = iTrait,
        effectType = "(Intercept)",
        entryType = "(Intercept)",
        environment = "(Intercept)"
      )
    }
    #fixed effect start
    fixTermObt = unlist(fixedTermModel[[iTrait]])
    if (length(grep("none", fixTermObt)) == 0) {
      for (iGroupFixed in fixTermObt) {
        # iGroupFixed=unlist(fixedTermModel[[iTrait]])[1]
        pick <- coef(mix)[["fixed"]][grepl(iGroupFixed, rownames(coef(mix)[["fixed"]]))]
        if (length(pick) > 1) {
          names(pick) <- rownames(coef(mix)[["fixed"]])[grepl(iGroupFixed, rownames(coef(mix)[["fixed"]]))]
          pick <- pick[which(pick != 0)]
        } else{
          if (pick == 0) {
            pick <- rep((-1) * mu, 2)
            names(pick) <- rep(rownames(coef(mix)[["fixed"]])[grepl(iGroupFixed, rownames(coef(mix)[["fixed"]]))], 2)
            pick <- pick[1]
          } else{
            pick <- rep(pick, 2)
            names(pick) <- rep(rownames(coef(mix)[["fixed"]])[grepl(iGroupFixed, rownames(coef(mix)[["fixed"]]))], 2)
            pick <- pick[1]
          }
        }
        blue <- pick + mu
        names(blue) <- names(pick)
        #blue[1] <- blue[1]-mu
        pev <- mix[["vcoeff"]][["fixed"]][grepl(iGroupFixed, rownames(coef(mix)[["fixed"]]))] * mix$sigma2
        if (length(pick) >= 1) {
          pev <- pev[which(pev != 0)]
        }
        if (is.matrix(pev)) {
          stdError <- (sqrt(Matrix::diag(pev)))
        } else{
          stdError <- sqrt(pev)
        }
        prov <- data.frame(
          designation = names(blue),
          predictedValue = blue,
          stdError = stdError,
          reliability = NA,
          trait = iTrait,
          effectType = iGroupFixed,
          environment = "(Intercept)"
        )
        for (iLabel in unique(unlist(fixedTerm))) {
          prov$designation <- gsub(paste0(iLabel, "_"), "", prov$designation)
        }
        # add additional entry type labels
        mydataSub[, "designationXXX"] <- apply(mydataSub[, unlist(strsplit(iGroupFixed, ":")), drop =
                                                           FALSE], 1, function(x) {
                                                             paste(x, collapse = ":")
                                                           })
        prov$entryType <- apply(data.frame(prov$designation), 1, function(x) {
          found <- which(mydataSub[, "designationXXX"] %in% x)
          if (length(found) > 0) {
            x2 <- paste(sort(unique(toupper(
              trimws(mydataSub[found, "entryType"])
            ))), collapse = "#")
          } else{
            x2 <- "unknown"
          }
          return(x2)
        })
        prov$entryType <- cgiarBase::replaceValues(prov$entryType,
                                                   Search = "",
                                                   Replace = "unknown")
        # save
        pp[[iGroupFixed]] <- prov
      }
    }#fixed effect finish
    #Start random effects
    x_option2 <- function(x, y, z) {
      switch(x,
             "none" = paste0(y),
             "Structure model_fa" = paste0("fa(", y, ",", z, ")"),
             "Structure model_diag" = paste0(y),
             "Structure model_us" = paste0(y),
             "Relationship structure_GenoA" = paste0("vm(", y, ", source = G)"),
             "Relationship structure_GenoA_mother" = paste0("vm(", y, ", source = Gm)"),
             "Relationship structure_GenoA_father" = paste0("vm(", y, ", source = Gf)"),
             "Relationship structure_Pedigree" = paste0("vm(", y, ", source = N)"),
             "Relationship structure_GenoAD" = paste0("vm(", y, ", source = Gad)"),
             "Relationship structure_GenoD" = paste0("vm(", y, ", source = Gd)"),
             stop("Invalid `x` value")
      )
    }
    #Random terms genetics
    subgroupGen = c()
    oriGen = c()
    relstrGen = c()
    only1 = which(lapply(randomTerm, length) == 1)
    if(length(only1)!=0){
    count = 1
    for (sb in 1:length(only1)) {
      if (randomTerm[[only1[sb]]] %in% c("designation", "mother", "father", "gid") == T) {
        subgroupGen[count] = x_option2(covMod[[only1[sb]]], randomTerm[[only1[sb]]], nFA[[only1[sb]]])
        oriGen[count] = randomTerm[[only1[sb]]]
        relstrGen[count] = covMod[[only1[sb]]]
        count = count + 1
      }
    }
    }
    #Random terms interaction
    subgroupInt = c()
    calcGenCorr=c()
    only2 = which(lapply(randomTerm, length) == 2)
    if(length(only2)!=0){
      for (sb2 in 1:length(only2)) {
        ed = all(randomTerm[[only2[sb2]]] %in% c("environment", "designation"))
        if (ed == T) {
          subgroupInt[sb2] = paste0(
            x_option2(covMod[[only2[sb2]]][1], randomTerm[[only2[sb2]]][1], nFA[[only2[sb2]]]),
            ":",
            x_option2(covMod[[only2[sb2]]][2], randomTerm[[only2[sb2]]][2], nFA[[only2[sb2]]])
          )
          if (covMod[[only2[sb2]]][1]=="Structure model_fa" | covMod[[only2[sb2]]][2]=="Structure model_fa"){ 
            calcGenCorr[sb2]=nFA[[only2[sb2]]]
          }else{
            calcGenCorr[sb2]=0
          }
        }
      }
    }
    
    calcGenCorr=calcGenCorr[calcGenCorr != 0]
    if (length(calcGenCorr)==1){
      env_levels=levels(mydataSub$environment)
      # Extract coefficients (loadings + specific variances)
      coefs <- summary(mix)$varcomp
      # Create matrix for loadings
      Lambda <- matrix(NA, nrow = length(env_levels), ncol = calcGenCorr,
                       dimnames = list(env_levels, paste0("FA", 1:calcGenCorr)))
      # Fill Lambda with loadings
      for (f in 1:calcGenCorr) {
        for (e in env_levels) {
          pat <- paste0(e,"!fa",f)
          val <- coefs[grep(pat, rownames(coefs)), "component"]
          if (length(val) == 1) {
            Lambda[e, f] <- val
          }
        }
      }
      # Extract specific variances (Ψ)
      psi <- numeric(length(env_levels))
      names(psi) <- env_levels
      for (e in env_levels) {
        pat <- paste0(e,"!var")
        val <- coefs[grep(pat, rownames(coefs)), "component"]
        if (length(val) == 1) {
          psi[e] <- val
        }
      }
      # Genetic covariance matrix: G = Lambda %*% t(Lambda) + diag(psi)
      G <- Lambda %*% t(Lambda) + diag(psi)
      rownames(G) <- colnames(G) <- env_levels
      # Genetic correlation matrix
      D <- diag(1 / sqrt(diag(G)))
      cor_mat <- D %*% G %*% D
      rownames(cor_mat) <- colnames(cor_mat) <- env_levels
      upper <- which(upper.tri(cor_mat), arr.ind = TRUE)
      prov <- data.frame(
        designation = rownames(cor_mat)[upper[,1]],
        predictedValue = cor_mat[upper],
        stdError = rep(NA, dim(upper)[1]),
        reliability = rep(NA, dim(upper)[1]),
        trait = iTrait,
        effectType = "GenCorrMat" ,
        environment = colnames(cor_mat)[upper[,2]],
        entryType = "unknown"
      )
      pp[["GenCorrMat"]] <- prov
    }
    
    #if (length(which(subgroupGen%in%c("mother","father")== T)>=2)) {
    if(all(unlist(lapply(c("mother", "father"), function(x) any(grepl(x, subgroupGen)))))){
      gca = summary(mix, coef = TRUE)$coef.random
      prov <- data.frame(
        designation = rownames(gca),
        predictedValue = gca[, 1],
        stdError = gca[, 2],
        reliability = rep(NA, dim(gca)[1]),
        trait = iTrait,
        effectType = "GCA/SCA" ,
        environment = "overall",
        entryType = "unknown"
      )
      prov$effectType[which(grepl("mother", prov$designation, ignore.case = TRUE))]="mother_GCA"
      prov$effectType[which(grepl("father", prov$designation, ignore.case = TRUE))]="father_GCA"
      if(length(which(subgroupGen%in%"designation"==T)!=0)){
        prov$effectType[which(grepl("designation", rownames(gca), ignore.case = TRUE))]="designation_SCA"
      }
      pp[["GCA/SCA"]] <- prov
    }
    
    groupingSub = c(subgroupGen, subgroupInt)
    if (length(groupingSub)!=0){
    for (iGroup in groupingSub) {
      #for( iGroup in names(groupingSub)){ # iGroup=groupingSub[2]
      
      blup = predict(mix, classify = iGroup)$pvals
      if (all(blup$status == "Aliased")) {
        statusmetrics = "Aliased estimation, problems with the model, please check!"
        stdError <- reliability <- rep(NA, dim(blup)[1])
        Vg <- NA
        lsdt <- NA
      } else{
        statusmetrics = mix$converge
        message(paste(" Calculating standar errors for",iTrait,iGroup,"predictions"))
        stdError <- blup$std.error # random effect was just one column
        if (iGroup %in% subgroupGen) {Vg <- var(blup$predicted.value) + ((stdError^2)/length(stdError))} else {Vg<-NA}
        reliability <- abs((Vg - (stdError^2)) / Vg) # reliability <- abs((Vg - Matrix::diag(pev))/Vg)
        lsdt <- qt(1 - 0.05 / 2, round(mix$nedf)) * predict(mix, classify =iGroup)$avsed
      }
      
      badRels <- which(reliability > 1)
      if (length(badRels) > 0) {
        reliability[badRels] <- 0.9999
      }
      badRels2 <- which(reliability < 0)
      if (length(badRels2) > 0) {
        reliability[badRels2] <- 0
      }
      
      if (iGroup %in% subgroupInt) {
        effTypeSub = "environment_designation"
        envTypeSub = blup$environment
        nameblup1 = iGroup
      } else{
        if (relstrGen[which(subgroupGen %in% iGroup == T)] == "Relationship structure_GenoA") {
          effTypeSub = paste0(oriGen[which(subgroupGen %in% iGroup == T)], "A")
        }
        if (relstrGen[which(subgroupGen %in% iGroup == T)] == "Relationship structure_GenoA_mother") {
          effTypeSub = paste0(oriGen[which(subgroupGen %in% iGroup == T)], "A")
        }
        if (relstrGen[which(subgroupGen %in% iGroup == T)] == "Relationship structure_GenoA_father") {
          effTypeSub = paste0(oriGen[which(subgroupGen %in% iGroup == T)], "A")
        }
        if (relstrGen[which(subgroupGen %in% iGroup == T)] == "Relationship structure_GenoD") {
          effTypeSub = paste0(oriGen[which(subgroupGen %in% iGroup == T)], "D")
        }
        if (relstrGen[which(subgroupGen %in% iGroup == T)] == "Relationship structure_GenoAD") {
          effTypeSub = paste0(oriGen[which(subgroupGen %in% iGroup == T)], "AD")
        }
        if (relstrGen[which(subgroupGen %in% iGroup == T)] == "Relationship structure_Pedigree") {
          effTypeSub = paste0(oriGen[which(subgroupGen %in% iGroup == T)], "Ped")
        }
        if (relstrGen[which(subgroupGen %in% iGroup == T)] == "none") {
          effTypeSub = paste0(oriGen[which(subgroupGen %in% iGroup == T)], "Idv")
        }
        envTypeSub = "(Intercept)"
        nameblup1 = grep(oriGen[which(subgroupGen %in% iGroup == T)], names(blup))
      }
      
      
      if (ncol(blup) > 4) {
        namesblup = c(nameblup1, names(blup)[3:5])
        blup = data.frame(cbind(paste0(blup[, 1], ":", blup[, 2]), blup[, 3], blup[, 4], blup[, 5]))
        names(blup) = namesblup
        blup$predicted.value = as.numeric(blup$predicted.value)
      }
      
      prov <- data.frame(
        designation = blup[, nameblup1],
        predictedValue = blup$predicted.value,
        stdError = stdError,
        reliability = reliability,
        trait = iTrait,
        effectType = effTypeSub ,
        environment = envTypeSub
      )
      
      # end of adding fixed effects
      sdP <- sd(prov[, "predictedValue"], na.rm = TRUE)
      cv <- (sd(prov[, "predictedValue"], na.rm = TRUE) / mean(prov[, "predictedValue"], na.rm = TRUE)) * 100      
      # add additional entry type labels
      mydataSub[, "designationXXX"] <- apply(mydataSub[, "designation", drop =
                                                         FALSE], 1, function(x) {
                                                           paste(x, collapse = ":")
                                                         })
      prov$entryType <- apply(data.frame(prov$designation), 1, function(x) {
        found <- which(mydataSub[, "designationXXX"] %in% x)
        if (length(found) > 0) {
          x2 <- paste(sort(unique(toupper(
            trimws(mydataSub[found, "entryType"])
          ))), collapse = "#")
          
        } else{
          x2 <- "unknown"
        }
        return(x2)
      })
      prov$entryType <- cgiarBase::replaceValues(prov$entryType,
                                                 Search = "",
                                                 Replace = "unknown")
      # save
      pp[[iGroup]] <- prov
      phenoDTfile$metrics <- rbind(
        phenoDTfile$metrics,
        data.frame(
          module = "mtaAsr",
          analysisId = mtaAnalysisId,
          trait = iTrait,
          environment = paste(unique(envTypeSub), collapse = "_"),
          parameter = c(paste(
            c("mean", "sd", "r2", "Var", "CV%", "LSD95%"),
            effTypeSub,
            sep = "_"
          )),
          method = c(
            "sum(x)/n",
            "sd",
            "(G-PEV)/G",
            "REML",
            "(sd/mean)*100",
            "t*avsed"
          ),
          value = c(
            mean(prov[, "predictedValue"], na.rm = TRUE),
            sdP,
            median(reliability),
            var(prov[, "predictedValue"], na.rm = TRUE),
            cv,
            lsdt
          ),
          stdError = c(
            NA,
            NA,
            sd(reliability, na.rm = TRUE) / sqrt(length(reliability)),
            NA,
            NA,
            NA
          )
        )
      )
    }
    }else{
      statusmetrics = mix$converge
      }
    kernels = paste(unlist(addG)[which(
      unlist(addG) %in% c(
        "Relationship structure_GenoA",
        "Relationship structure_GenoD",
        "Relationship structure_Pedigree",
        "Relationship structure_GenoAD"
      ) == T
    )], collapse = ",")
    if (length(kernels) == 0) { kernels = "none" }
    ## save the modeling used
    currentModeling <- data.frame(
      module = "MtaAsr",
      analysisId = mtaAnalysisId,
      trait = iTrait,
      environment = c(rep("across", 4), "designation"),
      parameter = c(
        "fixedFormula",
        "randomFormula",
        "family",
        "convergence",
        rep("kernels", 1)
      ),
      value = c(
        fixedTermSub,
        randomTermSub,
        traitFamily[iTrait],
        statusmetrics,
        kernels
      )
    )
    phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[, colnames(phenoDTfile$modeling)])
    ## save the environments used goodFields
    currentModeling <- data.frame(
      module = "MtaAsr",
      analysisId = mtaAnalysisId,
      trait = iTrait,
      environment = allEnvironments,
      parameter = "includedInMta",
      value = ifelse(
        allEnvironments %in% unique(mydataSub$environment),
        TRUE,
        FALSE
      )
    )
    phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[, colnames(phenoDTfile$modeling)])
    
    }else{
      # if model failed
      message(paste("Mixed model failed for trait",iTrait,". Aggregating and assuming h2 = 0 \n"))
      
      Ve <- NA
      means <- aggregate(predictedValue ~ designation,
                         FUN = mean,
                         data = mydataSub)
      means$environment <- "(Intercept)"
      means$stdError <- sd(means$predictedValue)
      means$reliability <- 1e-6
      means$trait <- iTrait
      means$effectType <- "designation"
      means$entryType <- "unknown"
      sdP <- sd(means$predictedValue, na.rm = TRUE)
      cv <- (sd(means$predictedValue, na.rm = TRUE) / mean(means$predictedValue, na.rm =TRUE)) * 100
      ## save metrics
      phenoDTfile$metrics <- rbind(
        phenoDTfile$metrics,
        data.frame(
          module = "mtaAsr",
          analysisId = mtaAnalysisId,
          trait = iTrait,
          environment = "across",
          parameter = c("mean", "sd", "r2", "Var", "CV%", "LSD95%"),
          method = c(
            "sum(x)/n",
            "sd",
            "(G-PEV)/G",
            "REML",
            "(sd/mean)*100",
            "t*avsed"
          ),
          value = c(
            mean(means$predictedValue, na.rm = TRUE),
            sdP,
            NA,
            NA,
            NA,
            NA
          ),
          stdError = NA
        )
      )
      currentModeling <- data.frame(
        module = "mtaAsr",
        analysisId = mtaAnalysisId,
        trait = iTrait,
        environment = "across",
        parameter = c(
          "fixedFormula",
          "randomFormula",
          "family",
          "designationEffectType"
        ),
        value = c("None", "None", "None", "mean")
      )
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[, colnames(phenoDTfile$modeling)])
      pp[["designation"]] <- means
    }
    predictionsTrait <- do.call(rbind, pp)
    ## add across env estimate for DESIGNATION effect type fitted
    match1 <- unlist(lapply(unlist(fixedTerm), function(x) {sum(as.numeric(grep("designation", x)))}))
    names(match1) <- unlist(lapply(unlist(fixedTerm), function(x) {paste(x, collapse = "_")}))
    match2 <- unlist(lapply(randomTerm, function(x) {sum(as.numeric(grep("designation", x)))}))
    names(match2) <- unlist(lapply(randomTerm, function(x) {
      paste(x, collapse = "_")
    }))
    match3 <- c(match1, match2)
    useForPreds <- unique(names(match3)[which(match3 > 0)])
    doublematch <- table(predictionsTrait$effectType,
                         predictionsTrait$environment)
    useForPreds <- rownames(doublematch)[unique(unlist(apply(data.frame(useForPreds), 1, function(x) {
      grep(x, rownames(doublematch))
    })))]
    interceptCheck <- sum(apply(data.frame(useForPreds), 1, function(x) {
      sum(as.numeric("(Intercept)" %in% colnames(doublematch)[which(doublematch[x, ] >
                                                                      0)]))
    }))
    if (length(useForPreds) > 0 & interceptCheck == 0) {
      # only if there was designation and no main effect exist then we aggregate
      provx <- predictionsTrait
      provx <- provx[which(provx$effectType %in% useForPreds), ]
      provx$designation <- apply(provx[, c("environment", "designation")], 1, function(x) {
        gsub(paste0(x[1], ":"), "", x[2])
      })
      provx <- aggregate(
        cbind(predictedValue, stdError, reliability) ~ designation + trait,
        FUN = mean,
        data = provx
      )
      provx$environment <- "(Intercept)"
      provx$effectType <- "designation"
      provx$entryType <- apply(data.frame(provx$designation), 1, function(x) {
        found <- which(mydataSub[, "designationXXX"] %in% x)
        if (length(found) > 0) {
          x2 <- paste(sort(unique(toupper(
            trimws(mydataSub[found, "entryType"])
          ))), collapse = "#")
          
        } else{
          x2 <- "unknown"
        }
        return(x2)
      })
      provx$entryType <- cgiarBase::replaceValues(provx$entryType, Search = "", Replace = "unknown")
      predictionsTrait <- rbind(predictionsTrait, provx[, colnames(predictionsTrait)])
    }
    ## add across env estimate for GID effect type fitted
    match1 <- unlist(lapply(unlist(fixedTerm), function(x) {sum(as.numeric(x == "gid"))}))
    names(match1) <- unlist(lapply(unlist(fixedTerm), function(x) {paste(x, collapse = "_")}))
    match2 <- unlist(lapply(randomTerm, function(x) {sum(as.numeric(x == "gid"))}))
    names(match2) <- unlist(lapply(randomTerm, function(x) {paste(x, collapse = "_")}))
    match3 <- c(match1, match2)
    useForPreds <- names(match3)[which(match3 > 0)]
    doublematch <- table(predictionsTrait$effectType,
                         predictionsTrait$environment)
    useForPreds <- rownames(doublematch)[unique(unlist(apply(data.frame(useForPreds), 1, function(x) {
      grep(x, rownames(doublematch))
    })))]
    interceptCheck <- sum(apply(data.frame(useForPreds), 1, function(x) {
      sum(as.numeric("(Intercept)" %in% colnames(doublematch)[which(doublematch[x, ] >
                                                                      0)]))
    }))
    if (length(useForPreds) > 0 & interceptCheck == 0) {
      # only if there was designation and no main effect exist then we aggregate
      provx <- predictionsTrait
      provx <- provx[which(provx$effectType %in% useForPreds), ]
      provx$designation <- apply(provx[, c("environment", "designation")], 1, function(x) {
        gsub(paste0(x[1], ":"), "", x[2])
      })
      provx <- aggregate(
        cbind(predictedValue, stdError, reliability) ~ designation + trait,
        FUN = mean,
        data = provx
      )
      provx$environment <- "(Intercept)"
      provx$effectType <- "gid"
      provx$entryType <- apply(data.frame(provx$designation), 1, function(x) {
        found <- which(mydataSub[, "designationXXX"] %in% x)
        if (length(found) > 0) {
          x2 <- paste(sort(unique(toupper(
            trimws(mydataSub[found, "entryType"])
          ))), collapse = "#")
          
        } else{
          x2 <- "unknown"
        }
        return(x2)
      })
      provx$entryType <- cgiarBase::replaceValues(provx$entryType, Search = "", Replace = "unknown")
      predictionsTrait <- rbind(predictionsTrait, provx[, colnames(predictionsTrait)])
    }
    phenoDTfile$metrics <- rbind(
      phenoDTfile$metrics,
      data.frame(
        module = "mtaAsr",
        analysisId = mtaAnalysisId,
        trait = iTrait,
        environment = "across",
        parameter = c("Var_residual", "nEnv", "nEntries"),
        method = c("REML", "n", "n"),
        value = c(Ve, length(goodFields), length(unique(
          mydataSub$designation
        ))),
        stdError = c(NA, NA, NA)
      )
    )
    predictionsList[[iTrait]] <- predictionsTrait
  }
  
  if (length(predictionsList) == 0) {
    stop(
      "There was no predictions to work with. Please look at your H2 boundaries. You may be discarding all envs.",
      call. = FALSE
    )
  }
  
  predictionsBind <- do.call(rbind, predictionsList)
  predictionsBind$analysisId <- mtaAnalysisId
  ## add timePoint of origin, stage and designation code
  entries <- unique(mydata[, "designation"])
  baseOrigin <- do.call(rbind, apply(data.frame(entries), 1, function(x) {
    out1 <- (sort(mydata[which(mydata$designation %in% x), "gid"], decreasing = FALSE))[1]
    out2 <- (sort(mydata[which(mydata$designation %in% x), "mother"], decreasing = FALSE))[1]
    out3 <- (sort(mydata[which(mydata$designation %in% x), "father"], decreasing = FALSE))[1]
    out4 <- paste(unique(sort(mydata[which(mydata$designation %in% x), "pipeline"], decreasing = FALSE)), collapse =
                    ", ")
    y <- data.frame(
      designation = x,
      gid = out1,
      mother = out2,
      father = out3,
      pipeline = out4
    )
    return(y)
  }))
  predictionsBind <- merge(predictionsBind,
                           baseOrigin,
                           by = "designation",
                           all.x = TRUE)
  predictionsBind$module <- "mtaAsr"
  rownames(predictionsBind) <- NULL
  
  if (!is.null(phenoDTfile$predictions)) {
    if ("effectType" %!in% colnames(phenoDTfile$predictions)) {
      phenoDTfile$predictions$effectType <- NA
    }
  }
  ##Adapt exports to GPCP model
  desigName = c(
    "designationA",
    "designationD",
    "designationAD",
    "designationPed",
    "designationIdv"
  )
  if (any(desigName %in% predictionsBind$effectType)) {
    more2 = desigName[which(desigName %in% predictionsBind$effectType == T)]
    avgDes=list()
    Bytraits=unique(predictionsBind$trait)
    for(tt in 1:length(Bytraits)){
      predictionsBindtmp=subset(predictionsBind,predictionsBind$trait==Bytraits[tt])
      avgDes[[tt]] <- predictionsBindtmp[which(predictionsBindtmp$effectType == more2[1]), ]
      if (length(more2) > 1) {
      # Filtrar
        df_sub <- predictionsBindtmp[predictionsBindtmp$effectType %in% more2, ]
        avgDes[[tt]]$predictedValue <- aggregate(
          predictedValue ~ designation,
          data = df_sub,
          FUN = function(x)
          mean(x, na.rm = TRUE)
        )[, 2]
        avgDes[[tt]]$stdError <- aggregate(
          stdError ~ designation,
          data = df_sub,
          FUN = function(x)
            mean(x, na.rm = TRUE)
        )[, 2]
        avgDes[[tt]]$reliability <- aggregate(
          reliability ~ designation,
          data = df_sub,
          FUN = function(x)
            mean(x, na.rm = TRUE)
        )[, 2]
        avgDes[[tt]]$effectType <- "designation"
      } else{
        avgDes[[tt]]$effectType <- "designation"
      }
    }  
    avgDes=do.call(rbind,avgDes)
    # Add the averaged designation rows
    predictionsBind <- rbind(predictionsBind, avgDes)
  }
    
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions, predictionsBind[, colnames(phenoDTfile$predictions)])
  newStatus <- data.frame(module = "mtaAsr",
                          analysisId = mtaAnalysisId,
                          analysisIdName = NA)
  phenoDTfile$status <- rbind(phenoDTfile$status, newStatus[, colnames(phenoDTfile$status)])
  ## add which data was used as input
  modeling <- data.frame(
    module = "mtaAsr",
    analysisId = mtaAnalysisId,
    trait = c("inputObject"),
    environment = "general",
    parameter = c("analysisId"),
    value = c(analysisId)
  )
  
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  message("Wrapping the results.")
  
  #save(phenoDTfile, file = "asrResult1.RData")
  
  return(phenoDTfile)
}



