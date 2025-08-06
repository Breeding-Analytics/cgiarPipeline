### 1. RR-BLUP / GBLUP - based GWAS model
gwas <- function (
    analysisId = NULL,
    analysisIdName = "",
    analysisIdForGenoModifications = NULL,
    phenoDTfile = NULL,
    trait = NULL,
    field = NULL,
    modelGwas = "gBLUP",
    maxIters = 50, 
    tolParInv = 1e-4,
    verbose = FALSE) {
  
  gwasAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the phenotype file", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the analysisId to be analyzed", call. = FALSE)}
  if(is.null(analysisIdForGenoModifications)){stop("Please provide the analysisId to be analyzed", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(phenoDTfile$data$geno)){ # modelType %in% c("gBLUP","rrBLUP") &
    stop("Please include marker information in your data structure to fit this model type.", call. = FALSE)
  }
  
  # if(modelGwas %in% c("gBLUP","rrBLUP") & !is.null(phenoDTfile$data$geno) & !is.null(phenoDTfile$metadata$geno)){
  if(!is.null(phenoDTfile$data$geno) & !is.null(phenoDTfile$metadata$geno)){
    Markers <- phenoDTfile$data$geno
    if(is.null(analysisIdForGenoModifications)){ # user didn't provide a modifications id
      if(length(which(adegenet::glNA(Markers) > 0)) > 0){stop("Markers have missing data and you have not provided a modifications table to impute the genotype data. Please go to the 'Markers QA/QC' module prior to run a gBLUP or rrBLUP model.", call. = FALSE)}
    }else{ # user provided a modifications Id
      modificationsMarkers <- phenoDTfile$modifications$geno
      theresMatch <- which(modificationsMarkers$analysisId %in% analysisIdForGenoModifications)
      if(length(theresMatch) > 0){ # there's a modification file after matching the Id
        modificationsMarkers <- modificationsMarkers[theresMatch,]
        #Markers <- phenoDTfile$data$geno_imp[[ as.character(as.integer(analysisIdForGenoModifications)) ]]
        Markers <- phenoDTfile$data$geno_imp[[grep(as.character(as.integer(analysisIdForGenoModifications)),names(phenoDTfile$data$geno_imp))]]  
      }else{ # there's no match of the modification file
        if(length(which(adegenet::glNA(Markers) > 0)) > 0){stop("Markers have missing data and your Id didn't have a match in the modifications table to impute the genotype data.", call. = FALSE)}
      }
    }
  }
  
  entryType <- unique(phenoDTfile$predictions[phenoDTfile$predictions$module %in% c("sta","mtaLmms"),"entryType"])
  entryType <- ifelse(length(entryType) > 1, sort(entryType)[1], entryType)
  
  for(iTrait in trait) {
    
    for(iField in field) {
      
      ### 2. Preparing marker data
      marker_map <- data.frame(SNP = Markers@loc.names,
                               CHR = Markers@chromosome,
                               BP = Markers@position)
      
      ### 1. Preparing phenotypic data
      "%>%" <- dplyr::`%>%`
      phenoData <- phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId == analysisId),]
      phenoData <- phenoData[which(phenoData$environment == iField),]
      phenoData <- phenoData %>%
        dplyr::select(ID = designation, trait, value = predictedValue) %>%
        tidyr::spread(key = trait, value = value) %>%
        dplyr::arrange(ID) %>%
        dplyr::filter(!is.na(!!dplyr::sym(iTrait))) %>%
        dplyr::select(c("ID", all_of(iTrait)))
      
      ### 3. Filtering for genotypes with phenotype and genotype data
      common_ids <- intersect(Markers@ind.names, phenoData$ID)
      markers <- as.matrix(Markers[which(Markers@ind.names %in% common_ids),])
      phenoData <- phenoData[phenoData$ID %in% common_ids, ]
      
      ### 4. Calculate GRM
      grm <- sommer::A.mat(markers, min.MAF = 0, return.imputed = TRUE) # monomorphic markers are removed!
      gmat_imputed <- grm$X

      ### 5. Remove filtered markers from marker map (at least monomorphic marker will be removed)
      ids_to_keep <- colnames(markers)
      marker_map <- marker_map[match(ids_to_keep, marker_map$SNP), , drop = FALSE]
      
      #
      meta <- marker_map %>%
        dplyr::group_by(CHR) %>%
        dplyr::summarise(N = length(CHR))
      
      CORR <- list()
      CORR[[1]] <- cor(markers[,1:meta$N[meta$CHR==1]])
      for(i in 2:nrow(meta)){
        dat = markers[,(meta$N[meta$CHR==(i-1)]+1):(meta$N[meta$CHR==(i-1)]+meta$N[meta$CHR==i])]
        CORR[[i]] = cor(dat)
      }
      z = list()
      for(i in 1:nrow(meta)){
        for(j in 1:nrow(CORR[[i]])){
          xx = as.vector(which(is.na(CORR[[i]][j,])))
        }
        if(length(xx)>0) {z[[i]] = CORR[[i]][-xx,-xx]} else {z[[i]] = CORR[[i]]}
      }
      corr = as.vector(NULL)
      for(i in 1:nrow(meta)){
        corr[i] = poolr::meff(poolr::mvnconv(z[[i]], target = "p", cov2cor = TRUE), method = "liji") # z[[i]]
      }
      Meff = sum(corr)
      p_threshold = (1-(1-0.05))/Meff # p-value to use as threshold for significance
      threshold = -log10(p_threshold)
      #
      
      common_rows <- intersect(rownames(markers), rownames(gmat_imputed))
      common_cols <- intersect(colnames(markers), colnames(gmat_imputed))
      gmat_imputed <- gmat_imputed[rownames(gmat_imputed) %in% common_rows, ]
      gmat_imputed <- gmat_imputed[, colnames(gmat_imputed) %in% common_cols]
      
      ### 6. Run GWAS using RR-BLUP - based GWAS.
      
      # Q: Any possible scenario where the number of fixed effects is > 1 in the BA analysis pipeline?
      k <- 1 # to be used for degrees of freedom (number of levels in fixed effects)
      
      #if (nrow(markers) >= ncol(markers)) {
      if (modelGwas %in% "rrBLUP") {
        mixRRBLUP <- sommer::mmer(as.formula(paste(iTrait,"~1")),
                                  random = ~sommer::vsr(gmat_imputed),
                                  rcov = ~units,
                                  nIters = maxIters,
                                  tolParInv = tolParInv,
                                  # emWeight = c(rep(1,4), seq(.9,.1,-.1), rep(.1,100)),  Q: @Eduardo: include or not?
                                  verbose = verbose,
                                  dateWarning = FALSE,                                  # @Eduardo: still give me a date warning.
                                  data = phenoData)
        
        a <- mixRRBLUP$U$`u:gmat_imputed`[[1]] # marker effects
        se_a <- sqrt(diag(kronecker(diag(ncol(gmat_imputed)), mixRRBLUP$sigma$`u:gmat_imputed`) -
                            mixRRBLUP$PevU$`u:gmat_imputed`[[1]]))
        
        t_stat <- a/se_a
        marker_map$P <- dt(t_stat, df = nrow(phenoData) - k - 1)
        
        marker_map <- marker_map[marker_map$P != 0, ]
        
      } else if(modelGwas %in% "gBLUP") {
        
        MMT <- tcrossprod(markers)
        MMTinv <- solve(MMT + diag(1e-6, ncol(MMT), ncol(MMT)))
        MTMMTinv <- t(markers) %*% MMTinv
        
        mixGBLUP <- sommer::mmer(as.formula(paste(iTrait,"~1")),
                                 random = ~sommer::vsr(ID, Gu = MMT),
                                 rcov = ~units,
                                 nIters = maxIters,
                                 tolParInv = tolParInv,
                                 # emWeight = c(rep(1,4), seq(.9,.1,-.1), rep(.1,100)),
                                 verbose = verbose,
                                 dateWarning = FALSE,
                                 data = phenoData)
        
        u <- matrix(mixGBLUP$U$`u:ID`[[1]], ncol = 1)
        rownames(u) <- names(mixGBLUP$U$`u:ID`[[1]])
        u <- as.matrix(u[match(colnames(MTMMTinv), rownames(u)), ])
        
        pev <- mixGBLUP$PevU$`u:ID`[[1]]
        pev <- pev[, match(colnames(MMT), colnames(pev))]
        pev <- pev[match(rownames(MMT), rownames(pev)), ]
        
        a_from_g <- MTMMTinv %*% u
        var_g <- kronecker(MMT, mixGBLUP$sigma$`u:ID`) - pev
        var_a_from_g <- t(markers) %*% MMTinv %*% (var_g) %*% t(MMTinv) %*% markers
        
        diag_var_a_from_g <- diag(var_a_from_g)
        keep_var_a_g <- (diag_var_a_from_g > 0)
        
        diag_var_a_from_g <- diag_var_a_from_g[keep_var_a_g]
        se_a_from_g <- sqrt(diag_var_a_from_g)
        a_from_g <- a_from_g[keep_var_a_g]
        marker_map <- marker_map[keep_var_a_g, ]
        
        t_stat_from_g <- a_from_g / se_a_from_g # t-statistic
        marker_map$P <- dt(t_stat_from_g, df = nrow(phenoData) - k - 1)
        
      }
      
      marker_map <- marker_map[!is.na(marker_map$CHR),]
      don <- marker_map %>%
        dplyr::group_by(CHR) %>%
        dplyr::summarise(chr_len = max(BP, na.rm = TRUE)) %>%
        dplyr::mutate(tot=cumsum(chr_len) - chr_len) %>%
        dplyr::left_join(marker_map, ., by=c("CHR" = "CHR")) %>%
        dplyr::arrange(CHR, BP) %>%
        dplyr::mutate(BPcum=BP+tot)
      
      don$BPcum <- as.numeric(don$BPcum)
      # don$text <- paste("SNP: ", don$SNP, "\nPosition: ", don$BP, "\nChromosome: ", don$CHR, "\n-log10(P) :", -log10(don$P) %>% round(2), sep="")
      
      don_sorted <- don[order(don$P), c(1:4)]
      don_sorted$"-log10(P)" <- -log10(don_sorted$P) # "-log10(P)"
      
      # output <- list(don, don_sorted) # phenoDTfile
      # pipeline <- unique(phenoDTfile$predictions$pipeline)
      don_sorted <- data.frame(module = "gwas", analysisId = gwasAnalysisId, pipeline = NA, trait = iTrait, gid = NA, designation = don_sorted$SNP,
                               mother = NA, father = NA, entryType = entryType,
                               environment = iField, predictedValue = don_sorted$`-log10(P)`,
                               stdError = don_sorted$CHR, reliability = don_sorted$BP,
                               effectType = NA)
      phenoDTfile$predictions <- rbind(phenoDTfile$predictions,don_sorted)
      #
      modeling <- data.frame(module="gwas", analysisId=gwasAnalysisId,trait=iTrait, environment=iField,
                             parameter="mta", value=TRUE)
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,modeling)
      #
      modeling <- data.frame(module="gwas", analysisId=gwasAnalysisId,trait=iTrait, environment=iField,
                             parameter=c("fixedFormula","randomFormula"),
                             value=c("predictedValue ~ 1", ifelse(modelGwas %in% "rrBLUP", "~ vsr(gmat_imputed)", "~ vsr(ID, Gu = MMT)")))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,modeling)
      #
      modeling <- data.frame(module="gwas",  analysisId=gwasAnalysisId, trait=iTrait, environment=iField,
                             parameter= c("Model"), value= ifelse(modelGwas %in% "rrBLUP", "rrBLUP", "gBLUP") )
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling)
      
      # loading the metrics
      metrics <- data.frame(module="gwas", analysisId=gwasAnalysisId,trait=iTrait, environment=iField,
                            parameter="threshold", method="Li and Ji (2005)", value=threshold, stdError=0)
      phenoDTfile$metrics <- rbind(phenoDTfile$metrics, metrics)
    }
  }
  phenoDTfile$status <- rbind(phenoDTfile$status, data.frame(module="gwas", analysisId=gwasAnalysisId, analysisIdName = analysisIdName))
  ## add which data was used as input
  modeling <- data.frame(module="gwas",  analysisId=gwasAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling)
  #
  #colnames(phenoDTfile$gwas)[length(colnames(phenoDTfile$gwas))] <- "-log10(P)"
  return(phenoDTfile)
}
