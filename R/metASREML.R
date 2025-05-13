metASREML <- function(
    phenoDTfile= NULL, analysisId=NULL,analysisIdgeno=NULL,gsca=FALSE,
    fixedTerm= NULL,  randomTerm=NULL, covMod=NULL, addG=NULL, nFA=NULL,
    envsToInclude=NULL, trait= NULL, traitFamily=NULL, useWeights=TRUE,
    calculateSE=TRUE, heritLB= 0.15,  heritUB= 0.95,
    meanLB=0, meanUB=Inf, maxIters=50,  verbose=TRUE
){
#save(phenoDTfile,analysisId,analysisIdgeno,gsca,fixedTerm,randomTerm,covMod, addG, nFA,envsToInclude, trait, traitFamily, useWeights,calculateSE, heritLB,  heritUB, meanLB, meanUB, maxIters,  verbose, file="NewAsr.RData")
  covMod<-lapply(covMod,gsub,pattern="\\.",replacement="")
  addG<-covMod
  ## THIS FUNCTION PERFORMS A MULT TRIAL ANALYSIS USING asreml
  mtaAnalysisId <- as.numeric(Sys.time())
  #fixedTerm<-as.list(fixedTerm)
  #randomTerm<-as.list(randomTerm)
  namesSeq <- function(x){
    nCharX <- nchar(x)
    maxZeros <- max(nCharX)
    nZeros <- abs(nCharX - maxZeros)
    zeros <- apply(data.frame(nZeros),1,function(x){paste(rep("0",x), collapse = "")})
    res <- paste0(zeros, as.character(x))
    return(res)
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
  if(is.null(traitFamily)){traitFamily <- rep("asr_gaussian(link = 'identity')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  #if(length(fixedTerm)==1 & unique(unlist(fixedTerm)) == "none"){fixedTerm <- "1"}
  traitsForExpCovariates <- unique(phenoDTfile$predictions[which(phenoDTfile$predictions$analysisId %in% analysisId),"trait"])
  ##########################################
  ##########################################
  ## EXTRACT POSSIBLE EXPLANATORY COVARIATES AND FORM KERNELS (30 lines)
  Weather <- cgiarPipeline::summaryWeather(phenoDTfile, wide=TRUE) # in form of covariates
  Weather <- apply(Weather,2,sommer::imputev)
  colnames(Weather) <- gsub(" ","",colnames(Weather))
  
    
  covars <- unique(unlist(addG))
  randomTermForCovars <- setdiff(unique(unlist(randomTerm)),c("environment","designation"))
  fixedTermForCovars <- setdiff(unique(unlist(fixedTerm)),c("environment","designation"))

  #if(covars!="none"){
    if( any( covars %in% c("Relationship structure_Geno","Relationship structure_Pedigree","Relationship structure_GenoAD") ) ){
      if(verbose){message("Checking and calculating kernels requested")}
      
      if(any(c("Relationship structure_Geno","Relationship structure_GenoAD")  %in% covars)){
        if(verbose){message(paste(" Marker kernel is requested"))}
        #qas <- which( phenoDTfile$status$module == "qaGeno" )
        #idqaGeno<-phenoDTfile$status$analysisId[qas]
        qas<-which( names(phenoDTfile$data$geno_imp)==analysisIdgeno )
        ## MARKER KERNEL
        Markers <- as.data.frame(phenoDTfile$data$geno_imp[qas]) # in form of covariates
       
        if(length(qas) > 0){
          if(length(which(is.na(Markers))) > 0){Markers <- apply(Markers,2,sommer::imputev)}
        }else{
          missing <- apply(Markers,2,sommer::propMissing)
          Markers <- apply(Markers[,which(missing < 0.9)],2,sommer::imputev)
        }

        ploidyFactor <- max(Markers)/2
        mydata <- phenoDTfile$predictions #
        mydata <- mydata[which(mydata$analysisId %in% analysisId),]
        
        if("Relationship structure_Geno" %in% covars){
          #G <- sommer::A.mat(Markers-ploidyFactor);
          males <- as.character(unique(mydata$designation))
          Markers <- Markers[intersect(rownames(Markers),c(males)),]
          G <- sommer::A.mat(as.matrix(Markers)-ploidyFactor)
          missing <- setdiff(c(males),rownames(G))
          A1m <- diag(mean(diag(G)),nrow = length(missing),ncol=length(missing))
          rownames(A1m) <- colnames(A1m) <- missing
          G <- sommer::adiag1(G,A1m) 
          G <- G + diag(1e-5,ncol(G),ncol(G))
          } # additive model
        if("Relationship structure_GenoAD" %in% covars){ # if genetic value is desired let's do a log marker model
          Markers <- apply(Markers+1,2,log)
          Gad <- sommer::A.mat(as.matrix(Markers))
          Gad <- Gad + diag(1e-5, ncol(Gad), ncol(Gad))
        } # additive + dominance model
      }
      
      ### WEATHER KERNEL
      #if("weather" %in% covars & !is.null(Weather)){
      #  #classify <- unique(unlist(randomTerm)[which(unlist(expCovariates) %in% "weather")])
      #  if(verbose){message(paste("   Weather kernel is requested"))}
      #  WeatherK <- Weather
      #  rownamesWeather <- rownames(WeatherK)
      #  WeatherK <- apply(WeatherK, 2, scale)
      #  WeatherK <- WeatherK[,which( !is.na(apply(WeatherK,2,var)) ), drop=FALSE]
      #  rownames(WeatherK) <- rownamesWeather
      #  W <- sommer::A.mat(WeatherK)
      #  W <- W + diag(1e-5, ncol(W), ncol(W))
      #  Wchol <- t(chol(W))
      #}
      
      ## PEDIGREE KERNEL
      if("Relationship structure_Pedigree" %in% covars){
        #classify <- unique(unlist(randomTerm)[which(unlist(expCovariates) %in% "pedigree")])
        if(verbose){message(paste("   Pedigree kernel is requested"))}
        Pedigree <- phenoDTfile$data$pedigree[,2:4]
        mydataX <-  phenoDTfile$predictions[which( phenoDTfile$predictions$analysisId %in% analysisId),]
        Pedigree<-Pedigree[!duplicated(Pedigree[,1]),]
        Pedigree <- Pedigree[which(Pedigree[,1] %in% unique(mydataX$designation) ), ]
        N <- asreml::ainverse(Pedigree)
      }
    }
 #}


  ##########################################
  ##########################################
  ## COMPLETE THE CLEANING PARAMETERS (7 lines)
  names(traitFamily) <- trait
  heritLB <- rep(heritLB,length(trait)); heritLB <- heritLB[1:length(trait)]; names(heritLB) <- trait
  heritUB <- rep(heritUB,length(trait)); heritUB <- heritUB[1:length(trait)]; names(heritUB) <- trait
  meanLB <- rep(meanLB,length(trait)); meanLB <- meanLB[1:length(trait)]; names(meanLB) <- trait
  meanUB <- rep(meanUB,length(trait)); meanUB <- meanUB[1:length(trait)]; names(meanUB) <- trait
  traitOrig <- trait ## ?????????
  #if(length(fixedTerm) == 0 | is.null(fixedTerm)){fixedTerm <- list("1")} # assign the intercept if there's no fixed effects
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
  myDataTraits <- fixedTermTrait <- randomTermTrait <- groupingTermTrait <- Mtrait <- envsTrait <- entryTypesTrait <- list()
  randomTermModel<-fixedTermModel<-list()
  x_option <- function(x,y,z) {
    switch(x,
           "none" = paste0("idv(",y,")"),
           "Structure model_fa" = paste0("fa(",y,",",z,")"),
           "Structure model_diag" = paste0("diag(",y,")"),
           "Structure model_us" = paste0("us(",y,")"),
           "Relationship structure_Geno"= paste0("vm(",y,",G)"),
           "Relationship structure_Pedigree"= paste0("vm(",y,",N)"),
           "Relationship structure_GenoAD"= paste0("vm(",y,",Gad)"),
           stop("Invalid `x` value")
    )
  }
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
      if(nrow(prov) > 0){ # if after filters there's still data for this trait we can continue and save the data
        if( var(prov[,"predictedValue"], na.rm = TRUE) > 0 ){ # check that there is variance
          # make new formula for this specific trait if it passed all the filters

          #fixed formula per trait
          fixedTermTrait[[iTrait]]<-paste("asreml::asreml(fixed = predictedValue ~ 1")
          fixedTermprov=unique(unlist(lapply(fixedTerm,paste,collapse=":")))
          fixedTermModel[[iTrait]]=fixedTermprov
          if(length(fixedTermprov) != 0 | !is.null(fixedTermprov)){
            for(iFixed in 1:length(fixedTermprov)){ # for each element in the list # iFixed=1
              fixedTermTrait[[iTrait]]<-paste(fixedTermTrait[[iTrait]],fixedTermprov[iFixed],sep=" + ")
            }
          }
          fixedTermTrait[[iTrait]]<-paste0(fixedTermTrait[[iTrait]],",")

          # random formula per trait
          randomTermTrait[[iTrait]]<-paste("random = ~")
          randomTermprov=list()
          for(w in 1:length(randomTerm)){
            if(length(covMod[[w]])==1){
              randomTermprov[[w]]=x_option(covMod[[w]],randomTerm[[w]],nFA[[w]])
            }else{
              randomTermprov[[w]]= c(x_option(covMod[[w]][1],randomTerm[[w]][1],nFA[[w]]),x_option(covMod[[w]][2],randomTerm[[w]][2],nFA[[w]]))
            }
          }
          randomTermprov=unique(unlist(lapply(randomTermprov,paste,collapse=":")))
          randomTermModel[[iTrait]]=randomTermprov
          if(length(randomTermprov) != 0 | !is.null(randomTermprov)){
            for(iRand in 1:length(randomTermprov)){ # for each element in the list # iFixed=1
              randomTermTrait[[iTrait]]<-paste(randomTermTrait[[iTrait]],randomTermprov[iRand],sep=" + ")
            }
          }
          randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],",")

          #if(covMod=="us_model" & addG=="none"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"us(environment):idv(designation),");intTerm="environment:designation"}
          #if(covMod=="us_model" & addG=="genoAD_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"us(environment):vm(designation,G),");intTerm="environment:vm(designation,G)"}
          #if(covMod=="us_model" & addG=="geno_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"us(environment):vm(designation,G),");intTerm="environment:vm(designation,G)"}
          #if(covMod=="us_model" & addG=="pedigree_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"us(environment):vm(designation,N),");intTerm="environment:vm(designation,N)"}

          #if(covMod=="diag_model" & addG=="none"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"diag(environment):idv(designation),");intTerm="environment:designation"}
          #if(covMod=="diag_model" & addG=="genoAD_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"diag(environment):vm(designation,G),");intTerm="environment:vm(designation,G)"}
          #if(covMod=="diag_model" & addG=="geno_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"diag(environment):vm(designation,G),");intTerm="environment:vm(designation,G)"}
          #if(covMod=="diag_model" & addG=="pedigree_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"diag(environment):vm(designation,N),");intTerm="environment:vm(designation,N)"}

          #if(covMod=="fa_model" & addG=="none"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],paste0("fa(environment,",nFA,"):idv(designation),"));intTerm=paste0("fa(environment,",nFA,"):designation")}
          #if(covMod=="fa_model" & addG=="genoAD_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],paste0("fa(environment,",nFA,"):vm(designation,G),"));intTerm=paste0("fa(environment,",nFA,"):vm(designation,G)")}
          #if(covMod=="fa_model" & addG=="geno_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],paste0("fa(environment,",nFA,"):vm(designation,G),"));intTerm=paste0("fa(environment,",nFA,"):vm(designation,G)")}
          #if(covMod=="fa_model" & addG=="pedigree_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],paste0("fa(environment,",nFA,"):vm(designation,N),"));intTerm=paste0("fa(environment,",nFA,"):vm(designation,N)")}

          #if(covMod=="none" & addG=="none"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"environment:designation,");intTerm="environment:designation"}
          #if(covMod=="none" & addG=="genoAD_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"idv(environment):vm(designation,G),");intTerm="environment:vm(designation,G)"}
          #if(covMod=="none" & addG=="geno_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"idv(environment):vm(designation,G),");intTerm="environment:vm(designation,G)"}
          #if(covMod=="none" & addG=="pedigree_model"){randomTermTrait[[iTrait]]<-paste0(randomTermTrait[[iTrait]],"idv(environment):vm(designation,N),");intTerm="environment:vm(designation,N)"}

          myDataTraits[[iTrait]] <- prov # dataset for this trait
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
    fixedTermSub<-fixedTermTrait[[iTrait]]
    randomTermSub<-randomTermTrait[[iTrait]]
    ## deregress if needed
    VarFull <- var(mydataSub[,"predictedValue"], na.rm = TRUE) # total variance
    effectTypeTrait <- phenoDTfile$modeling[which(phenoDTfile$modeling$analysisId == analysisId & phenoDTfile$modeling$trait == iTrait & phenoDTfile$modeling$parameter == "designationEffectType"),"value"]
    if(names(sort(table(effectTypeTrait), decreasing = TRUE))[1] == "BLUP"){ # if STA was BLUPs deregress
      mydataSub$predictedValue <- mydataSub$predictedValue/mydataSub$reliability
    }
    ## calculate weights
    mydataSub=mydataSub[with(mydataSub, order(environment)), ] # sort by environments
    mydataSub$w <- 1/(mydataSub$stdError^2) # add weights column
    # warnin messages in weights use
      if(verbose){message("   Using weights in the analysis. Residual variance will be fixed to 1.")  }
    #trasform to factor variables
    labfactor<-c("environment","module","analysisId","pipeline","trait","gid","designation","mother","father","entryType","effectType","year","location")
    tranfactfixed<-which(unlist(fixedTerm)%in%labfactor==T)
    tranfactrandom<-which(unique(unlist(randomTerm))%in%labfactor==T)
    #if(length(tranfactfixed)!=0){mydataSub[,unlist(fixedTerm)[tranfactfixed]]=as.factor(mydataSub[,unlist(fixedTerm)[tranfactfixed]])}
    #if(length(tranfactrandom)!=0){mydataSub[,unlist(randomTerm)[tranfactrandom]]=as.factor(data.frame(mydataSub[,unlist(randomTerm)[tranfactrandom]]))}
    if(length(tranfactfixed)!=0){mydataSub[,unlist(fixedTerm)[tranfactfixed]]=lapply(mydataSub[unlist(fixedTerm)[tranfactfixed]], as.factor)}
    if(length(tranfactrandom)!=0){mydataSub[,unique(unlist(randomTerm))[tranfactrandom]]=lapply(mydataSub[unlist(randomTerm)[tranfactrandom]], as.factor)}

    asreml::asreml.options(workspace = 30e7, pworkspace = 200e7, trace=T,ai.sing=T)
    ## model fit
    modelfit=paste0("mix<-",fixedTermSub, randomTermSub,"weights = w,family = ",traitFamily[iTrait],	",na.action = na.method(x='include'),data=mydataSub,maxit = ",maxIters,")")
	  mix <- try(eval(parse(text=modelfit)),silent=TRUE)

	if(!inherits(mix,"try-error") ){
	  upval=0
	  if (!mix$converge & upval<=10) { upval=upval+1; eval(parse(text='mix<-asreml::update.asreml(mix)')) }
	}
    # print(mix$VarDf)
    pp <- list()
    if(!inherits(mix,"try-error") ){
      # get variance components
      ss <- summary(mix)$varcomp
      Ve <- ss["units!R","component"]
      mu<-coef(mix)[["fixed"]]["(Intercept)",]
      Ci<-mix[["vcoeff"]][["fixed"]][1]*mix$sigma2
      #mu <- mix$coefMME[mix$ndxCoefficients$`(Intercept)`]
      #Ci <- solve(mix$C)
      if(length(mu) > 0){
        pp[["(Intercept)"]] <- data.frame(designation="(Intercept)", predictedValue=mu, stdError=sqrt(Ci), reliability=NA,
                                          trait=iTrait, effectType="(Intercept)", entryType="(Intercept)", environment="(Intercept)" )
      }
      for(iGroupFixed in unlist(fixedTermModel[[iTrait]])){ # iGroupFixed=unlist(fixedTermModel[[iTrait]])[1]
        pick <- coef(mix)[["fixed"]][grepl(iGroupFixed,rownames(coef(mix)[["fixed"]]))]
        if (length(pick)>1){
          names(pick)<-rownames(coef(mix)[["fixed"]])[grepl(iGroupFixed,rownames(coef(mix)[["fixed"]]))]
          pick <- pick[which(pick!=0)]
        }else{
          if(pick==0){
            pick<-rep((-1)*mu,2)
            names(pick)<-rep(rownames(coef(mix)[["fixed"]])[grepl(iGroupFixed,rownames(coef(mix)[["fixed"]]))],2)
            pick<-pick[1]
          }else{
            pick<-rep(pick,2)
            names(pick)<-rep(rownames(coef(mix)[["fixed"]])[grepl(iGroupFixed,rownames(coef(mix)[["fixed"]]))],2)
            pick<-pick[1]
          }
        }
        # shouldBeOne <- which(pick == 0)
        # if(length(shouldBeOne) > 0){pick[shouldBeOne] = 1}
        blue <- pick+mu; names(blue) <- names(pick); #blue[1] <- blue[1]-mu
        #start <- sum(mix$EDdf[1:(which(mix$EDdf$Term == iGroupFixed) - 1),"Model"]) # we don't add a one because we need the intercept
        #nEffects <- mix$EDdf[which(mix$EDdf$Term == iGroupFixed),"Effective"]#length(blue)
        #pev <- Ci[start:(start+nEffects-1),start:(start+nEffects-1)]
        pev <- mix[["vcoeff"]][["fixed"]][grepl(iGroupFixed,rownames(coef(mix)[["fixed"]]))]*mix$sigma2
        if (length(pick)>=1){pev <- pev[which(pev!=0)]}
        if(is.matrix(pev)){ stdError <- (sqrt(Matrix::diag(pev)))}else{stdError <- sqrt(pev)}
        prov <- data.frame(designation=names(blue), predictedValue=blue, stdError=stdError, reliability=NA,
                           trait=iTrait, effectType=iGroupFixed, environment="(Intercept)" )
        for(iLabel in unique(unlist(fixedTerm))){
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
      }

      subgroupGen=c() 
      only1=which(lapply(randomTerm,length)==1)
      for (sb in 1:length(only1)){
        subgroupGen[sb]=randomTerm[[only1[sb]]]
      }
      subgroupGen=subgroupGen[which(subgroupGen%in%c("designation","mother","father","gid")==T)]
      
      x_option2 <- function(x,y,z) {
        switch(x,
               "none" = paste0(y),
               "Structure model_fa" = paste0("fa(",y,",",z,")"),
               "Structure model_diag" = paste0(y),
               "Structure model_us" = paste0(y),
               "Relationship structure_Geno"= paste0("vm(",y,",G)"),
               "Relationship structure_Pedigree"= paste0("vm(",y,",N)"),
               "Relationship structure_GenoAD"= paste0("vm(",y,",Gad)"),
               stop("Invalid `x` value")
        )
      }
      
      subgroupInt=c() 
      only2=which(lapply(randomTerm,length)==2)
      for (sb2 in 1:length(only2)){
        ed=all(randomTerm[[only2[sb2]]]%in%c("environment","designation"))
        if(ed==T){subgroupInt[sb2]=paste0(x_option2(covMod[[only2[sb2]]][1],randomTerm[[only2[sb2]]][1],nFA[[only2[sb2]]]),":",x_option2(covMod[[only2[sb2]]][2],randomTerm[[only2[sb2]]][2],nFA[[only2[sb2]]]))}
      }
      
      if(gsca==T){
        gca=summary(mix,coef=TRUE)$coef.random
        prov <- data.frame(designation=rownames(gca), predictedValue=gca[,1], stdError=gca[,2], reliability=rep(NA,dim(gca)[1]),
                           trait=iTrait, effectType="GCA/SCA" , environment="GCA/SCA", entryType="unknown")
        pp[["GCA/SCA"]] <- prov
        }
      
      groupingSub=c(subgroupGen,subgroupInt)
      for(iGroup in groupingSub){
        #for( iGroup in names(groupingSub)){ # iGroup=groupingSub[2]
        if(iGroup%in%subgroupGen){
          Vg <- ss[iGroup,"component"]
        }else{Vg<-NA}
        
        blup=predict(mix,classify=iGroup)$pvals
        if (all(blup$status=="Aliased")){
          statusmetrics="Aliased estimation, problems with the model, please check!"
          stdError <- reliability <- rep(NA,dim(blup)[1])
        }else{
          statusmetrics=mix$converge
          message(paste("   Calculating standar errors for",iTrait, iGroup,"predictions"))
          stdError <- blup$std.error # random effect was just one column
          reliability <- abs((Vg - (stdError^2)) /Vg) # reliability <- abs((Vg - Matrix::diag(pev))/Vg)
        }

        badRels <- which(reliability > 1); if(length(badRels) > 0){reliability[badRels] <- 0.9999}
        badRels2 <- which(reliability < 0); if(length(badRels2) > 0){reliability[badRels2] <- 0}
        if(iGroup==subgroupInt){effTypeSub="environment_designation"}else{effTypeSub=iGroup}
        if(iGroup%in%c("designation","mother","father","gid")){envTypeSub="(Intercept)"}else{envTypeSub=blup$environment}

        if(ncol(blup)>4){
          namesblup=c(iGroup,names(blup)[3:5])
          blup=data.frame(cbind(paste0(blup[,1],":",blup[,2]),blup[,3],blup[,4],blup[,5]))
          names(blup)=namesblup
          blup$predicted.value=as.numeric(blup$predicted.value)
        }

          prov <- data.frame(designation=blup[,iGroup], predictedValue=blup$predicted.value, stdError=stdError, reliability=reliability,
                             trait=iTrait, effectType=effTypeSub , environment=envTypeSub)

          # end of adding fixed effects
          sdP <- sd(prov[,"predictedValue"],na.rm=TRUE)
          cv <- (sd(prov[,"predictedValue"],na.rm=TRUE)/mean(prov[,"predictedValue"],na.rm=TRUE))*100
          # add additional entry type labels
          mydataSub[,"designationXXX"] <- apply(mydataSub[,"designation",drop=FALSE],1,function(x){paste(x,collapse = ":")})
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
                                       data.frame(module="mtaAsr",analysisId=mtaAnalysisId, trait= iTrait,
                                                  environment=paste(unique(envTypeSub), collapse = "_"),
                                                  parameter=c( paste(c("mean","sd", "r2","Var"),effTypeSub,sep="_") ),
                                                  method=c("sum(x)/n","sd","(G-PEV)/G","REML"),
                                                  value=c(mean(prov[,"predictedValue"], na.rm=TRUE), sdP, median(reliability), var(prov[,"predictedValue"], na.rm=TRUE) ),
                                                  stdError=c(NA,NA,sd(reliability, na.rm = TRUE)/sqrt(length(reliability)),NA )
                                       )
          )
      }

      kernels=paste(unlist(addG)[which(unlist(addG)%in%"Relationship structure_Geno","Relationship structure_Pedigree","Relationship structure_GenoAD"==T)],sep=",")
      if(length(kernels)==0){kernels="none"}
      ## save the modeling used
      currentModeling <- data.frame(module="MtaAsr", analysisId=mtaAnalysisId,trait=iTrait, environment=c(rep("across",4), "designation"),
                                    parameter=c("fixedFormula","randomFormula","family","convergence",rep("kernels",1)),
                                    value=c(fixedTermSub,randomTermSub,traitFamily[iTrait],statusmetrics,kernels ))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
      ## save the environments used goodFields
      currentModeling <- data.frame(module="MtaAsr", analysisId=mtaAnalysisId,trait=iTrait, environment=allEnvironments,
                                    parameter="includedInMta",
                                    value=ifelse(allEnvironments%in%unique(mydataSub$environment), TRUE, FALSE))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )

    }else{ # if model failed
      if(verbose){ cat(paste("Mixed model failed for trait",iTrait,". Aggregating and assuming h2 = 0 \n"))}
      means <- aggregate(predictedValue ~ designation, FUN=mean, data=mydataSub)
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
                                   data.frame(module="mtaAsr",analysisId=mtaAnalysisId, trait=iTrait,
                                              environment="across",
                                              parameter=c("mean","sd", "r2","Var_designation","Var_residual"),
                                              method=c("sum(x)/n","sd","(G-PEV)/G","REML","REML"),
                                              value=c(mean(means$predictedValue, na.rm=TRUE), sdP, NA, NA, NA ),
                                              stdError=NA
                                   )
      )
      currentModeling <- data.frame(module="mtaAsr", analysisId=mtaAnalysisId,trait=iTrait, environment="across",
                                    parameter=c("fixedFormula","randomFormula","family","designationEffectType"),
                                    value=c("None","None","None","mean"))
      phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
      pp[["designation"]] <- means
    }

    predictionsTrait <- do.call(rbind,pp)
    #############################################################
    ## add across env estimate for DESIGNATION effect type fitted
    #############################################################
    match1 <- unlist(lapply(unlist(fixedTerm),function(x){sum(as.numeric(x=="designation"))}))
    names(match1) <- unlist(lapply(unlist(fixedTerm),function(x){paste(x,collapse = "_")}))
    match2 <- unlist(lapply(unlist(randomTerm),function(x){sum(as.numeric(x=="designation"))}))
    names(match2)<-unlist(lapply(unlist(randomTerm),function(x){paste(x,collapse = "_")}))
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
    match1 <- unlist(lapply(unlist(fixedTerm),function(x){sum(as.numeric(x=="gid"))}))
    names(match1) <- unlist(lapply(unlist(fixedTerm),function(x){paste(x,collapse = "_")}))
    match2 <- unlist(lapply(unlist(randomTerm),function(x){sum(as.numeric(x=="gid"))}))
    names(match2)<-unlist(lapply(unlist(randomTerm),function(x){paste(x,collapse = "_")}))
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
    phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                 data.frame(module="mtaAsr",analysisId=mtaAnalysisId, trait= iTrait, environment="across",
                                            parameter=c("Var_residual","nEnv","nEntries"),
                                            method=c("REML","n","n"),
                                            value=c( Ve, length(goodFields), length(unique(mydataSub$designation)) ),
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
  predictionsBind$module <- "mtaAsr"; rownames(predictionsBind) <- NULL
  # print(head(predictionsBind))
  #########################################
  ## update databases
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(!is.null(phenoDTfile$predictions)){
    if("effectType" %!in% colnames(phenoDTfile$predictions) ){
      phenoDTfile$predictions$effectType <- NA
    }
  }
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions,
                                   predictionsBind[,colnames(phenoDTfile$predictions)])
  newStatus <- data.frame(module="mtaAsr", analysisId=mtaAnalysisId, analysisIdName=NA)
  phenoDTfile$status <- rbind( phenoDTfile$status, newStatus[,colnames(phenoDTfile$status)] )
  ## add which data was used as input
  modeling <- data.frame(module="mtaAsr",  analysisId=mtaAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId ))

  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)
}
