resMod <- function(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, # per trait
    genoUnit = c("designation"),
    addRes=FALSE,
    verbose=TRUE
){
  staAnalysisId <- as.numeric(Sys.time())
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  fixedTerm <- unique(c("1", "rep","iBlock", genoUnit))
  
  ###################################
  # loading the dataset
  mydata <- phenoDTfile$data$pheno # extract relevant data for sta
  originalColumns <- colnames(mydata)
  myped <- phenoDTfile$data$pedigree
  ### change column names for mapping
  paramsPheno <- phenoDTfile$metadata$pheno
  paramsPheno <- paramsPheno[which(paramsPheno$parameter != "trait"),]
  colnames(mydata) <- cgiarBase::replaceValues(colnames(mydata), Search = paramsPheno$value, Replace = paramsPheno$parameter )
  paramsPed <- phenoDTfile$metadata$pedigree
  # make sure there is no columns called mother and father in the phenotype dataset
  badPaternalCols <- which(colnames(mydata) %in% c("mother","father"))
  if(length(badPaternalCols) > 0){mydata <- mydata[,-badPaternalCols, drop=FALSE]}
  motherColumn <- which(paramsPed$parameter == "mother")
  if(length(motherColumn) > 0){
    if(paramsPed$value[motherColumn] == ""){
      myped$mother <- NA; paramsPed$value[motherColumn] <- "mother"
    }
  }else{
    myped$mother <- NA; paramsPed$value[motherColumn] <- "mother"
    paramsPed <- rbind(paramsPed, data.frame(parameter="mother",value="mother"))
  }
  fatherColumn <- which(paramsPed$parameter == "father")
  if(length(fatherColumn) > 0){
    if(paramsPed$value[fatherColumn] == ""){
      myped$father <- NA; paramsPed$value[fatherColumn] <- "father"
    }
  }else{
    myped$father <- NA; paramsPed$value[fatherColumn] <- "father"
    paramsPed <- rbind(paramsPed, data.frame(parameter="father",value="father"))
  }

  
  myped <- cgiarBase::nrm2(pedData= myped, verbose=FALSE,returnMatrix=FALSE,
                           indivCol = paramsPed[paramsPed$parameter=="designation","value"],
                           damCol = paramsPed[paramsPed$parameter=="mother","value"],
                           sireCol = paramsPed[paramsPed$parameter=="father","value"]
  )
  colnames(myped) <- c("designation","mother","father")
  ### make sure all expected columns are present
  required_mapping <- c("stage", "pipeline", "country", "year", "season", "location", "trial", "environment", "rep", "iBlock", "row", "col", "designation", "gid", "entryType", "trait")
  for(iRequired in required_mapping){
    if(iRequired %in% colnames(mydata)){}else{mydata[,iRequired] <- NA}
  }
  if (nrow(mydata) < 2) stop("Not enough phenotypic data is available to perform a single trial analysis. Please add the phenotypic data to your data object.", call. = FALSE)
  if( length(setdiff(setdiff(fixedTerm,"1"),c( colnames(mydata), colnames(myped) ) )) > 0 ){stop(paste("column(s):", paste(setdiff(setdiff(fixedTerm,"1"),colnames(mydata)), collapse = ","),"couldn't be found."), call. = FALSE)}
  mydata$rowindex <- 1:nrow(mydata)

  # merge mother and father information
  if(!is.null(myped)){
    if(nrow(myped) > 0){
      mydata <- merge(mydata, myped, by="designation", all.x = TRUE)
      mydata <- mydata[with(mydata, order(rowindex)), ]
    }else{mydata$mother <- NA; mydata$father <- NA}
  }else{mydata$mother <- NA; mydata$father <- NA}
  # move the genotype columns to factor

  if(is.null(analysisId)){ # user doesn't want to use modifications
    stop("Please provide an analysisId from phenotype QA", call. = FALSE)
  }else{
    cleaning <- phenoDTfile$modifications$pheno # extract outliers
    cleaning <- cleaning[which(cleaning$analysisId %in% analysisId),]
  }
  # remove traits that are not actually present in the dataset
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% colnames(mydata)){
      if(verbose){
        cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))
      }
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  if(length(trait)==0){stop("None of the traits specified are available. Please double check", call. = FALSE)}
  traitTypes <- unlist(lapply(mydata, class))[trait]
  numericTraits <- which(traitTypes %in% c("numeric","integer"))
  if(length(numericTraits)==0){stop("None of the traits specified are numeric in nature. Please double check", call. = FALSE)}
  trait <- trait[numericTraits]
  
  # check if experimental design factor filtering is done
  if(nrow(cleaning[which(cleaning$module == "qaDesign"),] > 0)){
    cleaningSubDes <- cleaning[which(cleaning$trait %in% c("rep","iBlock")),]
    
    for (iDes in unique(cleaningSubDes$trait)){
      outDes <- which(mydata$rowindex %in% cleaningSubDes[which(cleaningSubDes$trait == iDes), "row"])
      mydata[outDes,iDes] <- NA
    }
  }
  
  if(addRes==TRUE){
    modifications <- phenoDTfile$modifications$pheno
    idQa <- phenoDTfile$status[which(phenoDTfile$status$module %in% c("qaMb")),"analysisId"];
    idQa <- idQa[length(idQa)]
    idQaType <- phenoDTfile$status[phenoDTfile$status$analysisId == idQa, "module"]
    modifications <- modifications[which( modifications$analysisId == idQa & modifications$module == idQaType),]
    traitsBoxplot <- unique(modifications[which( modifications$analysisId == idQa),"trait"])
    for(iTrait in traitsBoxplot){#iTrait=traitsBoxplot[1]
      if(nrow(modifications) > 0){
        mydata[which(mydata$rowindex %in% unique(modifications$row[which(modifications$trait == iTrait)])),iTrait]=NA
      }
    }
  }
  #####################################
  # single trial analysis
  fixedFormulaForFixedModel <- randomFormulaForFixedModel <- NULL
  # fields <- as.character(na.omit(unique(mydata$environment)))

  designationColumns <- paramsPed[which(paramsPed$parameter %in%  genoUnit),"value"]
  fieldsL <- list()
  for(igenoUnit in genoUnit){
    fieldsL[[igenoUnit]] <- names(which(apply(table(mydata[,"environment"],mydata[,igenoUnit]),1,sum)>3))
  }
  fields <- Reduce(intersect,fieldsL)


  if(length(fields) == nrow(mydata)){
    stop("The number of environment levels is equal to the number of records/rows in the dataset.
          This means that probably you didn't select the right columns to define the environment column.
          Please match again your raw file checking carefully the columns defining the environment.", call.=FALSE )
  }
  predictionsList <- list()
  counter=1
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    
    for(iField in fields){ # iField = fields[1]# "ARH1_2016"
      if(verbose){cat(paste("Analyzing field", iField,"\n"))}
      # subset data
      mydataSub <- droplevels(mydata[which(as.character(mydata$environment) %in% iField),])
      mydataSub$trait <- as.numeric(mydataSub[,iTrait])
      # make factors
      for(iEd in c("environment","rep","iBlock")){
        #if(iEd %in% c("row","col")){mydataSub[,iEd] <- as.numeric(mydataSub[,iEd])}
        # if(iEd %in% colnames(mydataSub)){
        mydataSub[,paste0(iEd,"F")] <-  as.factor(mydataSub[,iEd])
        # }else{  mydataSub[,paste0(iEd,"F")] <- NA; mydataSub[,iEd] <- NA   }
      }
      for(iName in c("designation","mother","father")){
        mydataSub[,iName] <- as.factor(mydataSub[,iName])
      }
      # check the genetic units
      nLevelsGenounit <- apply(data.frame(genoUnit),1,function(x){length(table(mydataSub[,x])) }); names(nLevelsGenounit) <- genoUnit
      genoUnitTraitField <- names(nLevelsGenounit)[which(nLevelsGenounit > 1)]
      if(length(genoUnitTraitField)==0){
        warning(paste("There is not enough levels in your specified genetic unit(s):", paste(genoUnit, collapse = ", "), "in environment", iField, ". Ignoring environment." ),call. = FALSE)
      }else{
        nLevelsFixedunit <- apply(data.frame(setdiff(fixedTerm,"1")),1,function(x){length(table(mydataSub[,x])) })
        names(nLevelsFixedunit) <- setdiff(fixedTerm,"1")
        badLevelsFixedunit <- names(nLevelsFixedunit)[which(nLevelsFixedunit <= 1)]
        fixedTermTraitField <- unique(c("1",setdiff(fixedTerm,badLevelsFixedunit)))
        # impute fixed effect columns if they are numeric
        toImpute <- unlist(lapply(mydataSub[,setdiff(fixedTermTraitField,"1")], class))
        keepToImpute <- which(toImpute %in% "numeric")
        if(length(keepToImpute) > 0){
          toImpute <- names(toImpute[keepToImpute])
          for(iImpute in toImpute){mydataSub[, iImpute] <- enhancer::imputev(mydataSub[, iImpute])}
        }
        # do analysis
        if(!is.na(var(mydataSub[,"trait"],na.rm=TRUE))){ # if there's variance
          if( var(mydataSub[,"trait"], na.rm = TRUE) > 0 ){
            # find best experimental design formula
            randomTerms <- c("repF", "iBlockF")
            min_levels <- c(repF = 2, iBlockF = 4)
            
            screened <- cgiarBase::screen_sta_random_effects(
              random_effects = randomTerms,
              dat = mydataSub,
              min_levels = min_levels
            )
            screened$summary
            newRandom <- if(length(screened$kept) > 0) screened$kept else NULL

            for(iGenoUnit in genoUnitTraitField){ # iGenoUnit <- genoUnitTraitField[1]
              myGeneticUnit <-  iGenoUnit
              fixF<-paste(iTrait," ~ 1 +", iGenoUnit)
              if(length(newRandom)>1){fixF2<-paste(" + repF + repF:iBlockF")}else{fixF2<-paste(" + ",newRandom)}
              mod1<-try(lm(as.formula(paste(fixF,fixF2)),data=mydataSub,na.action=na.exclude),silent=TRUE)              #hasta aqui todo corre
              if(!inherits(mod1,"try-error") ){
                datmp=data.frame(mydataSub[,c("rowindex","trait","designation","repF","iBlockF")])
                names(datmp)[1]="idRowR"
                names(datmp)[2]="Observed"
                #if(length(newRandom)==1){datmp$iBlockF=NA}
                pp<-cbind(datmp,trait=iTrait,environment=iField,fitted=fitted(mod1),residuals=residuals(mod1),cookD=cooks.distance(mod1))
                #pp$thr<-4/nrow(pp)
                pp$leverage <- hatvalues(mod1)
                pp$stdt<- rstudent(mod1)
                pp$N<- nobs(mod1)
                pp$p <- length(coef(mod1))
                pp$outlier="GOOD"
                criterios <- (pp$cook > (4/nobs(mod1))) + (abs(pp$stdt) > 3.5) +(pp$leverage > (2*length(coef(mod1)))/nobs(mod1) )
                pp$outlier[criterios >= 2] <- "outlier"
                pp$outlierFinal<-pp$outlier
                predictionsList[[counter]]<-pp
                counter=counter+1
              }
            }
          }
        }
      }
    }
  }
  resTmp<-do.call(rbind,predictionsList)
  resTmp[which(resTmp$cookD=="NaN"),"cookD"]<-NA
  return(resTmp)
}
