staLMM <- function(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, # per trait
    traitFamily=NULL,
    fixedTerm=c("1"),
    genoUnit = c("designation"),
    maxit=50,
    returnFixedGeno=TRUE,
    verbose=TRUE
){
  ## THIS FUNCTION CALCULATES A SINGLE TRIAL ANALYSIS FOR MULTIPLE FIELDS AND TRAITS
  ## IS USED IN THE BANAL APP UNDER THE GENETIC EVALUATION MODULES
  staAnalysisId <- as.numeric(Sys.time())#gsub(" ","-",gsub("[[:punct:]]", "-", as.character(Sys.time())) )
  if(is.null(phenoDTfile)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  # if(is.null(analysisId)){stop("Please provide analysisIds to be considered in cleaning", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  if(is.null(traitFamily)){traitFamily <- rep("gaussian(link = 'identity')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  names(traitFamily) <- trait
  # fixedTerm <- unique(c("1",fixedTerm,"designation"))
  fixedTerm <- unique(c("1", fixedTerm, genoUnit))

  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(is.null(phenoDTfile$metrics)){
    provMet <- as.data.frame(matrix(nrow=0, ncol=7))
    colnames(provMet) <- c("module","analysisId","trait","environment","parameter","value","stdError")
    phenoDTfile$metrics <- provMet
  }
  if(is.null(phenoDTfile$predictions)){
    provMet <- as.data.frame(matrix(nrow=0, ncol=13))
    colnames(provMet) <- c("module"   ,      "analysisId" ,    "pipeline"   ,    "trait"   ,       "gid"  ,          "designation"   ,"mother"  ,       "father"    ,     "entryType" ,     "environment" ,   "predictedValue" ,"stdError"  ,  "reliability"  )
    phenoDTfile$predictions <- provMet
  }
  if(is.null(phenoDTfile$modeling)){
    provMet <- as.data.frame(matrix(nrow=0, ncol=6))
    colnames(provMet) <- c("module" ,     "analysisId"  ,"trait"  ,     "environment" ,"parameter"  , "value" )
    phenoDTfile$modeling <- provMet
  }
  if(is.null(phenoDTfile$status)){
    provMet <- as.data.frame(matrix(nrow=0, ncol=3))
    colnames(provMet) <- c("module" ,     "analysisId" , "analysisIdName" )
    phenoDTfile$status <- provMet# data.frame(module="sta", analysisId=staAnalysisId, analysisIdName=NA)
  }
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
  # colnames(myped) <- cgiarBase::replaceValues(colnames(myped), Search = paramsPed$value, Replace = paramsPed$parameter )
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

  # if(length(intersect(paramsPed$value, colnames(myped)))  <  3 ){
  #   stop("Metadata for pedigree (mapping) and pedigree frame do not match. Please reupload and map your pedigree information.", call. = FALSE)
  # }

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
    cleaningSubDes <- cleaning[which(cleaning$trait %in% c("row","col","rep","iBlock")),]
    
    for (iDes in unique(cleaningSubDes$trait)){
      outDes <- which(mydata$rowindex %in% cleaningSubDes[which(cleaningSubDes$trait == iDes), "row"])
      mydata[outDes,iDes] <- NA
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
  predictionsList <- list(); columnsToAdd <- character(); counter=1
  library(LMMsolver)
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    mydata[,paste(iTrait,"residual",sep="-")] <- NA

    for(iField in fields){ # iField = fields[1]# "ARH1_2016"
      if(verbose){cat(paste("Analyzing field", iField,"\n"))}
      # subset data
      mydataSub <- droplevels(mydata[which(as.character(mydata$environment) %in% iField),])
      mydataSub$trait <- as.numeric(mydataSub[,iTrait])
      # make factors
      for(iEd in c("environment","trial","row","col","rep","iBlock")){
        if(iEd %in% c("row","col")){mydataSub[,iEd] <- as.numeric(mydataSub[,iEd])}
        # if(iEd %in% colnames(mydataSub)){
        mydataSub[,paste0(iEd,"F")] <-  as.factor(mydataSub[,iEd])
        # }else{  mydataSub[,paste0(iEd,"F")] <- NA; mydataSub[,iEd] <- NA   }
      }
      for(iName in c("designation","mother","father")){
        mydataSub[,iName] <- as.factor(mydataSub[,iName])
      }
      # remove outliers
      if(!is.null(nrow(cleaning))){ # if there's outliers
        cleaningSub <- cleaning[which(cleaning$trait %in% iTrait),]
        out <- which(mydataSub$rowindex %in% cleaningSub$row )
      }else{out <- numeric()}
      # check the genetic units
      nLevelsGenounit <- apply(data.frame(genoUnit),1,function(x){length(table(mydataSub[,x])) }); names(nLevelsGenounit) <- genoUnit
      genoUnitTraitField <- names(nLevelsGenounit)[which(nLevelsGenounit > 1)]
      if(length(genoUnitTraitField)==0){

        warning(paste("There is not enough levels in your specified genetic unit(s):", paste(genoUnit, collapse = ", "), "in environment", iField, ". Ignoring environment." ),call. = FALSE)

      }else{
        nLevelsFixedunit <- apply(data.frame(setdiff(fixedTerm,"1")),1,function(x){length(table(mydataSub[,x])) }); names(nLevelsFixedunit) <- setdiff(fixedTerm,"1")
        badLevelsFixedunit <- names(nLevelsFixedunit)[which(nLevelsFixedunit <= 1)]
        fixedTermTraitField <- unique(c("1",setdiff(fixedTerm,badLevelsFixedunit)))
        # impute fixed effect columns if they are numeric
        toImpute <- unlist(lapply(mydataSub[,setdiff(fixedTermTraitField,"1")], class))
        keepToImpute <- which(toImpute %in% "numeric")
        if(length(keepToImpute) > 0){
          toImpute <- names(toImpute[keepToImpute])
          for(iImpute in toImpute){mydataSub[, iImpute] <- enhancer::imputev(mydataSub[, iImpute])}
        }
        if(length(out) > 0){mydataSub[out,"trait"] <- NA} # set outliers to NA
        # do analysis
        if(!is.na(var(mydataSub[,"trait"],na.rm=TRUE))){ # if there's variance
          if( var(mydataSub[,"trait"], na.rm = TRUE) > 0 ){
            gridCheck <- with(mydataSub,table(row,col))
            if(nrow(gridCheck) > 1){
              # try to fix assuming they jusy have the rows and cols by replicate
              badRecords <- which(gridCheck > 1, arr.ind = TRUE)
              if(nrow(badRecords) > 0){
                mydataSub <- cgiarBase::fixCoords(mydataSub, rowcoord = "row", colcoord = "col", rep="rep")
              }
              # make the check once again in case that didn't help
              gridCheck <- with(mydataSub,table(row,col))
              badRecords <- which(gridCheck > 1, arr.ind = TRUE)
              if(nrow(badRecords) > 0){
                if(verbose){cat("Replicated records assigned to same row and column in the same trial. Ignoring row and column information for this trial. \n")}
                mydataSub$row <- mydataSub$rowF <- NA
                mydataSub$col <- mydataSub$colF <- NA
              }
            }

            # find best experimental design formula
            mde <- cgiarBase::asremlFormula(fixed=as.formula(paste("trait","~ 1")),
                                            random=~ at(environmentF):rowF + at(environmentF):colF + at(environmentF):trialF + at(environmentF):repF + at(environmentF):iBlockF,
                                            rcov=~at(environmentF):id(rowF):id(colF),
                                            dat=droplevels(mydataSub[which(!is.na(mydataSub[,"trait"])),]),

                                            minRandomLevels=list(rowF= 3, colF=3, trialF=2,repF=2, iBlockF=4),
                                            minResidualLevels=list(rowF=5, colF=5),

                                            exchangeRandomEffects=list(rowF="colF", colF="rowF"),

                                            exchangeResidualEffects=list(rowF="colF", colF="rowF"),

                                            customRandomLevels=NULL, customResidualLevels=NULL,

                                            xCoordinate= "rowF",yCoordinate ="colF",
                                            doubleConstraintRandom=c("rowF","colF"), verbose=verbose)

            factorsFitted <- unlist(lapply(mde$used$environmentF,length))
            factorsFittedGreater <- which(factorsFitted > 0)

            newRandom <- ifelse(length(factorsFittedGreater) > 0, TRUE, FALSE  )
            if(newRandom){newRandom <- names(factorsFitted)[factorsFittedGreater]}else{newRandom<-NULL}
            #
            if((length(mde$used$environmentF$rowF) == 0) & (length(mde$used$environmentF$colF) == 0)){
              newSpline <- NULL
            }else{
              newSpline = as.formula(paste("~spl2D(x1 = row, x2 = col, nseg = c(",min(c(round(nrow(gridCheck)/2), 10)),",", min(c(round(ncol(gridCheck)/2), 10)),") )"))
            }

            for(iGenoUnit in genoUnitTraitField){ # iGenoUnit <- genoUnitTraitField[1]
              myGeneticUnit <-  iGenoUnit
              randomTermForRanModel <- c(newRandom,genoUnit)
              fixedTermForRanModel <- setdiff(fixedTermTraitField,randomTermForRanModel)
              fixedFormulaForRanModel <- paste("trait ~",paste(fixedTermForRanModel, collapse = " + "))
              #
              ranran <- paste(c(myGeneticUnit, unique(intersect(randomTermForRanModel,newRandom)) ), collapse = " + ")
              randomFormulaForRanModel <- paste("~",ranran)
              # at least one condition met: replicated random terms, or replicated fixed terms
              if( (length(factorsFittedGreater) > 0) | (median(table(mydataSub[,iGenoUnit]), na.rm=TRUE) > 1.5) ){

                mixRandom <- try( # first model with genotypes as random
                  LMMsolver::LMMsolve(fixed =as.formula(fixedFormulaForRanModel),
                                      random = as.formula(randomFormulaForRanModel),
                                      spline = newSpline, #trace = TRUE,
                                      family = eval(parse(text = traitFamily[iTrait])),
                                      data = droplevels(mydataSub[which(!is.na(mydataSub[,iGenoUnit])),]), maxit = maxit),
                  silent = TRUE
                );  # mixRandom$VarDf
                # if random model run only keep variance components that were greater than zero and fit again
                # with genotypes as fixed
                if(!inherits(mixRandom,"try-error") ){ # if random model runs well try the fixed effect model  # & ((length(factorsFittedGreater) > 0) )
                  ## save residuals
                  provMydataSub <- droplevels(mydataSub[which(!is.na(mydataSub[,iGenoUnit])),])
                  whereResidualGoes <- provMydataSub[which(!is.na(provMydataSub$trait)),"rowindex"]

                  columnsToAdd <- unique(c(columnsToAdd, paste(iTrait,"residual",sep="-")))
                  mydata[whereResidualGoes,paste(iTrait,"residual",sep="-")] <- mixRandom$residuals[,1]
                  sm <- summary(mixRandom, which = "variances")
                  newRanran <- setdiff( (sm[,1])[which(sm[,2] >0.05)] , c("residual",genoUnitTraitField,"s(row, col)"))
                  removedTerms <- setdiff( setdiff(sm[,1],c("residual",genoUnitTraitField) ) , newRanran ) # to inform users that these were removed

                  ranran <- paste("~",paste(c(newRanran), collapse=" + "))
                  if(ranran=="~ "){randomFormulaForFixedModel=NULL}else{randomFormulaForFixedModel <- as.formula(ranran)}
                  rownames(sm) <- NULL

                  otherFixed <- setdiff(fixedTermTraitField,genoUnitTraitField)
                  fixedFormulaForFixedModel <- paste("trait ~",paste(c(iGenoUnit,otherFixed), collapse = " + "))
                  if(verbose){
                    print(sm)
                    cat(paste(fixedFormulaForFixedModel, "\n"))
                    cat(paste(randomFormulaForFixedModel,"\n"))
                  }
                  mixFixed <- try( # second model with genotype as fixed
                    LMMsolver::LMMsolve(fixed =as.formula(fixedFormulaForFixedModel),
                                        random = randomFormulaForFixedModel,
                                        spline = newSpline, #trace = TRUE,
                                        family = eval(parse(text = traitFamily[iTrait])),
                                        data = droplevels(mydataSub[which(!is.na(mydataSub[,iGenoUnit])),]), maxit = maxit),
                    silent = TRUE
                  ) # mixFixed$VarDf

                  ##################################################################
                  ## if model fails try the simplest model, no random and no spatial
                  if(inherits(mixFixed,"try-error") ){
                    mixFixed <- try( # urgency model, genotypes as fixed but no spatial at all
                      LMMsolver::LMMsolve(fixed =as.formula(fixedFormulaForFixedModel),
                                          family = eval(parse(text = traitFamily[iTrait])),
                                          data = droplevels(mydataSub[which(!is.na(mydataSub[,iGenoUnit])),]), maxit = maxit),
                      silent = TRUE
                    )
                    if( inherits(mixFixed,"try-error") ){
                      if(verbose){cat(paste("Fixed effects models failed. Returning deregressed BLUPs \n"))}
                      currentModeling <- data.frame(module="sta", analysisId=staAnalysisId,trait=iTrait,environment=iField,
                                                    parameter=c("fixedFormula","randomFormula","spatialFormula","family","designationEffectType"),
                                                    value=c( fixedFormulaForRanModel,randomFormulaForRanModel,
                                                             as.character(newSpline)[2],traitFamily[iTrait],"BLUP"))
                      phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[,colnames(phenoDTfile$modeling)])
                      predictedValue <-  mixRandom$coefMME[mixRandom$ndxCoefficients[[iGenoUnit]]] +  mixRandom$coefMME[mixRandom$ndxCoefficients$`(Intercept)`]
                      dims <- mixRandom$EDdf
                      start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) + 1 # we add the one when is random
                      
                      Cinv <- as.matrix(solve(mixRandom$C))
                      pev <- Cinv[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                      
                      #Build X
                      idxIntercept <- mixRandom$ndxCoefficients$`(Intercept)`   # single position
                      idxGenos     <- mixRandom$ndxCoefficients[[iGenoUnit]]    # one per genotype
                      
                      p <- length(mixRandom$coefMME)                 # total β + u length
                      Xpred <- matrix(0, nrow = length(idxGenos), ncol = p)
                      Xpred[ ,idxIntercept] <- 1                                # intercept column
                      for (j in seq_along(idxGenos)) Xpred[j, idxGenos[j]] <- 1  # the right uᵢ
                      
                      stdError <- stdErrorRandom <- sqrt(diag(Xpred %*% Cinv %*% t(Xpred))) # use XC⁻¹X'
                      
                      badSEs <- which( stdError < (sd(predictedValue, na.rm = TRUE)/100) )
                      if(length(badSEs) > 0){stdError[badSEs] <- sd(predictedValue, na.rm = TRUE)}

                      designation <- gsub(paste0(iGenoUnit,"_"),"", names(mixRandom$ndxCoefficients[[iGenoUnit]]))
                      pp <- data.frame(designation,predictedValue,stdError)
                      pp$trait <- iTrait
                      pp$environmentF <- iField
                      pp$entryType <- apply(data.frame(pp$designation),1,function(x){found <-which(mydataSub$designation %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "##"),"unlabeled"); return(x2)})
                      pp$effectType <- "designation"
                      if(iGenoUnit != "designation"){pp$entryType <- paste(iGenoUnit, pp$entryType, sep = "##" )}
                      ## heritabilities
                      ss = mixRandom$VarDf#summary(mixRandom, which = "variances")
                      rownames(ss) <- ss$VarComp
                      vg <- ss[iGenoUnit,2]; vr <- ss["residual",2]
                      cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
                      cv2 <- (sqrt(vr)/mean(as.numeric(mydataSub[,iTrait]),na.rm=TRUE))*100
                      ## reliability
                      A <- diag(nrow(pev))
                      G <- A*vg # G matrix
                      R2 = (G - pev)/G
                      pp$reliability <- diag(R2)
                      pp$predictedValue <- pp$predictedValue/pp$reliability
                      badRels <- which(pp$reliability > 1); if(length(badRels) > 0){pp$reliability[badRels] <- 0.9999}
                      badRels2 <- which(pp$reliability < 0); if(length(badRels2) > 0){pp$reliability[badRels2] <- 0}
                      predictionsList[[counter]] <- pp
                      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                                   data.frame(module="sta",analysisId=staAnalysisId, trait=iTrait, environment=iField,
                                                              parameter=c(paste(c("plotH2","CV", "r2",paste0("V_",as.character(ss$VarComp)),"mean"), iGenoUnit, sep="_"),"CV_environment"),
                                                              method= c(paste( c("vg/(vg+ve)","sd/mu","(G-PEV)/G",rep("REML", nrow(ss)),"sum/n"), iGenoUnit, sep = "-" ),"sqrt(MSE)/GM"),
                                                              value=c(vg/(vg+vr), cv, mean(pp$reliability), ss$Variance, mean(pp$predictedValue,na.rm=TRUE),cv2),
                                                              stdError=c(NA,NA,sd(pp$reliability, na.rm = TRUE)/sqrt(length(pp$reliability)),rep(NA,nrow(ss)), NA, NA)
                                                   )
                      )
                      counter=counter+1
                    }else{ # something still went wrong, return means
                      if(verbose){cat(paste("No design to fit or singularities encountered in the random model, aggregating and assuming h2 = 0 \n"))}
                      pp <- aggregate(as.formula(paste("trait ~", iGenoUnit)), FUN=mean, data=mydataSub)
                      colnames(pp)[1:2] <- c("designation","predictedValue")
                      pp$stdError <- sd(pp$predictedValue) #1
                      pp$reliability <- 1e-6
                      pp$trait <- iTrait
                      pp$environmentF <- iField
                      pp$entryType <- apply(data.frame(pp$designation),1,function(x){found <-which(mydataSub$designation %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unlabeled"); return(x2)})
                      pp$effectType <- "designation"
                      if(iGenoUnit != "designation"){pp$entryType <- paste(iGenoUnit, pp$entryType, sep = "##" )}
                      predictionsList[[counter]] <- pp;
                      cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
                      phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                                   data.frame(module="sta",analysisId=staAnalysisId, trait=iTrait, environment=iField,
                                                              parameter=c(paste( c("plotH2","CV", "r2","V_designation","V_residual","mean"), iGenoUnit, sep="_"),"CV_environment"),
                                                              method=c(paste( c("vg/(vg+ve)","sd/mu","(G-PEV)/G","REML","REML","sum/n"), iGenoUnit, sep = "-"),"sqrt(MSE)/GM") ,
                                                              value=c(0, cv, 0, 0, 0, mean(pp$predictedValue,na.rm=TRUE),0 ),
                                                              stdError=c(NA,NA,NA,NA, NA, NA, NA)
                                                   )
                      )
                      currentModeling <- data.frame(module="sta", analysisId=staAnalysisId,trait=iTrait,environment=iField,
                                                    parameter=c("fixedFormula","randomFormula","spatialFormula","family","designationEffectType"),
                                                    value=c( ifelse(returnFixedGeno, fixedFormulaForFixedModel, fixedFormulaForRanModel),
                                                             ifelse(returnFixedGeno, as.character(randomFormulaForFixedModel)[2], randomFormulaForRanModel ),
                                                             as.character(newSpline)[2],traitFamily[iTrait],ifelse(returnFixedGeno,"BLUE","BLUP")))
                      phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[,colnames(phenoDTfile$modeling)])
                      counter=counter+1
                    }
                  }else{ # fixed model run well
                    currentModeling <- data.frame(module="sta", analysisId=staAnalysisId,trait=iTrait,environment=iField,
                                                  parameter=c("fixedFormula","randomFormula","randomTermsRemoved","spatialFormula","family","designationEffectType"),
                                                  value=c( ifelse(returnFixedGeno, fixedFormulaForFixedModel, fixedFormulaForRanModel),
                                                           ifelse(returnFixedGeno, as.character(randomFormulaForFixedModel)[2], randomFormulaForRanModel ),
                                                           paste(removedTerms, collapse = "+"),
                                                           as.character(newSpline)[2],traitFamily[iTrait],ifelse(returnFixedGeno,"BLUE","BLUP")))
                    phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[,colnames(phenoDTfile$modeling)])

                    if(returnFixedGeno){ # user wants fixed effect predictions for genotype
                      shouldBeOne <- which(mixFixed$ndxCoefficients[[iGenoUnit]] == 0)
                      if(length(shouldBeOne) > 0){mixFixed$ndxCoefficients[[iGenoUnit]][shouldBeOne] = 1}
                      predictedValue <- mixFixed$coefMME[mixFixed$ndxCoefficients[[iGenoUnit]]] + mixFixed$coefMME[mixFixed$ndxCoefficients$`(Intercept)`]
                      if(length(shouldBeOne) > 0){predictedValue[1] = mixFixed$coefMME[mixFixed$ndxCoefficients$`(Intercept)`]} # adjust the value for first entry
                      dims <- mixFixed$EDdf
                      start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) # we don't add a one because we need the intercept
                      
                      Cinv <- as.matrix(solve(mixFixed$C))
                     
                      #Build X
                      idxIntercept <- mixFixed$ndxCoefficients$`(Intercept)`   # single position
                      idxGenos     <- mixFixed$ndxCoefficients[[iGenoUnit]]    # one per genotype
                      
                      p <- length(mixFixed$coefMME)                 # total β + u length
                      Xpred <- matrix(0, nrow = length(idxGenos), ncol = p)
                      Xpred[ ,idxIntercept] <- 1                                # intercept column
                      for (j in seq_along(idxGenos)) Xpred[j, idxGenos[j]] <- 1  # the right uᵢ
                      
                      stdError <- sqrt(diag(Xpred %*% Cinv %*% t(Xpred))) # use XC⁻¹X'
                      
                      badSEs <- which( stdError < (sd(predictedValue, na.rm = TRUE)/100) )
                      if(length(badSEs) > 0){stdError[badSEs] <- sd(predictedValue, na.rm = TRUE)}
                      # just for reliability calculation
                      dims <- mixRandom$EDdf
                      start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) + 1 # we add the one when is random
                      pev <- as.matrix(solve(mixRandom$C))[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                      # stdErrorRandom <- (sqrt(diag(pevRandom)))
                    }else{ # user wants random effect predictions for genotype
                      predictedValue <- mixRandom$coefMME[mixRandom$ndxCoefficients[[iGenoUnit]]] +  mixRandom$coefMME[mixRandom$ndxCoefficients$`(Intercept)`]
                      dims <- mixRandom$EDdf
                      start <- sum(dims[1:(which(dims$Term == iGenoUnit) - 1),"Model"]) + 1 # we add the one when is random
                      
                      Cinv <- as.matrix(solve(mixRandom$C))
                      pev <- Cinv[start:(start+length(predictedValue)-1),start:(start+length(predictedValue)-1)]
                      
                      #Build X
                      idxIntercept <- mixRandom$ndxCoefficients$`(Intercept)`   # single position
                      idxGenos     <- mixRandom$ndxCoefficients[[iGenoUnit]]    # one per genotype
                      
                      p <- length(mixRandom$coefMME)                 # total β + u length
                      Xpred <- matrix(0, nrow = length(idxGenos), ncol = p)
                      Xpred[ ,idxIntercept] <- 1                                # intercept column
                      for (j in seq_along(idxGenos)) Xpred[j, idxGenos[j]] <- 1  # the right uᵢ
                      
                      stdError <- stdErrorRandom <- sqrt(diag(Xpred %*% Cinv %*% t(Xpred))) # use XC⁻¹X'
                      
                      # move to std deviation if model is wrong and stdError is close to zero
                      badSEs <- which( stdError < (sd(predictedValue, na.rm = TRUE)/100) )
                      if(length(badSEs) > 0){stdError[badSEs] <- sd(predictedValue, na.rm = TRUE)}
                    }
                    designation <- gsub(paste0(iGenoUnit,"_"),"", names(mixRandom$ndxCoefficients[[iGenoUnit]]))
                    pp <- data.frame(designation,predictedValue,stdError)
                    pp$trait <- iTrait
                    pp$environmentF <- iField
                    pp$entryType <- apply(data.frame(pp$designation),1,function(x){found <-which(mydataSub$designation %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unlabeled"); return(x2)})
                    pp$effectType <- "designation"
                    if(iGenoUnit != "designation"){pp$entryType <- paste(iGenoUnit, pp$entryType, sep = "##" )}
                    ## heritabilities
                    ss = mixRandom$VarDf
                    rownames(ss) <- ss$VarComp
                    vg <- ss[iGenoUnit,2]; vr <- ss["residual",2]
                    cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
                    cv2 <- (sqrt(vr)/mean(as.numeric(mydataSub[,iTrait]),na.rm=TRUE))*100
                    ## reliability
                    A <- diag(nrow(pev))
                    G <- A*vg # G matrix
                    R2 = (G - pev)/G
                    pp$reliability <- diag(R2)
                    badRels <- which(pp$reliability > 1); if(length(badRels) > 0){pp$reliability[badRels] <- 0.9999}
                    badRels2 <- which(pp$reliability < 0); if(length(badRels2) > 0){pp$reliability[badRels2] <- 0}
                    phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                                 data.frame(module="sta",analysisId=staAnalysisId, trait=iTrait, environment=iField,
                                                            parameter=c(paste( c("plotH2","CV", "r2",paste0("V_",as.character(ss$VarComp)),"mean"), iGenoUnit, sep="_"),"CV_environment"),
                                                            method=c(paste( c("vg/(vg+ve)","sd/mu","(G-PEV)/G",rep("REML",nrow(ss)),"sum/n"), iGenoUnit, sep = "-" ),"sqrt(MSE)/GM"),
                                                            value=c(vg/(vg+vr), cv, mean(pp$reliability), ss$Variance, mean(pp$predictedValue,na.rm=TRUE), cv2),
                                                            stdError=c(NA,NA,sd(pp$reliability, na.rm = TRUE)/sqrt(length(pp$reliability)),rep(NA, nrow(ss)), NA, NA)
                                                 )
                    )
                    predictionsList[[counter]] <- pp
                    counter=counter+1

                  }
                  ##################################################################
                }else{ # if there was singularities we just take means and assigna h2 of zero
                  if(verbose){cat(paste("No design to fit or singularities encountered in the random model, aggregating and assuming h2 = 0 \n"))}
                  pp <- aggregate(as.formula(paste("trait ~", iGenoUnit)), FUN=mean, data=mydataSub)
                  colnames(pp)[1:2] <- c("designation","predictedValue")
                  pp$stdError <- sd(pp$predictedValue) #1
                  pp$reliability <- 1e-6
                  pp$trait <- iTrait
                  pp$environmentF <- iField
                  pp$entryType <- apply(data.frame(pp$designation),1,function(x){found <-which(mydataSub$designation %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unlabeled"); return(x2)})
                  pp$effectType <- "designation"
                  if(iGenoUnit != "designation"){pp$entryType <- paste(iGenoUnit, pp$entryType, sep = "##" )}
                  predictionsList[[counter]] <- pp;
                  cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
                  phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                               data.frame(module="sta",analysisId=staAnalysisId, trait=iTrait, environment=iField,
                                                          parameter=c(paste( c("plotH2","CV", "r2","V_designation","V_residual","mean"), iGenoUnit, sep="_"),"CV_environment"),
                                                          method=c(paste( c("vg/(vg+ve)","sd/mu","(G-PEV)/G","REML","REML","sum/n"), iGenoUnit, sep = "-"),"sqrt(MSE)/GM"),
                                                          value=c(0, cv, 0, 0, 0, mean(pp$predictedValue,na.rm=TRUE),0 ),
                                                          stdError=c(NA,NA,NA,NA, NA, NA, NA)
                                               )
                  )
                  currentModeling <- data.frame(module="sta", analysisId=staAnalysisId,trait=iTrait,environment=iField,
                                                parameter=c("fixedFormula","randomFormula","spatialFormula","family","designationEffectType"),
                                                value=c( ifelse(returnFixedGeno, fixedFormulaForFixedModel, fixedFormulaForRanModel),
                                                         ifelse(returnFixedGeno, as.character(randomFormulaForFixedModel)[2], randomFormulaForRanModel ),
                                                         as.character(newSpline)[2],traitFamily[iTrait],ifelse(returnFixedGeno,"BLUE","BLUP")))
                  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[,colnames(phenoDTfile$modeling)])
                  counter=counter+1
                } # end of is mixed model run well

              }else{

                if(verbose){cat(paste("No design to fit, aggregating and assuming h2 = 0 \n"))}
                pp <- aggregate(as.formula(paste("trait ~", iGenoUnit)), FUN=mean, data=mydataSub)
                colnames(pp)[1:2] <- c("designation","predictedValue")
                pp$stdError <- sd(pp$predictedValue)  # 1
                pp$reliability <- 1e-6
                pp$trait <- iTrait
                pp$environmentF <- iField
                pp$entryType <- apply(data.frame(pp$designation),1,function(x){found <-which(mydataSub$designation %in% x); x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unlabeled"); return(x2)})
                pp$effectType <- "designation"
                if(iGenoUnit != "designation"){pp$entryType <- paste(iGenoUnit, pp$entryType, sep = "##" )}
                predictionsList[[counter]] <- pp;
                cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100
                phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                             data.frame(module="sta",analysisId=staAnalysisId, trait=iTrait, environment=iField,
                                                        parameter= c(paste( c("plotH2","CV", "r2","V_designation","V_residual", "mean"), iGenoUnit, sep="_"),"CV_environment"),
                                                        method=c(paste( c("vg/(vg+ve)","sd/mu","(G-PEV)/G","REML","REML","sum/n"), sep = "-"),"sqrt(MSE)/GM") ,
                                                        value=c(0, cv, 0, 0, 0, mean(pp$predictedValue,na.rm=TRUE),0 ), stdError=c(NA,NA,NA,NA,NA,NA,NA)
                                             )
                )
                currentModeling <- data.frame(module="sta", analysisId=staAnalysisId,trait=iTrait,environment=iField, parameter=c("fixedFormula","randomFormula","spatialFormula","family","designationEffectType"), value=c("None","None","None","None","mean"))
                phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[,colnames(phenoDTfile$modeling)])
                counter=counter+1

              }

            }


          }else{

            for(iGenoUnit in genoUnitTraitField){ # iGenoUnit = genoUnitTraitField[1]

              if(verbose){
                cat(paste("No design to fit, aggregating for predicted values, std. errors assumed equal to std. deviation of the trial. In addition assuming h2 = 0 for the trial \n"))
              }
              pp <- aggregate(as.formula(paste("trait ~", iGenoUnit)), FUN=mean, data=mydataSub)
              colnames(pp)[1:2] <- c("designation","predictedValue")
              pp$stdError <- sd(pp$predictedValue)  # 1
              pp$trait <- iTrait
              pp$environmentF <- iField
              pp$entryType <- apply(data.frame(pp$designation),1,function(x){
                found <-which(mydataSub$designation %in% x);
                x2 <- ifelse(length(found) > 0, paste(sort(unique(toupper(trimws(mydataSub[found,"entryType"])))), collapse = "#"),"unlabeled");
                return(x2)
              }
              )
              pp$effectType <- "designation"
              if(iGenoUnit != "designation"){pp$entryType <- paste(iGenoUnit, pp$entryType, sep = "##" )}
              pp$reliability <- 1e-6
              predictionsList[[counter]] <- pp;
              cv <- (sd(pp$predictedValue,na.rm=TRUE)/mean(pp$predictedValue,na.rm=TRUE))*100

              phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                           data.frame(module="sta",analysisId=staAnalysisId, trait=iTrait, environment=iField,
                                                      parameter= c(paste( c("plotH2","CV", "r2","V_designation","V_residual","mean"), iGenoUnit, sep="_"),"CV_environment") ,
                                                      method=c(paste( c("vg/(vg+ve)","sd/mu","(G-PEV)/G","REML","REML","sum/n"), iGenoUnit, sep = "-"),"sqrt(MSE)/GM") ,
                                                      value=c(0, cv, 0,0,0, mean(pp$predictedValue,na.rm=TRUE),0), stdError=c(NA,NA,NA,NA, NA, NA,NA)
                                           )
              )

              currentModeling <- data.frame(module="sta", analysisId=staAnalysisId,trait=iTrait,environment=iField, parameter=c("fixedFormula","randomFormula","spatialFormula","family","designationEffectType"), value=c("None","None","None","None","mean"))
              phenoDTfile$modeling <- rbind(phenoDTfile$modeling, currentModeling[,colnames(phenoDTfile$modeling)])
              counter=counter+1

            }
          }
        }

      }



    }
  }
  predictionsBind <- do.call(rbind, predictionsList)  

  if(nrow(predictionsBind) == 0){
    stop( "No predictions to work with.",call. = FALSE)
  }

  predictionsBind$analysisId <- staAnalysisId
  predictionsBind$module <- "sta"
  colnames(predictionsBind) <- cgiarBase::replaceValues(Source=colnames(predictionsBind), Search=c("designation","environmentF"), Replace=c("designation","environment"))

  ##########################################
  ## add timePoint of origin and stage and designation code
  entries <- unique(predictionsBind[,"designation"])
  baseOrigin <- do.call(rbind, apply(data.frame(entries),1,function(x){
    out1 <- (sort(mydata[which(mydata$designation %in% x),"gid"], decreasing = FALSE))[1]
    out2 <- (sort(mydata[which(mydata$designation %in% x),"mother"], decreasing = FALSE))[1]
    out3 <- (sort(mydata[which(mydata$designation %in% x),"father"], decreasing = FALSE))[1]
    out4 <- paste(unique(sort(mydata[which(mydata$designation %in% x),"pipeline"], decreasing = FALSE)),collapse=", ")
    y <- data.frame(designation=x,gid=out1,mother=out2,father=out3,pipeline=out4)
    return(y)
  }))
  predictionsBind <- merge(predictionsBind,baseOrigin, by="designation", all.x=TRUE)

  ##########################################
  ## update data tables
  '%!in%' <- function(x,y)!('%in%'(x,y))

  if( all( c("predictions" %in% names(phenoDTfile), "effectType" %!in% colnames(predictionsBind) ) ) ){
    phenoDTfile$predictions$effectType <- NA
  }
  phenoDTfile$predictions <- rbind(phenoDTfile$predictions, predictionsBind[,colnames(phenoDTfile$predictions)] )
  phenoDTfile$predictions$stdError = as.numeric(phenoDTfile$predictions$stdError)
  newStatus <- data.frame(module="sta", analysisId=staAnalysisId, analysisIdName=NA)
  phenoDTfile$status <- rbind( phenoDTfile$status, newStatus[,colnames(phenoDTfile$status)])
  ### change column names back for mapping
  colnames(mydata) <- cgiarBase::replaceValues(colnames(mydata), Replace = paramsPheno$value,  Search= paramsPheno$parameter )
  ##
  #phenoDTfile$data$pheno <- cbind(phenoDTfile$data$pheno, mydata[,columnsToAdd]) #mydata[,unique(c(originalColumns,columnsToAdd))]#mydata[,-which(colnames(mydata) %in% c("mother","father") )]
  phenoDTfile$data$pheno[, columnsToAdd] <- mydata[, columnsToAdd, drop = FALSE]
  ## add which analysisId was used as input
  modeling <- data.frame(module="sta",  analysisId=staAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling)
  return(phenoDTfile)
}
