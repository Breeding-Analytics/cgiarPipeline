rggPiepho <- function(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, # per trait
    yearsToUse=NULL,
    environmentToUse=NULL,
    entryTypeToUse=NULL,
    effectTypeToUse=NULL,
    verbose=TRUE,
    traitFamily=NULL,
    sampleN=50, # max number of individuals per environment to sample
    bootstrappingN=10,
    forceRules=TRUE
){
  ## THIS FUNCTION CALCULATES THE REALIZED GENETIC GAIN FOR SOME TRAITS
  ## IS USED IN THE BANAL APP UNDER THE METRICS MODULES
  fixedTerm="yearOfOrigin"
  gTerm <- "designation"
  rggAnalysisId <- as.numeric(Sys.time())
  if(is.null(traitFamily)){traitFamily <- rep("quasi(link = 'identity', variance = 'constant')", length(trait))}
  if(length(traitFamily) != length(trait)){stop("Trait distributions should have the same length than traits to be analyzed.", call. = FALSE)}
  names(traitFamily) <- trait
  if(is.null(phenoDTfile)){stop("Please provide the name of the analysis to locate the predictions", call. = FALSE)}
  if(is.null(analysisId)){stop("Please provide the ID of the analysis to use as input", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be analyzed", call. = FALSE)}
  ############################
  # loading the dataset
  '%!in%' <- function(x,y)!('%in%'(x,y))
  if(!is.null(phenoDTfile$predictions)){
    if("effectType" %!in% colnames(phenoDTfile$predictions) ){
      phenoDTfile$predictions$effectType <- "general"
    }
  }
  mydata <- phenoDTfile$predictions
  mydata <- mydata[which(mydata$analysisId %in% analysisId),]
  if(nrow(mydata)==0){stop("No match for this analysisId. Please correct.", call. = FALSE)}
  ## add male, female and yearOfOrigin columns
  myPed <- phenoDTfile$data$pedigree
  paramsPed <- phenoDTfile$metadata$pedigree
  colnames(myPed) <- cgiarBase::replaceValues(colnames(myPed), Search = paramsPed$value, Replace = paramsPed$parameter )
  myPed <- unique(myPed[,c(gTerm,fixedTerm)])
  if(is.null(myPed) || (nrow(myPed) == 0 ) ){stop("yearOfOrigin column was not matched in your original file. Please correct.", call. = FALSE)}
  mydata <- merge(mydata, myPed[,c(gTerm,fixedTerm)], by=gTerm, all.x=TRUE )
  mydata <- mydata[which(!is.na(mydata$yearOfOrigin)),]
  # add the other available columns to the dataset
  metaPheno <- phenoDTfile$metadata$pheno[which(phenoDTfile$metadata$pheno$parameter %in% c("environment","year","season","country","location","trial")),]

  otherMetaCols <- unique(phenoDTfile$data$pheno[,metaPheno$value,drop=FALSE])
  colnames(otherMetaCols) <- cgiarBase::replaceValues(Source = colnames(otherMetaCols), Search = metaPheno$value, Replace = metaPheno$parameter )
  otherMetaCols <- otherMetaCols[which(!duplicated(otherMetaCols[,"environment"])),,drop=FALSE] # we do this in case the users didn't define the environment properly
  mydata <- merge(mydata, otherMetaCols, by="environment", all.x = TRUE)
  ##
  if(!is.null(yearsToUse)){ # reduce the dataset
    yearsToUse <- as.numeric(as.character(yearsToUse))
    mydata <- mydata[which(mydata$yearOfOrigin %in% yearsToUse),]
  }
  if(!is.null(environmentToUse)){ # reduce the dataset
    environmentToUse <- as.character(environmentToUse)
    mydata <- mydata[which(mydata$environment %in% environmentToUse),]
  }
  if(!is.null(entryTypeToUse)){ # reduce the dataset
    entryTypeToUse <- as.character(entryTypeToUse)
    mydata <- mydata[which(mydata$entryType %in% entryTypeToUse),]
  }
  if(!is.null(effectTypeToUse)){ # reduce the dataset
    effectTypeToUse <- as.character(effectTypeToUse)
    mydata <- mydata[which(mydata$effectType %in% effectTypeToUse),]
  }
  if(nrow(mydata) == 0){stop("No data to work with with the specified parameters. You may want to check the yearsToUse parameter. Maybe you have not mapped the 'yearOfOrigin' column in the Data Retrieval tad under the 'Pedigree' section.",call. = FALSE)}

  if(forceRules){ # if we enforce ABI rules for minimum years of data
    if(length(unique(na.omit(mydata[,fixedTerm]))) < 5){
      warning("Less than 5 years of data have been detected. Make sure you have mapped the 'yearOfOrigin' column in the Data Retrieval tad under the 'Pedigree' section. If yes, interpret results with caution ", call. = FALSE)
      less_than_five = TRUE
    }else{less_than_five = FALSE}
  }else{
    if(length(unique(na.omit(mydata[,fixedTerm]))) <= 1){stop("Only one year of data. Realized genetic gain analysis cannot proceed.Maybe you have not mapped the 'yearOfOrigin' column in the Data Retrieval tad under the 'Pedigree' section. ", call. = FALSE)}
  }
  # define wether we should deregress or not
  modelingInput <- phenoDTfile$modeling
  modelingInput <- modelingInput[which(modelingInput$analysisId == analysisId),]
  designationEffectType <- modelingInput[which(modelingInput$parameter == "randomFormula"),"value"]
  if(length(grep("designation", designationEffectType)) > 0){
    deregress=TRUE
  }else{ # BLUE
    deregress=FALSE
  }
  # if(unique(modelingInput$module) == "sta"){
  #   designationEffectType <- modelingInput[which(modelingInput$parameter == "randomFormula"),"value"]
  #   if(length(grep("designation", designationEffectType)) > 0){
  #     deregress=TRUE
  #   }else{ # BLUE
  #     deregress=FALSE
  #   }
  # }else{ # mtaLmms
  #   designationEffectType <- modelingInput[which(modelingInput$parameter == "randomFormula"),"value"]
  #   if(length(grep("grp\\(designation\\)", designationEffectType)) > 0){
  #     deregress=TRUE
  #   }else{ # BLUE
  #     deregress=FALSE
  #   }
  # }
  
  # remove traits that are not actually present in the dataset
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% unique(mydata$trait)){
      if(verbose){
        cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))
      }
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  trait <- setdiff(trait,traitToRemove)
  if(length(trait)==0){stop("None of the traits specified are available. Please double check", call. = FALSE)}
  ############################
  ## gg analysis
  counter=1
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){
      cat(paste("Analyzing trait", iTrait,"\n"))
    }
    # subset data
    mydataSub <- droplevels(mydata[which(mydata$trait == iTrait),])
    mydataSub$designation <- as.factor(mydataSub$designation)
    mydataSub$predictedValue.d <- mydataSub$predictedValue/mydataSub$rel
    mydataSub$yearOfOrigin <- as.numeric(mydataSub$yearOfOrigin)
    for(iMetaCol in colnames(otherMetaCols)){
      mydataSub[,iMetaCol] <- as.factor(mydataSub[,iMetaCol])
    }
    val <- stdError <- list()
    for(iBoot in 1:bootstrappingN){ # iBoot=1
      ll <- split(mydataSub, mydataSub$environment)
      ll <- lapply(ll, function(x){y <- x[sample(1:nrow(x), min(c(nrow(x), sampleN)) ),]; return(y)})
      mydataSub2 <- do.call(rbind, ll)

      # do analysis
      if(!is.na(var(mydataSub2[,"predictedValue"],na.rm=TRUE))){ # if there's variance
        if( var(mydataSub2[,"predictedValue"], na.rm = TRUE) > 0 ){
          checks <- as.character(unique(mydataSub2[grep("check",mydataSub2[,"entryType"], ignore.case = TRUE),gTerm]))
          ranran <- "~NULL"
          if(deregress){
            mydataSub2$predictedValue <- mydataSub2$predictedValue.d
            if(verbose){
              print("Deregressing predicted values using the reliability. Assuming you are providing BLUPs.")
            }
          }else{
            if(verbose){
              print("Using predicted values directly. Assuming you are providing BLUEs.")
            }
          }
          # mydataSub2[,"designationLocationYear"] <- paste0(mydataSub2[,"designation"],mydataSub2[,"location"], mydataSub2[,"year"])
          # mydataSub2[,"designationLocation"] <- paste0(mydataSub2[,"designation"],mydataSub2[,"location"])
          mydataSub2[,"designationEnvironment"] <- paste0(mydataSub2[,"designation"],mydataSub2[,"environment"])
          # mydataSub2[,"designationYear"] <- paste0(mydataSub2[,"designation"],mydataSub2[,"year"])
          # mydataSub2[,"locationYear"] <- paste0(mydataSub2[,"location"],mydataSub2[,"year"])
          if(length(which(metaPheno$parameter == "year")) > 0){
            if(length(table(mydataSub2$year)) > 1){
              fix <- paste("predictedValue ~ year + yearOfOrigin + environment")
            }else{
              fix <- paste("predictedValue ~ yearOfOrigin + environment")
            }
          }else{
            fix <- paste("predictedValue ~ yearOfOrigin + environment")
          }

          random <- "~designation + designationEnvironment"
          ranres <- "~units"#"~dsum(~units | environment)"
          mydataSub2=mydataSub2[with(mydataSub2, order(environment)), ]
          mydataSub2$w <- 1/(mydataSub2$stdError)
          # plot(mydataSub2$predictedValue ~ mydataSub2$yearOfOrigin)
          mix <- try(
            LMMsolver::LMMsolve(fixed =as.formula(fix),
                                random = as.formula(random),
                                residual=NULL, #as.formula(ranres),
                                weights = "w",
                                ginverse = NULL,
                                # group = groupTrait,
                                family = eval(parse(text = traitFamily[iTrait])),
                                data = mydataSub2, maxit = 30),
            silent = TRUE
          )
          mixCoeff <- coef(mix, se=TRUE)
          b0 <- mixCoeff$`(Intercept)`[1,"value"]
          b1 <- mixCoeff$yearOfOrigin[1,"value"] #*ifelse(deregress,deregressWeight,1)
          seb1 <- mixCoeff$yearOfOrigin[1,"se"]
          seb0 <- mixCoeff$`(Intercept)`[1,"se"]
          baseline <- b0 + ( b1*min(as.numeric(mydataSub2[which(mydataSub2$trait == iTrait),fixedTerm]) , na.rm=TRUE ))
          baselineAverageYear <- b0 + ( b1*min(as.numeric(mydataSub2[which(mydataSub2$trait == iTrait),fixedTerm]) , na.rm=TRUE ))
          b1Perc <- round(( b1 /baseline) * 100,3)
          b1PercAverageYear <- round(( b1 /baselineAverageYear) * 100,3)
          b1PercSe <- round((seb1/baseline) * 100,3)
          b1PercSeAverageYear <- round((seb1/baselineAverageYear) * 100,3)
          if(length(which(metaPheno$parameter == "year")) > 0){
            ngt <- baseline + mixCoeff$year$value
            ngtSe <- mixCoeff$year$se
          }else{
            ngt <- NA
            ngtSe <- NA
          }
          r2 <- NA # sm$r.squared
          pv <- NA# 1 - pf(sm$fstatistic[1], df1=sm$fstatistic[2], df2=sm$fstatistic[3])

          gg.y1<- sort(unique(mydataSub2[,fixedTerm]), decreasing = FALSE)[1]
          gg.yn <- sort(unique(mydataSub2[,fixedTerm]), decreasing = TRUE)[1]
          ntrial <- phenoDTfile$predictions # number of trials
          ntrial <- ntrial[which(ntrial$analysisId ==analysisId),]
          ntrial <- ntrial[which(ntrial$trait ==iTrait),]
          ntrial <- unique(ntrial$environment)
          ntrial <- length(ntrial[which(ntrial %in% environmentToUse)] )
          val[[iBoot]] <- c(b1,b0, b1Perc, b1PercAverageYear, r2, pv, ntrial,gg.y1,gg.yn, ngt  )
          stdError[[iBoot]] <- c(seb1,seb0,b1PercSe,b1PercSeAverageYear,0,0,0,0,0, ngtSe)
          counter=counter+1
        }
      }
    }
    # save the modeling table
    currentModeling <- data.frame(module="rgg", analysisId=rggAnalysisId,trait=iTrait, environment="across",
                                  parameter=c("deregression","lessThanFiveYears","fixedFormula","randomFormula","residualFormula","family"),
                                  value=c(deregress, less_than_five, fix, random,"~units", traitFamily[iTrait]))
    phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
    # save parameters
    val <- apply(do.call(rbind, val),2, function(x){median(x, na.rm=TRUE)})
    stdError <- apply(do.call(rbind, stdError),2, function(x){median(x, na.rm=TRUE)})
    phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                 data.frame(module="rgg",analysisId=rggAnalysisId, trait=iTrait, environment=c(rep("across",9), mixCoeff$year$coef),
                                            parameter=c("ggSlope","ggInter", "gg%(first.year)","gg%(average.year)","r2","pVal","nTrial","initialYear","lastYear",rep("nonGeneticTrend",length(mixCoeff$year$coef) )),
                                            method=ifelse(deregress,"piephoDeregress","piepho"),
                                            value=val, stdError=stdError
                                 )
    )
    myPreds <- data.frame(module = "rgg", analysisId = rggAnalysisId, pipeline = NA, trait = iTrait,
                          gid = NA, designation = gsub("designation_","",mixCoeff$designation$coef),
                          mother = NA, father = NA, entryType = NA, effectType=NA,
                          environment = "across", predictedValue = mixCoeff$designation$value + baseline,
                          stdError = mixCoeff$designation$se, reliability = NA)
    phenoDTfile$predictions <- rbind(phenoDTfile$predictions, myPreds[,colnames(phenoDTfile$predictions)])

  }
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  newStatus <- data.frame(module="rgg", analysisId=rggAnalysisId, analysisIdName=NA)
  phenoDTfile$status <- rbind( phenoDTfile$status,newStatus[,colnames(phenoDTfile$status)] )
  ## add which data was used as input
  modeling <- data.frame(module="rgg",  analysisId=rggAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)#
}
