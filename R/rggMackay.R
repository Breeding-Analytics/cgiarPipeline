rggMackay <- function(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, # per trait
    deregressWeight=1,
    partition=FALSE,
    environmentToUse=NULL,
    yearsToUse=NULL,
    entryTypeToUse=NULL,
    effectTypeToUse=NULL,
    verbose=TRUE,
    forceRules=TRUE,
    propTopIndsPerYear=1 # by default we pick all inds per year
){
  ## THIS FUNCTION CALCULATES THE REALIZED GENETIC GAIN FOR SOME TRAITS
  ## IS USED IN THE BANAL APP UNDER THE METRICS MODULES
  fixedTerm="yearOfOrigin"
  gTerm <- "designation"
  rggAnalysisId <- as.numeric(Sys.time())
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
  # add male, female and yearOfOrigin columns
  myPed <- phenoDTfile$data$pedigree
  paramsPed <- phenoDTfile$metadata$pedigree
  colnames(myPed) <- cgiarBase::replaceValues(colnames(myPed), Search = paramsPed$value, Replace = paramsPed$parameter )
  myPed <- unique(myPed[,c(gTerm,fixedTerm)])
  if(is.null(myPed) || (nrow(myPed) == 0 ) ){stop("yearOfOrigin column was not matched in your original file. Please correct.", call. = FALSE)}
  mydata <- merge(mydata, myPed[,c(gTerm,fixedTerm)], by=gTerm, all.x=TRUE )
  mydata <- mydata[which(!is.na(mydata$yearOfOrigin)),]
  if(!is.null(yearsToUse)){ # reduce the dataset by years selected
    yearsToUse <- as.numeric(as.character(yearsToUse))
    mydata <- mydata[which(mydata$yearOfOrigin %in% yearsToUse),]
  }
  if(!is.null(environmentToUse)){ # reduce the dataset
    environmentToUse <- as.character(environmentToUse)
    mydata <- mydata[which(mydata$environment %in% environmentToUse),]
  }
  if(!is.null(entryTypeToUse)){ # reduce the dataset by entry types selected
    entryTypeToUse <- as.character(entryTypeToUse)
    mydata <- mydata[which(mydata$entryType %in% entryTypeToUse),]
  }
  if(!is.null(effectTypeToUse)){ # reduce the dataset by entry types selected
    effectTypeToUse <- as.character(effectTypeToUse)
    mydata <- mydata[which(mydata$effectType %in% effectTypeToUse),]
  }
  # reduce dataset by top entries selected
  if(nrow(mydata) == 0){stop("No data to work with with the specified parameters. You may want to check the yearsToUse parameter. Maybe you have not mapped the 'yearOfOrigin' column in the Data Retrieval tad under the 'Pedigree' section.",call. = FALSE)}
  if(forceRules){ 
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
  
  #Get intercept
  preds_this <- phenoDTfile$predictions[phenoDTfile$predictions$analysisId == analysisId, , drop = FALSE]
  
  # Intercept μ by trait
  mu_by_trait <- tapply(
    preds_this$predictedValue[
      preds_this$effectType == "(Intercept)" & preds_this$environment == "(Intercept)"
    ],
    preds_this$trait[
      preds_this$effectType == "(Intercept)" & preds_this$environment == "(Intercept)"
    ],
    mean
  )
  
  #Helper functions to use it to "clean" blups from fixed effects
  split_interaction <- function(s) {
    x <- trimws(unlist(strsplit(s, "[:_]", perl = TRUE)))
    x[nzchar(x)]
  }
  
  get_fixed_main_terms <- function(fixed_formula_chr) {
    if (!length(fixed_formula_chr)) return(character())
    ff <- try(stats::as.formula(fixed_formula_chr[1]), silent = TRUE)
    if (inherits(ff, "try-error")) return(character())
    tt <- attr(stats::terms(ff), "term.labels")   # includes interactions
    # keep ONLY main effects (drop any term with '_')
    setdiff(tt[!grepl(":", tt, fixed = TRUE)], character(0))
  }
  
  get_random_interactions_with_designation <- function(random_formula_chr) {
    if (!length(random_formula_chr)) return(list())
    rf <- random_formula_chr[1]
    grp_inside <- regmatches(rf, gregexpr("grp\\(([^)]*)\\)", rf, perl = TRUE))[[1]]
    if (!length(grp_inside)) return(list())
    insides <- gsub("^grp\\(|\\)$", "", grp_inside)
    raw_terms <- unlist(strsplit(insides, "\\+")); raw_terms <- trimws(raw_terms)
    term_lists <- lapply(raw_terms, split_interaction)
    # keep vectors that contain 'designation' AND at least one other factor
    Filter(function(v) any(v == "designation") && length(v) >= 2, term_lists)
  }
  
  get_fixed_blues_table <- function(trait_name, effect_name) {
    subset(preds_this,
           trait == trait_name &
             effectType == effect_name &
             environment == "(Intercept)",
           select = c("designation","predictedValue"))
  }
  
  choose_levels_for_factor <- function(mydataSub, v, blues_tbl) {
    if (v %in% names(mydataSub)) {
      levs <- unique(as.character(mydataSub[[v]]))
      levs <- intersect(levs, as.character(blues_tbl$designation))
      if (length(levs)) return(levs)
    }
    as.character(blues_tbl$designation)
  }

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
    ## subset to top individuals declared by user
    yearSplit <- split(mydataSub, mydataSub$yearOfOrigin)
    yearSplit <- lapply(yearSplit, function(x){
      y <- x[with(x, order(-predictedValue)), ]
      y <- y[1:round(propTopIndsPerYear*nrow(y)),]
      return(y)
    })
    mydataSub <- do.call(rbind, yearSplit)
    ###
    mydataSub$environment <- as.factor(mydataSub$environment)
    mydataSub$designation <- as.factor(mydataSub$designation)
    
    #Get intercept and fixed effect part from BLUPs
    ## pull the actual formulas MTA stored for this trait
    fixed_formula_chr  <- modelingInput$value[modelingInput$trait == iTrait &
                                                modelingInput$parameter == "fixedFormula"]
    random_formula_chr <- modelingInput$value[modelingInput$trait == iTrait &
                                                modelingInput$parameter == "randomFormula"]
    
    fixed_main <- get_fixed_main_terms(fixed_formula_chr)
    rand_ints  <- get_random_interactions_with_designation(random_formula_chr)
    
    has_main_desig <- grepl("grp\\([^)]*\\bdesignation\\b(?![^)]*[:_])", random_formula_chr[1],
                            perl = TRUE)
    
    
    # collect all non-designation factors that co-occur with designation in random
    co_factors <- unique(unlist(lapply(rand_ints, function(v) setdiff(v, "designation"))))
    if (!length(co_factors)) co_factors <- character()
    
    # overlap with fixed MAIN effects (mirrors MTA's feToAdd intersection by names)
    overlap_fixed <- intersect(co_factors, fixed_main)
    
    # intercept μ for this trait
    mu <- as.numeric(mu_by_trait[iTrait]); if (!is.finite(mu)) mu <- 0
    
    sum_dev <- 0
    
    if (!has_main_desig && length(overlap_fixed)) {
      for (v in overlap_fixed) {
        blues_v <- get_fixed_blues_table(iTrait, v)
        if (!nrow(blues_v)) next
        # MTA writes fixed BLUEs including μ; subtract μ to get deviations
        blues_v$dev <- blues_v$predictedValue - mu
        levs <- choose_levels_for_factor(mydataSub, v, blues_v)
        if (!length(levs)) next
        dev_v <- mean(blues_v$dev[match(levs, blues_v$designation)], na.rm = TRUE)
        if (is.finite(dev_v)) sum_dev <- sum_dev + dev_v
      }
    }
    
    fixed_part <- mu + sum_dev
    centered <- mydataSub$predictedValue - fixed_part
    mydataSub$predictedValue.d <- fixed_part + centered/mydataSub$rel
    
    for(iFt in fixedTerm){mydataSub[,iFt] <- as.numeric(mydataSub[,iFt])}

    # do analysis
    if(!is.na(var(mydataSub[,"predictedValue"],na.rm=TRUE))){ # if there's variance
      if( var(mydataSub[,"predictedValue"], na.rm = TRUE) > 0 ){
        checks <- as.character(unique(mydataSub[grep("check",mydataSub[,"entryType"], ignore.case = TRUE),gTerm]))
        ranran <- "~NULL"
        if(deregress){
          mydataSub$predictedValue <- mydataSub$predictedValue.d
          if(verbose){
            print("Deregressing predicted values using the reliability. We detected that you are providing BLUPs.")
          }
        }else{
          if(verbose){
            print("Using predicted values directly. We detected that you are providing BLUEs.")
          }
        }
        fix <- paste("predictedValue ~",paste(fixedTerm, collapse=" + "))
        ranres <- "~units"#"~dsum(~units | environment)"
        mydataSub=mydataSub[with(mydataSub, order(environment)), ]
        mydataSub$w <- 1/(mydataSub$stdError)
        # remove extreme outliers or influential points
        hh<-split(mydataSub,mydataSub[,fixedTerm])
        #hh <- lapply(hh,function(x){
        #  outlier <- boxplot.stats(x=x[, "predictedValue"],coef=1.5 )$out
        #  bad <- which(x$predictedValue %in% outlier)
        #  if(length(bad) >0){out <- x[-bad,]}else{out<-x}
        #  return(out)
        #})

        if(partition){
          p1 <- p2 <- p2b <- p3 <- p3b <- p4 <- p5 <- p6 <- p7 <- p8 <- numeric();cc <- 1
          for(u in 1:(length(hh))){
            for(w in 1:u){
              if(u != w){
                mydataSub2 <- do.call(rbind,hh[c(u,w)])
                mix <- lm(as.formula(fix), data=mydataSub2)
                sm <- summary(mix)
                p1[cc] <- sm$coefficients[2,1]*ifelse(deregress,deregressWeight,1)
                baselineFirstYear <- mix$coefficients[1] + ( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1))*min(as.numeric(mydataSub[which(mydataSub$trait == iTrait),fixedTerm]), na.rm=TRUE ))
                baselineAverageYear <- mix$coefficients[1] + ( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1))*mean(as.numeric(mydataSub[which(mydataSub$trait == iTrait),fixedTerm]) , na.rm=TRUE ))
                p2[cc] <- round(( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1)) /baselineFirstYear) * 100,3)
                p2b[cc] <- round(( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1)) /baselineAverageYear) * 100,3)
                p3[cc] <- round((sm$coefficients[2,2]/baselineFirstYear) * 100,3)
                p3b[cc] <- round((sm$coefficients[2,2]/baselineAverageYear) * 100,3)
                p4[cc] <- sm$coefficients[1,1]
                p5[cc] <- sm$coefficients[2,2]
                p6[cc] <- sm$coefficients[1,2]
                p7[cc] <- sm$r.squared
                p8[cc] <- 1 - pf(sm$fstatistic[1], df1=sm$fstatistic[2], df2=sm$fstatistic[3])
                cc <- cc+1
              }
            }
          }
          gg <- median(p1, na.rm=TRUE); ggPercentage <- median(p2, na.rm=TRUE); ggPercentageAverageYear <- median(p2b, na.rm=TRUE);
          seGgPercentage <- median(p3, na.rm=TRUE); seGgPercentageAverageYear <- median(p3b, na.rm=TRUE)
          inter <- median(p4, na.rm=TRUE); seb1 <- median(p5, na.rm=TRUE); seb0<- median(p6, na.rm=TRUE)
          r2 <- median(p7, na.rm=TRUE); pv <- median(p8, na.rm=TRUE)
        }else{
          mydataSub <- do.call(rbind, hh)
          mix <- lm(as.formula(fix), data=mydataSub)
          sm <- summary(mix)
          gg <- sm$coefficients[2,1]*ifelse(deregress,deregressWeight,1)
          baselineFirstYear <- mix$coefficients[1] + ( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1))*min(as.numeric(mydataSub[which(mydataSub$trait == iTrait),fixedTerm]) , na.rm=TRUE ))
          baselineAverageYear <- mix$coefficients[1] + ( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1))*mean(as.numeric(mydataSub[which(mydataSub$trait == iTrait),fixedTerm]) , na.rm=TRUE ))
          ggPercentage <- round(( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1)) /baselineFirstYear) * 100,3)
          ggPercentageAverageYear <- round(( (mix$coefficients[2]*ifelse(deregress,deregressWeight,1)) /baselineAverageYear) * 100,3)
          seGgPercentage <- round((sm$coefficients[2,2]/baselineFirstYear) * 100,3)
          seGgPercentageAverageYear <- round((sm$coefficients[2,2]/baselineAverageYear) * 100,3)
          inter <- sm$coefficients[1,1]
          seb1 <- sm$coefficients[2,2]
          seb0 <- sm$coefficients[1,2]
          r2 <- sm$r.squared
          pv <- 1 - pf(sm$fstatistic[1], df1=sm$fstatistic[2], df2=sm$fstatistic[3])
        }
        gg.y1<- sort(unique(mydataSub[,fixedTerm]), decreasing = FALSE)[1] # first year
        gg.yn <- sort(unique(mydataSub[,fixedTerm]), decreasing = TRUE)[1] # last year
        ntrial <- phenoDTfile$predictions # number of trials
        ntrial <- ntrial[which(ntrial$analysisId ==analysisId),]
        ntrial <- ntrial[which(ntrial$trait ==iTrait),]
        ntrial <- ntrial[which(ntrial$effectType =="environment"),]
        ntrial <- length(unique(ntrial$designation))
        phenoDTfile$metrics <- rbind(phenoDTfile$metrics,
                                     data.frame(module="rgg",analysisId=rggAnalysisId, trait=iTrait, environment="across",
                                                parameter=c("ggSlope","ggInter", "gg%(first.year)","gg%(average.year)","r2","pVal","nTrial","initialYear","lastYear"), method=ifelse(deregress,"blup+dereg","mackay"),
                                                value=c(gg,inter, ggPercentage,ggPercentageAverageYear, r2, pv, ntrial,gg.y1,gg.yn  ),
                                                stdError=c(seb1,seb0,seGgPercentage,seGgPercentageAverageYear,0,0,0,0,0)
                                     )
        )
        currentModeling <- data.frame(module="rgg", analysisId=rggAnalysisId,trait=iTrait, environment="across",
                                      parameter=c("deregression","partitionedModel","lessThanFiveYears"), value=c(deregress, partition,less_than_five))
        phenoDTfile$modeling <- rbind(phenoDTfile$modeling,currentModeling[,colnames(phenoDTfile$modeling)] )
        myPreds <- mydataSub[,colnames(phenoDTfile$predictions)]
        myPreds$module <- "rgg"
        myPreds$analysisId <- rggAnalysisId
        phenoDTfile$predictions <- rbind(phenoDTfile$predictions, myPreds)
        counter=counter+1
      }
    }
  }
  #########################################
  ## update databases
  ## write the parameters to the parameter database
  newStatus <- data.frame(module="rgg", analysisId=rggAnalysisId, analysisIdName=NA)
  phenoDTfile$status <- rbind( phenoDTfile$status, newStatus[,colnames(phenoDTfile$status)])
  ## add which data was used as input
  modeling <- data.frame(module="rgg",  analysisId=rggAnalysisId, trait=c("inputObject"), environment="general",
                         parameter= c("analysisId"), value= c(analysisId ))
  phenoDTfile$modeling <- rbind(phenoDTfile$modeling, modeling[, colnames(phenoDTfile$modeling)])
  return(phenoDTfile)#
}
