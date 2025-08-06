traitTransformation <- function(
    object= NULL,
    trait=NULL, # per trait
    transformation = NULL,
    verbose=FALSE
){
  if(is.null(object)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(trait)){stop("Please provide traits to be converted", call. = FALSE)}
  if(length(trait)==0){stop("Please select at least one trait conversion to apply this function.", call. = FALSE)}
  if(length(trait) != length(transformation)){stop("The vector of transformations required need to be as many as the number of traits and viceversa.", call. = FALSE)}
  if(any(paste(trait, transformation, sep="_") %in% colnames(object$data$pheno))){
    stop(paste0("Column(s) ", paste0("'", paste(paste(trait, transformation, sep="_")[which(paste(trait, transformation, sep="_") %in% colnames(object$data$pheno) == TRUE)], collapse="', '"), "'")," already exist."), call. = FALSE)
  }    
  cbrt <- function(x){return(x^(1/3))}
  traitOrig <- trait
  transAnalysisId <- as.numeric(Sys.time())
  ###################################
  # loading the dataset
  mydata <- object$data$pheno
  traitToRemove <- character()
  for(k in 1:length(trait)){
    if (!trait[k] %in% colnames(mydata)){
      if(verbose){
        cat(paste0("'", trait[k], "' is not a column in the given dataset. It will be removed from trait list \n"))
      }
      # trait <- trait[-k]
      traitToRemove <- c(traitToRemove,trait[k])
    }
  }
  traitPresent <- setdiff(trait,traitToRemove)
  trait <- trait[which(traitOrig %in% traitPresent)]
  transformation <- transformation[which(traitOrig %in% traitPresent)]
  #####################################
  # transformation
  counter <- 1
  for(iTrait in trait){ # iTrait=trait[1]
    if(verbose){cat(paste("Analyzing trait", iTrait,"\n"))}
    mydata[,paste(iTrait,transformation[counter], sep = "_")] <- do.call(transformation[counter],list(mydata[,iTrait]))
    counter <- counter + 1
  }
  object$data$pheno <- mydata
  object$metadata$pheno <- rbind( object$metadata$pheno, data.frame(parameter="trait", value=paste(trait, transformation, sep="_")))
  ##########################################
  ## update databases
  ## status
  newStatus <- data.frame(module="transP", analysisId=transAnalysisId, analysisIdName="")
  object$status <- rbind( object$status, newStatus[,colnames(object$status)])
  # modeling
  currentModeling <- data.frame(module="transP", analysisId=transAnalysisId,trait=trait, environment=NA,
                                parameter=paste(trait, transformation, sep="_"),
                                value=paste0(transformation,"(",trait, ")"))
  object$modeling <- rbind(object$modeling,currentModeling )
  return(object)
}

freeFunction <- function(object, formula, newName = NULL){
  if(is.null(object)){stop("Please provide the name of the file to be used for analysis", call. = FALSE)}
  if(is.null(formula)){stop("Please provide formula to be used", call. = FALSE)}
  if(newName == ""){stop("Please provide a name for calculated trait", call. = FALSE)}
  if(is.null(newName)){stop("Please provide a name for calculated trait", call. = FALSE)}
  if(newName %in% colnames(object$data$pheno)){
    stop(paste0("Column '", newName,"' already exist. Please provide a different name for calculated trait."), call. = FALSE)
  }
  transAnalysisId <- as.numeric(Sys.time())
  ###################################
  # loading the dataset
  mydata <- object$data$pheno
  f <- eval(parse(text = formula )) # transform text into a function
  # now get variables
  dep <- deparse(f)
  dep <- dep[grep("function",dep)]
  dep <- gsub("function","cbind",dep)
  responsef <- as.formula(paste(dep,"~1"))
  mfna <- try(model.frame(responsef, data = mydata, na.action = na.pass), silent = TRUE)
  if(inherits(mfna,"try-error")) {
    stop("We could not find the variables specified in your function.", call. = FALSE)
  }
  mfna <- as.matrix(mfna) # matrix format to avoid a single column
  colnames(mfna) <- all.vars(responsef) # assign colnames
  mfna <- as.data.frame(mfna) # return back to data.frame format
  # run the function and get the new variable
  newVar <-  try( do.call(f, as.list(mfna) ), silent = TRUE )
  if(inherits(newVar,"try-error")) {
    stop(paste("Computation failed with the following error message: \n\n",newVar[[1]]), call. = FALSE)
  }
  # if(is.null(newName) ){
  #   newName <- paste0(paste(colnames(mfna),collapse = "_"),"_new")
  # }
  # if( newName == "" ){
  #   newName <- paste0(paste(colnames(mfna),collapse = "_"),"_new")
  # }
  mydata[,newName] <- as.vector(newVar)
  # store updated data
  object$data$pheno <- mydata
  object$metadata$pheno <- rbind( object$metadata$pheno, data.frame(parameter="trait", value=newName))
  ##########################################
  ## update databases
  ## status
  newStatus <- data.frame(module="transF", analysisId=transAnalysisId, analysisIdName="")
  object$status <- rbind( object$status, newStatus[,colnames(object$status)])
  # modeling
  currentModeling <- data.frame(module="transF", analysisId=transAnalysisId,trait="inputObject", environment="general",
                                parameter=newName,
                                value=formula)
  object$modeling <- rbind(object$modeling,currentModeling )
  
  return(object)
}

