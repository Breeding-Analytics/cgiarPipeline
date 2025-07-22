#accumulative_temp
#accumulative_rh
#accumulative_ppt
#accumulative_dpv (posible)
# library(nasapower)
# library(pacman)

# p_load(nasapower,readr,tidyverse,plotly)
#
# # head(result$data$pheno)
# # unique(result$data$pheno[,c("Env","Year")])
# #
# # # Upload pheno data -------------------------------------------------------------
# # bioflow <- read_csv("D:/OneDrive - CGIAR/2024/HACKATHON FEB 2024/bioflow.csv")
# # n<-names(bioflow %>% dplyr::select(-all_of(trait)))
# # #unique(bioflow$Field_Location)
#
# # Parameters --------------------------------------------------------------
# trait<-c("Yield_Mg_ha")
# #f_location<-c("NYH3")
# #y<-c(2021)
# LAT<- -27.48
# LONG<-151.81
#
# obj<-data()
# df<-obj$data$pheno
#
# date_planted<-min(df$Date_Planted)
# date_harvest<-max(df$Date_Harvested)
#
# temporal<-"hourly" #daily,monthly

# Funtion for download Climate from NASA  ---------------------------------
nasaPowerExtraction <- function(LAT,LONG,date_planted,date_harvest,environments, temporal="hourly"){

# Climate data ------------------------------------------------------------
  wthList <- metaList <- list()
  for(iEnv in 1:length(environments)){ # iEnv=1
    if(temporal == "hourly"){
      prov <-nasapower::get_power(community = "ag", #c("ag", "re", "sb")
                                  lonlat = c(LONG[iEnv],LAT[iEnv]), #Decimal degrees
                                  pars = c("RH2M","T2M","PRECTOTCORR"),# "T2M_MAX","T2M_MIN"
                                  dates = c(date_planted[iEnv], date_harvest[iEnv]),#YYYY-MM-DD
                                  temporal_api = temporal) %>%
        dplyr::mutate(datetime=ISOdate(YEAR, MO, DY,HR),
                      date=as.Date(datetime))
    }else if(temporal == "daily"){
      prov <-nasapower::get_power(community = "ag", #c("ag", "re", "sb")
                                  lonlat = c(LONG[iEnv],LAT[iEnv]), #Decimal degrees
                                  pars = c("RH2M","T2M","PRECTOTCORR"),# "T2M_MAX","T2M_MIN"
                                  dates = c(date_planted[iEnv], date_harvest[iEnv]),#YYYY-MM-DD
                                  temporal_api = temporal) %>%
        dplyr::mutate(datetime=ISOdate(YEAR, MM, DD),
                      date=as.Date(datetime))
      colnames(prov) <- cgiarBase::replaceValues(colnames(prov), Search = c("MM","DD"), Replace = c("MO","DY") )
    }else if(temporal == "monthly"){

      timeDiff1 <- ( as.Date("2022/12/31") - date_harvest[iEnv]  )# if users selects a date newer than available, go to latest available
      if(timeDiff1 < 0){
        date_harvest[iEnv] <- date_harvest[iEnv] + timeDiff1 - 1
        date_planted[iEnv] <- date_planted[iEnv] + timeDiff1 - 1
      }
      timeDiff2 <- date_harvest[iEnv] - date_planted[iEnv]
      if(timeDiff1 < 365){
        date_planted[iEnv] <- date_planted[iEnv] - (365 - timeDiff2)
      }
      prov <-nasapower::get_power(community = "ag", #c("ag", "re", "sb")
                                  lonlat = c(LONG[iEnv],LAT[iEnv]), #Decimal degrees
                                  pars = c("RH2M","T2M","PRECTOTCORR"),# "T2M_MAX","T2M_MIN"
                                  dates = c(date_planted[iEnv], date_harvest[iEnv]),#YYYY-MM-DD
                                  temporal_api = temporal) #%>%
        # dplyr::mutate(datetime=ISOdate(YEAR, MM, DD),
        #               date=as.Date(datetime))
      provList <- split(prov, prov$PARAMETER)
      provB <- lapply(provList, function(x){ # reshape by parameter
        reshape(x[,setdiff(colnames(x), c("PARAMETER","ANN"))], idvar = c("LON","LAT","YEAR"), varying = list(4:15),
                v.names = as.character(unique(x[,"PARAMETER"])) , direction = "long",
                timevar = "MM" )
      })
      prov <- provB[[1]]
      if(length(provB) > 1){
        for(i in 2:length(provB)){
          prov <- merge(prov, provB[[i]], by=c("LON","LAT","YEAR","MM"), all.x = TRUE )
        }
      }
      prov$datetime <- ISOdate(prov$YEAR, prov$MM, 1)
      prov$date <- as.Date(prov$datetime)

    }

    prov$environment <- environments[iEnv]
    # metadata
    meta <- data.frame(parameter=c("trait","trait","trait","latitude","longitude","year","month","day", "date","environment"),
                       value= c("RH2M","T2M","PRECTOTCORR","LAT","LON","YEAR","MO","DY","date","environment")
                       )
    # meta <- data.frame( environment=environments[iEnv], trait=c("RH2M","T2M","PRECTOTCORR",   "RH2M","T2M","PRECTOTCORR",    "latitude", "longitude", "plantingDate","harvestingDate"),
    #             parameter=c(rep("mean",3), rep("sd",3), c("coordinate", "coordinate", "date", "date")),
    #             value= c( apply(prov[,c("RH2M","T2M","PRECTOTCORR")],2,mean, na.rm=TRUE), apply(prov[,c("RH2M","T2M","PRECTOTCORR")],2,sd, na.rm=TRUE), c(LAT[iEnv], LONG[iEnv], date_planted[iEnv], date_harvest[iEnv] ) )
    # )
    # save
    wthList[[iEnv]] <- as.data.frame(prov)
    metaList[[iEnv]] <-  meta
  }
  WTH <- unique(do.call(rbind,wthList))
  META <- unique(do.call(rbind, metaList))

# descriptive  --------------------------------------------------------
  # descriptive<-summary(dplyr::select(WTH,RH2M,T2M,PRECTOTCORR))

# Outputs -----------------------------------------------------------------
  output<- list(data = WTH,
                metadata = META)
  return(output)
}



summaryWeather <- function(object, wide=FALSE){

  if(is.null(object$data$weather)){
    out <- as.data.frame(matrix(nrow=0, ncol=4));  colnames(out) <- c("environment", "trait", "parameter" ,  "value")
    ### summarize environmental indices
    data2 <- object$data$pheno
    metadata2 <- object$metadata$pheno
    if(!is.null(metadata2)){
      #check if each column has all missing values
      all_miss <- apply(data2, 2, function(x) all(is.na(x)))
      #display columns with all missing values
      all_miss_var <- names(all_miss[all_miss>0])

      traits2 <- setdiff(metadata2[which(metadata2$parameter == "trait"),"value"],all_miss_var)
      environ2 <- metadata2[which(metadata2$parameter == "environment"),"value"]
      provList2 <- list()
      for(iTrait2 in traits2){ # iTrait2 = traits2[1]
        data2$traitUsed <- data2[,iTrait2]
        ei <- aggregate(as.formula(paste("traitUsed","~environment")), data=data2,FUN=mean, na.rm=TRUE);
        colnames(ei)[2] <- "value"
        ei$trait <- iTrait2
        ei$parameter <- "mean"
        #
        ei2 <- ei
        ei2$value <- ei$value - mean(ei$value, na.rm=TRUE)
        ei2$parameter <- "envIndex" # paste0(iTrait,"-envIndex")
        provList2[[iTrait2]] <- rbind(ei,ei2)
      }
      out <- rbind(out, do.call(rbind,provList2) )
    }
  }else{
    data <- object$data$weather
    metadata <- object$metadata$weather

    traits <- metadata[which(metadata$parameter == "trait"),"value"]
    environ <- metadata[which(metadata$parameter == "environment"),"value"]
    lat <- metadata[which(metadata$parameter == "latitude"),"value"]
    lon <- metadata[which(metadata$parameter == "longitude"),"value"]
    ## summarize traits
    provList <- list()
    for(iTrait in c(traits, lat, lon) ){ # iTrait = traits[1]
      # mean
      prov <- aggregate(as.formula( paste(iTrait, "~", "environment") ), FUN=mean, na.rm=TRUE, data=data)
      colnames(prov)[2] <- "value"
      prov$trait <- iTrait
      prov$parameter <- "mean"
      # sd
      prov2 <- aggregate(as.formula( paste(iTrait, "~", "environment") ), FUN=sd, na.rm=TRUE, data=data)
      colnames(prov2)[2] <- "value"
      prov2$trait <- iTrait
      prov2$parameter <- "sd"
      #
      provList[[iTrait]] <- rbind(prov,prov2)
    }
    out <- do.call(rbind,provList)

    ### summarize environmental indices
    data2 <- object$data$pheno
    metadata2 <- object$metadata$pheno
    if(!is.null(metadata2)){
      #check if each column has all missing values
      all_miss <- apply(data2, 2, function(x) all(is.na(x)))
      #display columns with all missing values
      all_miss_var <- names(all_miss[all_miss>0])

      traits2 <- setdiff(metadata2[which(metadata2$parameter == "trait"),"value"],all_miss_var)
      environ2 <- metadata2[which(metadata2$parameter == "environment"),"value"]
      provList2 <- list()
      for(iTrait2 in traits2){ # iTrait2 = traits2[1]
        ei <- aggregate(as.formula(paste(iTrait2,"~environment")), data=data2,FUN=mean, na.rm=TRUE);
        colnames(ei)[2] <- "value"
        ei$trait <- iTrait2
        ei$parameter <- "mean"
        #
        ei2 <- ei
        ei2$value <- ei$value - mean(ei$value, na.rm=TRUE)
        ei2$parameter <- "envIndex" # paste0(iTrait,"-envIndex")
        provList2[[iTrait2]] <- rbind(ei,ei2)
      }
      out <- rbind(out, do.call(rbind,provList2) )
    }

  }

  if(wide){
    out$tp <- paste(out$trait, out$parameter)
    out <- reshape(out[,c("environment","tp","value")], direction = "wide", idvar = "environment",
                       timevar = "tp", v.names = "value", sep= "_")
    rownames(out) <- out$environment
    if(nrow(out) > 1){
      Z = model.matrix(~environment-1, data=out); colnames(Z) <- gsub("ironment","",colnames(Z))
    }else{
      Z <- matrix(1,1,1); colnames(Z) <- paste0("env",out$environment)
    }

    out <- cbind(out,Z)
    # colnames(out) <- gsub("[[:punct:]]", "", colnames(out) )
    colnames(out) <- gsub("[^[:alnum:]]","",colnames(out))
    out <- out[,-1]
  }

  # meta <- data.frame( environment=environments[iEnv], trait=c("RH2M","T2M","PRECTOTCORR",   "RH2M","T2M","PRECTOTCORR",    "latitude", "longitude", "plantingDate","harvestingDate"),
  #             parameter=c(rep("mean",3), rep("sd",3), c("coordinate", "coordinate", "date", "date")),
  #             value= c( apply(prov[,c("RH2M","T2M","PRECTOTCORR")],2,mean, na.rm=TRUE), apply(prov[,c("RH2M","T2M","PRECTOTCORR")],2,sd, na.rm=TRUE), c(LAT[iEnv], LONG[iEnv], date_planted[iEnv], date_harvest[iEnv] ) )
  # )
  return(out)

}





# d<-np(LAT,LONG,date_planted,date_harvest) #test example
#
# d$TS_RH
#
#
#
#
# # Filter by... ------------------------------------------------------------
# bioflow <- bioflow %>%
#   dplyr::select(all_of(n),trait) %>%
#   dplyr::filter(Year %in% y) %>%
#   #dplyr::filter(Field_Location %in% f_location) %>%
#   dplyr::mutate("pl_date.{trait}":=Date_Planted,
#                 "ev_date.{trait}":=Date_Harvested,
#                 "days.{trait}"   :=Date_Harvested-Date_Planted)
# #table(bioflow$days.Yield_Mg_ha)
# arg<-unique(bioflow$days.Yield_Mg_ha)
# # max(bioflow$Date_Harvested)
# # min(bioflow$Date_Planted)
#
# # Intern parameters -------------------------------------------------------
# date1<-min(bioflow$Date_Planted)
# date2<-max(bioflow$Date_Harvested)
#
# # Upload climate data  ----------------------------------------------------
# WTH<-nasapower::get_power(community = "ag", #c("ag", "re", "sb")
#                       lonlat = c(LONG,LAT), #Decimal degrees
#                       pars = c("RH2M","T2M","PRECTOTCORR"),# "T2M", "PRECTOTCORR","T2M_MAX","T2M_MIN"),
#                       dates = c(date1, date2),#YYYY-MM-DD
#                       temporal_api = "hourly")
# wth<-WTH %>%
#   dplyr::mutate(date=ISOdate(YEAR, MO, DY,HR),
#          date=as.Date(date),"days.{trait}":=date-date1) %>%
#   dplyr::filter(days.Yield_Mg_ha>=min(bioflow$days.Yield_Mg_ha))
#
# #View(wth)
# #names(wth)
#
#
#
# # b_clim <- wth%>%
# #   #dplyr::filter(date<=date2) %>%
# #   mutate(days.Yield_Mg_ha=as.factor(days.Yield_Mg_ha)) %>%
# #   ggplot(aes(x=days.Yield_Mg_ha,y=RH2M,colour=days.Yield_Mg_ha))+
# #   geom_boxplot()+
# #   # geom_label_repel(
# #   #   mapping = aes(label = ifelse(best, Lines, NA)),
# #   #   size = 2,
# #   #   max.overlaps = 50,
# #   #   alpha = 0.9, fill = "yellow")+
# #   #facet_wrap(~TRT,nrow = 2)+
# #   stat_summary(
# #     fun = median,
# #     geom = 'line',
# #     aes(group = days.Yield_Mg_ha , colour = days.Yield_Mg_ha),
# #     position = position_dodge(width = 0.9) #this has to be added
# #   )+xlab(" ")+
# #   ylab("Relative Humidity (%)")+
# #   labs(title="",colour="")+
# #   theme(legend.position = "none",axis.text.x=element_text(angle = 60))
# #
# # b_clim





