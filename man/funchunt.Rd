\name{funchunt}
\alias{funchunt}

\title{
subset of a population
}
\description{
Optimization for select subset of a population
}
\usage{
funchunt(datosgen,dir_fileGen,dir_filePhen,dir_fileDist,size1,w)
}

\arguments{
  \item{datosgen}{
Marker data (data.frame).
}
\item{dir_fileGen}{
Character value, "none" or "used", for know what type of data we will used in the analysis.
}
\item{dir_filePhen}{
Character value, "none" or "used", for know what type of data we will used in the analysis.
}
\item{dir_file_Dist}{
Character value, "none" or "used", for know what type of data we will used in the analysis.
}
\item{size1}{
Numeric value between (0,1) for determinate the size of subset in terms of percent.
}
\item{w}{
Vector with 0's and 1's of size 15, for determinate with algorithm of optimization we will used for select the subset.
}

}
\details{
Can construct cores based on genetic marker data, phenotypic traits or precomputed distance matrices, optimizing one of many provided evaluation measures depending on the precise purpose of the core (e.g. high diversity, representativeness, or allelic richness). 

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{mds_res}{a dataframe with the 3 first factor from the multidimensional scaling (MDS) analysis and the identification of the individuals in the subset.}
\item{statsF}{a dataframe with the basic diversity statistics for population and for subset: He:expected heterocigozity, sd_He:standard deviation of the expected heterocigozity, Ho:observed heterocigozity, sd_Ho:standard deviation of the observed heterocigozity, CV:coeficient of variance , SH_Index: shannon index, sd_SH: standard deviation of the shannon index.}
\item{perctCP12}{a numeric vector with the values of variability explained by each factor from the MDS analysis. }

}
}

\references{
Core Hunter 3, Guy Davenport and Herman De Beukelaer.
}

\examples{
####=========================================####

data(DT_example)
mydata<-as.data.frame(result$data$geno)
paropt=5
ww=rep(0,15)
ww[which(c("ENMR","ENCE","ENGD","EN","ANMR","ANCE","ANGD","AN","EEMR","EECE","EEGD","EE","SH","HE","CV")==paropt)]=1

final<- funchunt(datosgen=mydata,dir_fileGen="used",dir_filePhen="none",dir_fileDist="none",size1=0.3,w=ww)


}

