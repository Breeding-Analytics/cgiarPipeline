\name{baseIndex}
\alias{baseIndex}

\title{
Base selection index pipeline 
}
\description{
Calculates the base selection index given a table of traits and weights in traits (see Details section).
}
\usage{
baseIndex(
  phenoDTfile, 
  analysisId,   
  analysisIdName,
  trait,  
  weights
)
}

\arguments{
  \item{phenoDTfile}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisId}{
A numeric id to link to the predictions to be used for the selection index.
}
\item{analysisId}{
A character id to personalized the analysis ID.
}
\item{trait}{
A character vector with the traits that should be included in the selection index. Once the predictions are subset to the requested analysisId an additional subset for the traits of interest is performed.  
}
\item{weights}{
A numeric vector with trait weights. 
}

}
\details{

The base index is just an average weight where a matrix of traits Y of dimensions n rows (individuals) and t columns (traits) is multiplied by weights to compute the vector m of total merit of individuals:

m = Y w

Then the function row binds a new table of predictions to the predictions table of the data object together with some additional information to the modeling, metrics and status tables. 

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Ceron-Rojas, J. J., & Crossa, J. (2018). Linear selection indices in modern plant breeding (p. 256). Springer Nature.

}

\examples{
####=========================================####
data(DT_example)

result$status
result$metadata$pheno

final <- baseIndex(
  phenoDTfile= result, 
  analysisId=result$status$analysisId[4],
  analysisIdName=result$status$analysisIdName[4],
  trait= c("Plant_Height_cm", "Ear_Height_cm", "Yield_Mg_ha"),
  weights = c(0, 0 , 1)
)

tail(final$predictions)

}

