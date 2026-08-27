\name{pgg}
\alias{pgg}

\title{
Predicted genetic gain pipeline.
}
\description{
Performs a predicted genetic gain estimation.
}
\usage{
pgg(
    phenoDTfile= NULL,
    analysisId=NULL,
    trait=NULL, # per trait
    by=NULL,
    percentage=10,
    cycleTime=NULL,
    verbose=TRUE
)
}

\arguments{
  \item{phenoDTfile}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisId}{
A numeric id to link to the predictions to be used for the analysis.
}
\item{trait}{
A character vector with the traits that should be analyzed.
}
\item{by}{
The column by which partitions of the data should be built to extract the metrics from breeders' equation. 
}
\item{percentage}{
Percentage of selected individuals to assume.
}
\item{cycleTime}{
Time in years to complete one breeding cycle (from crossing until the selection of the new parents). Can be several years (e.g. 5) or less than a year (e.g. 0.5). Used to express the predicted genetic gain as a per-year rate (PGG/year = PGG / cycleTime). Defaults to 1 if not provided.
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}

}
\details{

In order to monitor the efficacy of genetic evaluation in the current cycle of selection the predicted genetic gain formula is used. This option aims to calculate the predicted genetic gain from the classical breeders' equation R = i*r*s being R the response to selection, i the selection intensity, r the selection accuracy, s the genetic standard deviation.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Lush, J. L. (2013). Animal breeding plans. Read Books Ltd.

Mrode, R. A. (2014). Linear models for the prediction of animal breeding values. Cabi.

}

\examples{
####=========================================####
data(DT_example)

result$status
result$metadata$pheno

final <- pgg(
    phenoDTfile= result,
    analysisId=result$status$analysisId[3],
    trait="Plant_Height_cm", # per trait
    by=NULL, # if NULL 'environment' will be assumed
    percentage=10,
    verbose=TRUE
)

tail(final$metrics)

}

