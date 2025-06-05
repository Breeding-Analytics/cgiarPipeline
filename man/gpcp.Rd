\name{gpcp}
\alias{gpcp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomic Prediction of Cross Performance Pipeline
}
\description{
Performs genomic prediction of cross performance based on Werner et al. (2023) together with optimal cross selection based on Woolliams et al. (2015).
}
\usage{
gpcp(
    phenoDTfile= NULL,
    analysisId=NULL,
    relDTfile= NULL,
    trait= NULL,
    environment="across",
    nCross=20,
    targetAngle=30,
    verbose=FALSE,
    maxRun=100,
    numberBest=100,
    entryType=NULL
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{phenoDTfile}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisId}{
A numeric id to link to the predictions to be used for the selection index.
}
\item{relDTfile}{
One of the following; GRM, NRM single-step relationship matrix.
}
\item{trait}{
A character vector with the trait that should be fitted.
}
\item{environment}{
If the user wants to use predictions from an specific environment. In NULL all are used.
}
\item{nCross}{
Number of top crosses to be selected.
}
\item{targetAngle}{
Target angle defining the trade-off between performance and diversity. Zero degrees is weighting strongly towards performance. Ninety degrees is weighting strongly towards diversity.
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}
\item{maxRun}{
Maximum number of runs (iterations) without change in the genetic algorithm.
}
\item{numberBest}{
The complexity and computation time of the algorithm scales up with greater number of individuals used for predicted crosses. This arguments applies a filter to only use the top N individuals for the trait of interest.
}
\item{entryType}{
Which entry types should be used in the algorithm.
}
}

\details{
The genetic evaluation approach we use known as 'two-step' first analyze trait by trait and trial by trial to remove the spatial noise from experiments using experimental factors like blocking and spatial coordinates. Each trial is one level of the environment column (defined when the user matches the expected columns to columns present in the initial phenotypic input file). Genotype is fitted as both, fixed and random. The user defines which should be returned in the predictions table. By default genotype (designation column) predictions and their standard errors are returned.
}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
\item{predictions}{the predictions table in this particular case returns in the stdError and the reliability columns the inbreeding coefficient and the OCS optimization value respectively.}
}
}
\references{

Werner, C. R., Gaynor, R. C., Sargent, D. J., Lillo, A., Gorjanc, G., & Hickey, J. M. (2023). Genomic selection strategies for clonally propagated crops. Theoretical and Applied Genetics, 136(4), 74.

Kinghorn, B. (1999). 19. Mate Selection for the tactical implementation of breeding programs. Proceedings of the Advancement of Animal Breeding and Genetics, 13, 130-133.

https://alphagenes.roslin.ed.ac.uk/wp/wp-content/uploads/2019/05/01_OptimalContributionSelection.pdf?x44213

Woolliams, J. A., Berg, P., Dagnachew, B. S., & Meuwissen, T. H. E. (2015). Genetic contributions and their optimization. Journal of Animal Breeding and Genetics, 132(2), 89-99.

}
\examples{
####=========================================####
data(DT_example)

result$status
result$metadata$pheno
table(result$predictions$trait, result$predictions$analysisId)

final <- gpcp(
    phenoDTfile= result, # analysis to be picked from predictions database
    analysisId=result$status$analysisId[4],
    relDTfile= "grm", # "nrm", "grm", "both"
    trait= "desireIndex", # per trait
    environment="across",
    nCross=20,
    targetAngle=30, # in radians
    verbose=FALSE,
    maxRun=100,
    numberBest=100,
    entryType=NULL
)

tail(final$predictions)

}

