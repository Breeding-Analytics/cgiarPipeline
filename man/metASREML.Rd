\name{metASREML}
\alias{metASREML}

\title{
Multi-environmental trial analysis pipeline using ASReml package
}
\description{
Perform the multi-trial genetic evaluation according to a custom-based formula.
}
\usage{
metASREML(
    phenoDTfile= NULL, analysisId=NULL,
    fixedTerm= list("1"),  randomTerm=NULL, covMod=NULL, addG=NULL, nFA=NULL,
    envsToInclude=NULL, trait= NULL, traitFamily=NULL, useWeights=TRUE,
    calculateSE=TRUE, heritLB= 0.15,  heritUB= 0.95,
    meanLB=0, meanUB=Inf, maxIters=50,  verbose=TRUE
)
}

\arguments{
  \item{phenoDTfile}{
The data object structure produced from bioflow which contains data, metadata, modifications and tables of results (metrics, predictions, modeling, status), in order to pick the necessary information to calculate an index.
}
\item{analysisId}{
A numeric id to link to the predictions to be used for the analysis.
}
\item{fixedTerm}{
A list specifying the terms to be considered fixed. If an element in the list has more than one term is considered an interaction between those terms.
}
\item{randomTerm}{
A list specifying the terms to be considered random. If an element in the list has more than one term is considered an interaction between those terms.
}
\item{covMod}{
The option of variance and covariance model to fit, the options are: "fa_model, "us_model", "diag_model", "none".
}
\item{addG}{
The option that add type of effect should be used to calculate the covariance between levels of designation, the optios are: "geno_model", "pedigree_model", "genoAD_model", "none".
}
\item{nFA}{
A numeric term that indicate the number of factors to be used in the model in the case that covMod option will be "fa_model".
}
\item{envsToInclude}{
A table of fields (rows) by traits (columns) to specify which fields should be included in the analysis for each trait.
}
\item{trait}{
A character vector with the traits that should be fitted.
}
\item{traitFamily}{
A character vector specifying the family distribution for each trait in case normal gaussian is not a good fit.
}
\item{useWeights}{
A TRUE or FALSE statement to specify if the weights using the standard error should be applied. The default is TRUE since we assume the two-stage model.
}
\item{calculateSE}{
A TRUE or FALSE statement to specify if the standard errorrs for the random effects should be calculated or not.
}
\item{heritLB}{
A numeric value or vector for trait heritability and reliablity lower-bound-threshold(s). Any trait by environment combination below this value will be excluded from the analysis.
}
\item{heritUB}{
A numeric value or vector for trait heritability and reliablity upper-bound-threshold(s). Any trait by environment combination above this value will be excluded from the analysis.
}
\item{meanLB}{
A numeric value or vector for trait mean lower-bound-threshold(s). Any trait by environment combination below this value will be excluded from the analysis.
}
\item{meanUB}{
A numeric value or vector for trait mean upper-bound-threshold(s). Any trait by environment combination above this value will be excluded from the analysis.
}
\item{maxIters}{
Maximum number of restricted maximum likelihood iterations to be run for each trait.
}
\item{verbose}{
Should the logs or warninig messages of the analysis be printed?
}

}
\details{
The core algorithm of the genetic evaluation using the two-step approach is the multi-trial analysis. This option aims to model fixed and random effects using the results from the single trial (weighted by the standard errors) analysis and optionally include any relationship matrix between levels of the random effects using the pedigree, or marker data. If the user wants to build trait-based relationship matrices is also possible.

}
\value{
If arguments are properly specified the function returns:
\describe{
\item{result}{the same input object with additional rows in the predictions, metrics, modeling and status tables.}
}
}

\references{

Finlay, K. W., & Wilkinson, G. N. (1963). The analysis of adaptation in a plant-breeding programme. Australian journal of agricultural research, 14(6), 742-754.

Henderson Jr, C. R. (1982). Analysis of covariance in the mixed model: higher-level, nonhomogeneous, and random regressions. Biometrics, 623-640.

Odegard, J., Indahl, U., Stranden, I., & Meuwissen, T. H. (2018). Large-scale genomic prediction using singular value decomposition of the genotype matrix. Genetics Selection Evolution, 50(1), 1-12.

}

\examples{
####=========================================####

data(DT_example)
result$status
result$metadata$pheno

final <- metASREML(
    phenoDTfile= result, analysisId=result$status$analysisId[3],
    fixedTerm= list("environment"),
    randomTerm=list("designation"), 
    covMod="none", addG=NULL, nFA=NULL,
    envsToInclude=NULL, trait= "Yield_Mg_ha", traitFamily=NULL, useWeights=TRUE,
    calculateSE=TRUE, heritLB= 0.15,  heritUB= 0.95,
    meanLB=0, meanUB=Inf, maxIters=50,  verbose=TRUE
)

final$status
tail(final$predictions,10)
tail(final$metrics,10)

pick <- which(final$predictions$analysisId == final$status$analysisId[length(final$status$analysisId)])
plot(final$predictions[pick,"predictedValue"], col=as.factor(final$predictions[pick,"environment"]))

}

