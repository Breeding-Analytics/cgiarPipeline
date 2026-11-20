\name{metSommer}
\alias{metSommer}
\title{Multi-environment analysis using sommer}
\description{
Fit second-stage genetic evaluation models without an ASReml license, using
SOA adjusted means, their standard errors, and optional relationship kernels.
}
\usage{
metSommer(phenoDTfile = NULL, analysisId = NULL, analysisIdgeno = NULL,
  fixedTerm = list("environment"), randomTerm = list("designation"),
  covMod = NULL, addG = NULL, nFA = NULL, envsToInclude = NULL,
  trait = NULL, traitFamily = NULL, useWeights = TRUE,
  calculateSE = TRUE, heritLB = 0.15, heritUB = 0.95,
  meanLB = 0, meanUB = Inf, maxIters = 50, verbose = TRUE)
}
\arguments{
\item{phenoDTfile}{Bioflow object containing data, metadata, status,
predictions, metrics and modeling tables.}
\item{analysisId}{One or more SOA/STA analysis IDs. Overlapping genotype by
environment predictions for a trait are rejected.}
\item{analysisIdgeno}{Genotype QA/QC stamp, required for genomic kernels.
The selected imputed genlight version or selected matrix modifications are used.}
\item{fixedTerm, randomTerm}{Lists of variable names. An element containing
several names denotes an interaction. Fixed \code{"none"} or \code{"1"}
denotes an intercept-only model. Numeric variables remain numeric.}
\item{covMod}{List parallel to \code{randomTerm}, with one structure per
variable: \code{"none"}, \code{"Structure model_diag"},
\code{"Structure model_us"}, \code{"Structure model_fa"},
\code{"Relationship structure_GenoA"}, \code{"Relationship structure_GenoD"},
\code{"Relationship structure_GenoAD"}, \code{"Relationship structure_Pedigree"},
or \code{"weatherInfo"}. Dotted \code{"none."} menu placeholders are accepted.}
\item{addG}{Compatibility alias for \code{covMod}, used only if \code{covMod}
is missing.}
\item{nFA}{List with one FA rank per random term; defaults to one.}
\item{envsToInclude}{Environment by trait matrix or data frame. Positive cells
include environments. Environment row names and trait column names are required.}
\item{trait}{Traits to analyze, including configured TPP aliases.}
\item{traitFamily}{One label per trait: \code{"gaussian"}, \code{"binomial"},
\code{"Gamma"} or \code{"poisson"}. The corresponding ASReml UI labels are
also accepted, without evaluating R code. Default is Gaussian.}
\item{useWeights}{Use precision weights \eqn{1/SE^2}. For Gaussian second-stage
fits the residual multiplier is fixed at one. Without weights it is estimated.}
\item{calculateSE}{Export prediction standard errors and reliabilities.
When false these output columns contain \code{NA}.}
\item{heritLB, heritUB, meanLB, meanUB}{Inclusive SOA reliability/heritability
and mean bounds, each scalar or length equal to the number of traits.
Heritability metrics must be available; if SOA mean metrics are unavailable,
environment means are computed from the selected SOA predictions.}
\item{maxIters}{Maximum fitting iterations.}
\item{verbose}{Print fitting progress. Warnings, including nonconvergence,
remain visible regardless of this setting.}
}
\details{
Requires sommer 4.5.25 or newer. Models use the public \code{mmes},
\code{vsm}, \code{dsm}, \code{usm} and \code{fam} APIs, with no ASReml calls.
Gaussian models use REML; the other families use PQL with their canonical
UI links (identity, logit, inverse and log, respectively).

Main effects, compound symmetry, diagonal, main plus diagonal,
Finlay-Wilkinson, factor analytic, additive/dominance and parental GCA/SCA
models are specified through the same list interface as \code{metASREML}.
\code{envIndex} is the centered environmental mean for each trait.
Pedigree kernels use mapped pedigree columns through pedigreemm, independently
of marker availability. Marker and weather kernels are scaled to mean
diagonal one and receive a small diagonal regularization. GenoAD retains
the ASReml pipeline's centered log-dosage kernel convention.

SOA BLUP inputs are deregressed only for stamps recorded as BLUP for the
corresponding trait; their reliabilities must be in \eqn{(0,1]}.
TPP trait and environment mappings in \code{metadata$tpp_analysis_config}
are applied before fitting, and output retains the TPP trait label.
Invalid inputs and empty filtered traits fail explicitly without partial results.

Merits use an equally weighted environment by entry prediction grid, including
unobserved cells, when covariates map uniquely to entries or environments.
For covariates varying within these cells, a warning announces observed-record
marginal predictions instead. Joint contrast PEVs account for covariance
between components and environments; standard errors are not averaged.
Reliability is bounded \eqn{1-PEV/G} on the latent genetic scale, excluding
non-genetic random terms. It is unavailable when no genetic random variance
is fitted, including models with only fixed genotype effects.
For non-Gaussian families merits are inverse-linked to the response scale
and their standard errors use the delta method. Fixed coefficients, component
BLUPs, parental GCA and regression slopes remain on the link scale.
Fixed genotype contrasts are exported as \code{designation_BLUE}, separately
from the genotype merits that retain these fixed effects. Finlay-Wilkinson
fits export both slopes and environment-specific merits.
REML/PQL estimates and convergence are engine-specific, not guaranteed to
be numerically identical to ASReml results.
}
\value{
The input object with additional \code{moaSommer} rows in \code{predictions},
\code{metrics}, \code{modeling} and \code{status}. The input is not modified
in place; no files or global variables are written.
}
\seealso{\code{\link{metASREML}}, \code{\link{metLMMsolver}}}
\examples{
\dontrun{
data(DT_example)
ids <- result$status$analysisId[result$status$module \%in\% c("soa", "sta")]
fit <- metSommer(result, analysisId = ids[1], trait = "Yield_Mg_ha",
  randomTerm = list("designation", c("designation", "environment")),
  covMod = list("none", c("none", "Structure model_diag")))
}
}
