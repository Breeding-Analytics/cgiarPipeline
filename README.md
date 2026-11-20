
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CGIAR pipelines

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

This package has been developed primarily for computing breeding
analytic pipelines such as single trial analysis, multi-trial analysis,
population structure, etc. The focus of the package is pipelines related
to the use and understanding of evolutionary forces. Is currently
structured by evolutionary forces:

Selection modules

- Genetic evaluation modules

- Selection history modules

Mutation

- Gene discovery modules

- Mutation history modules

Gene flow and drift modules

- Gene frequencies modules

- Drift history modules

QC & Transformation modules

Please go to the package help pages to understand better each of the
modules and methods behind.

## Scope

The current scope is limited to enable Biometrical Genetics methods
recommended by the Quantitative Genetics community of the CGIAR
\[Excellence in Breeding (EiB) and collaborator centers, now Accelerated
Breeding Initiative (ABI)\] to standardize the methods and key
performance indicators (KPIs) used across the CGIAR for decision making.
For our collaborators, please keep this in mind when submitting pull
requests for new methods since these may not be accepted if the are not
under the current recommended approaches (see contribution process
below).

The target audience is the breeders and geneticists from the CGIAR but
this is a package that runs behind the scenes of the bioflow interface.

## Roles and responsabilities

In this project we have two different roles

Maintainers: Have the permission and responsibility to commit, push, and
merge changes. Contributor: Have the permission and responsibility to
commit, push using forked versions.

The current list of maintainers and contributors is the following:

Mainteiners: Khaled Al-Sham’aa, Ibnou Dieng, Eduardo Covarrubias, Angela
Pacheco

Contributor: Bert De Boeck, Raul Eizaguirre, Abhishek Rathore, Roma Das,
Keith Gardner, Fernando Toledo, Juan Burgueno, Aubin Amagnide, Luis
Delgado Munoz, Christian Cadena, Christian Werner, Sergio Cruz, Alaine
Gulles, Justine Bonifacio.

## Installation

You can install the development version of cgiarPipeline like this:

``` r
# devtools::install_github("Breeding-Analytics/cgiarPipeline")
```

You can clone the GitHub repository in your computer if you would like
to add something or make a change as

``` r
# git clone https://github.com/Breeding-Analytics/cgiarPipeline.git
```

## Multi-environment analysis with sommer

`metSommer()` provides the list-based model interface used by `metASREML()`
without requiring an ASReml license. Install sommer >= 4.5.25 and this version
of cgiarPipeline before using the sommer tab in bioflow.

```r
# Main + environment-specific genetic variance
# fit <- metSommer(result, analysisId = soa_stamp, trait = "Yield_Mg_ha",
#   fixedTerm = list("environment"),
#   randomTerm = list("designation", c("designation", "environment")),
#   covMod = list("none", c("none", "Structure model_diag")))
```

The function supports weighted REML, Gaussian/binomial/Gamma/Poisson families
(non-Gaussian fits use PQL), FA and unstructured covariance, random regressions,
marker/pedigree/weather kernels, parental GCA/SCA, and configured TPP traits.
Results are appended under `moaSommer` and can be used by `indexDesire()`.
See `?metSommer` for weighting, prediction-scale and filtering conventions.

## Contribution process

In this collaborative development model, anyone can fork an existing
repository and push changes to their personal fork. You do not need
permission from the original repository to push to your own fork. The
changes can be merged into the original repository by the project
maintainer. When you open a pull request proposing changes from your
user-owned fork to a branch in the source (upstream) repository, you can
allow anyone with push access to the upstream repository to make changes
to your pull request. This model is popular with open-source projects as
it reduces the amount of friction for new contributors and allows people
to work independently without upfront coordination.

For more details, please check the following guidelines:
<https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project>

### Contributing an R function or pipeline?

If you wish to contribute a function that can be used across modules
(e.g., a function to summarize data) please make your pull request to
the cgiarBase repository.

If you wish to contribute a pipeline in the form of an R function, make
sure that it can use the data structure proposed (an example can be
found calling data(example) in any of the packages) and submit your new
pipeline function to the cgiarPipeline repository.

If you wish to contribute a shiny module please make sure that it can
use the data structure proposed and make your pull request to the
bioflow repository.

### How is a contribution review and accepted?

Priority issues and functionalities will be posted in the confluence
space for internal and external collaborators interested in contributing
(link).

The submitted functions will be tested using a sample of multiple
datasets collected by the different centers to ensure that new functions
and interfaces perform well across a variety of scenarios.

### Types of contributions that will be accepted (check list)

Only contributions that met the following criteria will be accepted:

1)  Functions that contribute to the main goal of bioflow (see scope
    section above).

2)  Functions that enable methods approved by the Quantitative Genetics
    community lead by the Accelerated Breeding Initiative (before EiB).

3)  Functions that are proven better than existing methods, robust
    enough to handle the variability of data encountered across the
    CGIAR (use sample datasets to test your own function).

### Open an issue before submitting a pull-request

For more information about GitHub Issues and Pull Request (PR)
templates, please check the following link:
<https://github.blog/2016-02-17-issue-and-pull-request-templates/>

### References

Bernardo Rex. 2010. Breeding for quantitative traits in plants. Second
edition. Stemma Press. 390 pp.

Gilmour et al. 1995. Average Information REML: An efficient algorithm
for variance parameter estimation in linear mixed models. Biometrics
51(4):1440-1450.

Kang et al. 2008. Efficient control of population structure in model
organism association mapping. Genetics 178:1709-1723.

Lee, D.-J., Durban, M., and Eilers, P.H.C. (2013). Efficient
two-dimensional smoothing with P-spline ANOVA mixed models and nested
bases. Computational Statistics and Data Analysis, 61, 22 - 37.

Lee et al. 2015. MTG2: An efficient algorithm for multivariate linear
mixed model analysis based on genomic information. Cold Spring Harbor.
doi: <http://dx.doi.org/10.1101/027201>.

Maier et al. 2015. Joint analysis of psychiatric disorders increases
accuracy of risk prediction for schizophrenia, bipolar disorder, and
major depressive disorder. Am J Hum Genet; 96(2):283-294.

Rodriguez-Alvarez, Maria Xose, et al. Correcting for spatial
heterogeneity in plant breeding experiments with P-splines. Spatial
Statistics 23 (2018): 52-71.

Searle. 1993. Applying the EM algorithm to calculating ML and REML
estimates of variance components. Paper invited for the 1993 American
Statistical Association Meeting, San Francisco.

Yu et al. 2006. A unified mixed-model method for association mapping
that accounts for multiple levels of relatedness. Genetics 38:203-208.

Tunnicliffe W. 1989. On the use of marginal likelihood in time series
model estimation. JRSS 51(1):15-27.

Zhang et al. 2010. Mixed linear model approach adapted for genome-wide
association studies. Nat. Genet. 42:355-360.

Jensen, J., Mantysaari, E. A., Madsen, P., and Thompson, R. (1997).
Residual maximum likelihood estimation of (co) variance components in
multivariate mixed linear models using average information. Journal of
the Indian Society of Agricultural Statistics, 49, 215-236.

Covarrubias-Pazaran G. Genome assisted prediction of quantitative traits
using the R package sommer. PLoS ONE 2016, 11(6):
<doi:10.1371/journal.pone.0156744>
