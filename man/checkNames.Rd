\name{checkNames}
\alias{checkNames}

\title{Check names of traits}
\description{
Check that trait names correspond with standard labels or CO numbers in 
\url{https://cropontology.org/}.
}
\usage{
checkNames(phenoDTfile, crop)
}
\arguments{
\item{crop}{A character string to identify the crop (\code{potato} or \code{sweetpotato}).}
\item{phenoDTfile}{The data object structure produced from bioflow which contains data, metadata, modifications and tables of results.}
}
\details{
The data frame must use the labels (lower or upper case) listed it \code{?st4gi::check.names.pt} and
\code{?st4gi::check.names.sp}.
}
\value{
It returns a data frame with all traits names in lower case, and a list of the
traits with names not recognized.
}

\examples{
####=========================================####
data(DT_example)

tmp <- checkNames(crop = 'potato', phenoDTfile = result)

}

