% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{size.ci.lc.stdmean.ws}
\alias{size.ci.lc.stdmean.ws}
\title{Sample size for a within-subjects standardized linear contrast of means
confidence interval}
\usage{
size.ci.lc.stdmean.ws(alpha, d, cor, w, q)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{d}{planning value of standardized linear contrast}

\item{cor}{planning value of average correlation between measurements}

\item{w}{desired confidence interval width}

\item{q}{vector of within-subjects contrast coefficients}
}
\value{
Returns the required sample size for each standardizer
}
\description{
Computes the sample size required to estimate two types of standardized
linear contrasts of population means (unweighted standardizer and single
level standardizer) with desired confidence interval precision in a
within-subjects design. For a conservatively large sample size, set the
standardized linear contrast of means planning value to the largest value
within a plausible range, and set the Pearson correlation planning value
to the smallest value within a plausible range.
}
\examples{
q <- c(.5, .5, -.5, -.5)
size.ci.lc.stdmean.ws(.05, 1, .7, .6, q)

# Should return:
#                            Sample size
# Unweighted standardizer:            26
# Single level standardizer:          35
 

}
\references{
\insertRef{Bonett2009}{statpsych}
}
