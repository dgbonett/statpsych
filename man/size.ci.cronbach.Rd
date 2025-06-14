% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{size.ci.cronbach}
\alias{size.ci.cronbach}
\title{Sample size for a Cronbach reliability confidence interval}
\usage{
size.ci.cronbach(alpha, rel, r, w)
}
\arguments{
\item{alpha}{alpha value for 1-alpha confidence}

\item{rel}{reliability planning value}

\item{r}{number of measurements (items, raters, forms)}

\item{w}{desired confidence interval width}
}
\value{
Returns the required sample size
}
\description{
Computes the sample size required to estimate a Cronbach reliability
with desired confidence interval precision. Set the reliability planning
value to the smallest value within a plausible range for a
conservatively large sample size.
}
\examples{
size.ci.cronbach(.05, .85, 5, .1)

# Should return:
# Sample size
#          89
 

}
\references{
\insertRef{Bonett2015}{statpsych}
}
