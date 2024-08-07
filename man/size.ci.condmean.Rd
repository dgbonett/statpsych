% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych2.R
\name{size.ci.condmean}
\alias{size.ci.condmean}
\title{Sample size for a conditional mean confidence interval}
\usage{
size.ci.condmean(alpha, evar, xvar, diff, w)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{evar}{planning value of within group (error) variance}

\item{xvar}{variance of fixed predictor variable}

\item{diff}{difference between x* and mean of x}

\item{w}{desired confidence interval width}
}
\value{
Returns the required total sample size
}
\description{
Computes the total sample size required to estimate a population conditional
mean of y at x = x* in a fixed-x simple linear regression model with desired
confidence interval precision. In an experimental design, the total sample
size would be allocated to the levels of the quantitative factor and it
might be necessary to increase the total sample size to achieve equal
sample sizes. Set the error variance planning value to the largest value
within a plausible range for a conservatively large sample size.
}
\examples{
size.ci.condmean(.05, 120, 125, 15, 5)

# Should return:
# Total sample size
#               210
 

}
