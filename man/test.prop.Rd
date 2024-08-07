% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych3.R
\name{test.prop}
\alias{test.prop}
\title{Hypothesis test for a proportion}
\usage{
test.prop(f, n, h)
}
\arguments{
\item{f}{number of participants who have the attribute}

\item{n}{sample size}

\item{h}{null hypothesis value of proportion}
}
\value{
Returns a 1-row matrix. The columns are:
\itemize{
\item Estimate - ML estimate of proportion
\item z - z test statistic
\item p - two-sided p-value
}
}
\description{
Computes a continuity-corrected z-test for a population proportion in a
1-group design. A confidence interval for a population proportion
is a recommended supplement to the z-test (see \link[statpsych]{ci.prop}).
}
\examples{
test.prop(9, 20, .2)

# Should return:
# Estimate        z          p
#     0.45 2.515576 0.01188379


}
\references{
\insertRef{Snedecor1980}{statpsych}
}
