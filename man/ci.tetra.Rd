% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych3.R
\name{ci.tetra}
\alias{ci.tetra}
\title{Confidence interval for a tetrachoric correlation}
\usage{
ci.tetra(alpha, f00, f01, f10, f11)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{f00}{number of participants with y = 0 and x = 0}

\item{f01}{number of participants with y = 0 and x = 1}

\item{f10}{number of participants with y = 1 and x = 0}

\item{f11}{number of participants with y = 1 and x = 1}
}
\value{
Returns a 1-row matrix. The columns are:
\itemize{
\item Estimate - estimate of tetrachoric approximation
\item SE - recovered standard error
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes a confidence interval for an approximation to the tetrachoric
correlation. This function requires the frequency counts from a 2 x 2
contingency table for two dichotomous variables. This measure of
association assumes both of the dichotomous variables are artificially
dichotomous. An approximate standard error is recovered from the
confidence interval.
}
\examples{
ci.tetra(.05, 46, 15, 54, 85)

# Should return:
#  Estimate         SE        LL        UL
# 0.5135167 0.09301703 0.3102345 0.6748546


}
\references{
\insertRef{Bonett2005}{statpsych}
}
