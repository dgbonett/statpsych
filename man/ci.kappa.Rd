% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych3.R
\name{ci.kappa}
\alias{ci.kappa}
\title{Confidence interval for two kappa reliability coefficients}
\usage{
ci.kappa(alpha, f00, f01, f10, f11)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{f00}{number of objects rated y = 0 and x = 0}

\item{f01}{number of objects rated y = 0 and x = 1}

\item{f10}{number of objects rated y = 1 and x = 0}

\item{f11}{number of objects rated y = 1 and x = 1}
}
\value{
Returns a 2-row matrix. The results in row 1 are for the intraclass
kappa. The results in row 2 are for Cohen's kappa. The columns are:
\itemize{
\item Estimate - estimate of interrater reliability
\item SE - standard error
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes confidence intervals for the intraclass kappa coefficient and
Cohen's kappa coefficient with two dichotomous ratings.
}
\examples{
ci.kappa(.05, 31, 12, 4, 58)

# Should return:
#               Estimate         SE        LL        UL
# IC kappa:    0.6736597 0.07479965 0.5270551 0.8202643
# Cohen kappa: 0.6756757 0.07344761 0.5317210 0.8196303


}
\references{
\insertRef{Fleiss2003}{statpsych}
}
