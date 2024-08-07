% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{ci.cv}
\alias{ci.cv}
\title{Confidence interval for a coefficient of variation}
\usage{
ci.cv(alpha, m, sd, n)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{m}{estimated mean}

\item{sd}{estimated standard deviation}

\item{n}{sample size}
}
\value{
Returns a 1-row matrix. The columns are:
\itemize{
\item Estimate - estimated coefficient of variation
\item SE - recovered standard error
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes a confidence interval for a population coefficient of variation
(standard deviation divided by mean). This confidence interval is the
reciprocal of a confidence interval for a standardized mean (see
\link[statpsych]{ci.stdmean}). An approximate standard error is recovered
from the confidence interval. The coefficient of variation assumes
ratio-scale scores.
}
\examples{
ci.cv(.05, 24.5, 3.65, 40)

# Should return:
#  Estimate        SE        LL       UL
# 0.1489796 0.01817373 0.1214381 0.1926778
 

}
