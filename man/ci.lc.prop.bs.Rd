% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych3.R
\name{ci.lc.prop.bs}
\alias{ci.lc.prop.bs}
\title{Confidence interval for a linear contrast of proportions in a between-
subjects design}
\usage{
ci.lc.prop.bs(alpha, f, n, v)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{f}{vector of frequency counts of participants who have the attribute}

\item{n}{vector of sample sizes}

\item{v}{vector of between-subjects contrast coefficients}
}
\value{
Returns a 1-row matrix. The columns are:
\itemize{
\item Estimate - adjusted estimate of proportion linear contrast
\item SE - adjusted standard error
\item z - z test statistic
\item p - two-sided p-value
\item LL - lower limit of the adjusted Wald confidence interval
\item UL - upper limit of the adjusted Wald confidence interval
}
}
\description{
Computes an adjusted Wald confidence interval for a linear contrast of
population proportions in a between-subjects design.
}
\examples{
f <- c(26, 24, 38)
n <- c(60, 60, 60)
v <- c(-.5, -.5, 1)
ci.lc.prop.bs(.05, f, n, v)

# Should return:
#  Estimate         SE        z           p         LL        UL
# 0.2119565 0.07602892 2.787841 0.005306059 0.06294259 0.3609705


}
\references{
\insertRef{Price2004}{statpsych}
}
