% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych3.R
\name{ci.yule}
\alias{ci.yule}
\title{Confidence intervals for generalized Yule coefficients}
\usage{
ci.yule(alpha, f00, f01, f10, f11)
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
\item Estimate - estimate of generalized Yule coefficient
\item SE - standard error
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes confidence intervals for four generalized Yule measures of
association (Yule Q, Yule Y, Digby H, and Bonett-Price Y*) using a
transformation of a confidence interval for an odds ratio with .5 added to
each cell frequency. This function requires the frequency counts from a
2 x 2 contingency table for two dichotomous variables. Digby H is sometimes
used as a crude approximation to the tetrachoric correlation. Yule Y is
equal to the phi coefficient only when all marginal frequencies are equal.
Bonett-Price Y* is a better approximation to the phi coefficient when the
marginal frequencies are not equal.
}
\examples{
ci.yule(.05, 229, 28, 96, 24)

# Should return:
#      Estimate         SE         LL        UL
# Q:  0.3430670 0.13280379 0.06247099 0.5734020
# Y:  0.1769015 0.07290438 0.03126603 0.3151817
# H:  0.2619244 0.10514465 0.04687994 0.4537659
# Y*: 0.1311480 0.05457236 0.02307188 0.2361941


}
\references{
\insertRef{Bonett2007}{statpsych}
}
