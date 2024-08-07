% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{ci.stdmean.strat}
\alias{ci.stdmean.strat}
\title{Confidence intervals for a 2-group standardized mean difference with
stratified sampling}
\usage{
ci.stdmean.strat(alpha, m1, m2, sd1, sd2, n1, n2, p1)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{m1}{estimated mean for group 1}

\item{m2}{estimated mean for group 2}

\item{sd1}{estimated standard deviation for group 1}

\item{sd2}{estimated standard deviation for group 2}

\item{n1}{sample size for group 1}

\item{n2}{sample size for group 2}

\item{p1}{proportion of total population in subpopulation 1}
}
\value{
Returns a 3-row matrix. The columns are:
\itemize{
\item Estimate - estimated standardized mean difference
\item adj Estimate - bias adjusted standardized mean difference estimate
\item SE - standard error
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes confidence intervals for a population standardized mean difference
in a 2-group nonexperimental design with stratified random sampling (a random
sample of a specified size from each subpopulation) using a square root
weighted variance standardizer or single group standard deviation
standardizer.  Equality of variances is not assumed.
}
\examples{
ci.stdmean.strat(.05, 33.2, 30.8, 10.5, 11.2, 200, 200, .533)

# Should return:
#                         Estimate  adj Estimate         SE         LL        UL
# Weighted standardizer: 0.2215549     0.2211371 0.10052057 0.02453817 0.4185716
# Group 1 standardizer:  0.2285714     0.2277089 0.10427785 0.02419059 0.4329523
# Group 2 standardizer:  0.2142857     0.2277089 0.09776049 0.02267868 0.4058927


}
\references{
\insertRef{Bonett2020a}{statpsych}
}
