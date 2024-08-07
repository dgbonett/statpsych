% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{power.mean.ps}
\alias{power.mean.ps}
\title{Approximates the power of a paired-samples t-test for a planned sample size}
\usage{
power.mean.ps(alpha, n, var1, var2, es, cor)
}
\arguments{
\item{alpha}{alpha level for hypothesis test}

\item{n}{planned sample size}

\item{var1}{planning value of measurement 1 variance}

\item{var2}{planning value of measurement 2 variance}

\item{es}{planning value of mean difference}

\item{cor}{planning value of correlation between measurements}
}
\value{
Returns the approximate power of the test
}
\description{
Computes the approximate power of a paired-samples t-test for a planned
sample size. For a conservatively low power approximation, set the variance
planning values to the largest values within their plausible ranges, set the
correlation planning value to the smallest value within its plausible range,
and set the effect size to a minimally interesting value. The variances of
the two measurements can be unequal.
}
\examples{
power.mean.ps(.05, 20, 10.0, 12.0, 2, .7)

# Should return:
#     Power
# 0.9074354


}
