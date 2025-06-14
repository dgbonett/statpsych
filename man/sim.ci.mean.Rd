% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{sim.ci.mean}
\alias{sim.ci.mean}
\title{Simulates confidence interval coverage probability for a mean}
\usage{
sim.ci.mean(alpha, n, dist, rep)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{n}{sample size}

\item{dist}{type of distribution (1, 2, 3, 4,or 5)
\itemize{
\item 1 = Gaussian (skewness = 0 and excess kurtosis = 0)
\item 2 = platykurtic (skewness = 0 and excess kurtosis = -1.2)
\item 3 = leptokurtic (skewness = 0 and excess kurtosis = 6)
\item 4 = moderate skew (skewness = 1 and excess kurtosis = 1.5)
\item 5 = large skew (skewness = 2 and excess kurtosis = 6)
}}

\item{rep}{number of Monte Carlo samples}
}
\value{
Returns a 1-row matrix. The columns are:
\itemize{
\item Coverage - probability of confidence interval including population mean
\item Lower Error - probability of lower limit greater than population mean
\item Upper Error - probability of upper limit less than population mean
\item Ave CI Width - average confidence interval width
}
}
\description{
Performs a computer simulation of the confidence interval performance for a
population mean. Sample data can be generated from five different population
distributions. All distributions are scaled to have a standard deviation
of 1.0.
}
\examples{
sim.ci.mean(.05, 10, 1, 5000)

# Should return (within sampling error):
# Coverage Lower Error Upper Error Ave CI Width
#   0.9484      0.0264      0.0252     1.392041

sim.ci.mean(.05, 40, 4, 1000)

# Should return (within sampling error):
# Coverage Lower Error Upper Error Ave CI Width
#  0.94722     0.01738      0.0354    0.6333067


}
