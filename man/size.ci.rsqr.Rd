% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych2.R
\name{size.ci.rsqr}
\alias{size.ci.rsqr}
\title{Sample size for a squared multiple correlation confidence interval}
\usage{
size.ci.rsqr(alpha, r2, s, w)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{r2}{planning value of squared multiple correlation}

\item{s}{number of predictor variables in model}

\item{w}{desired confidence interval width}
}
\value{
Returns the required sample size
}
\description{
Computes the sample size required to estimate a population squared multiple
correlation in a random-x regression model with desired confidence interval
precision. Set the planning value of the squared multiple correlation to 1/3
for a conservatively large sample size.
}
\examples{
size.ci.rsqr(.05, .25, 5, .2)

# Should return:
# Sample size
#         214
 

}
