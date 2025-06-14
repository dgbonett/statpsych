% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych3.R
\name{ci.cramer}
\alias{ci.cramer}
\title{Confidence interval for Cramer's V}
\usage{
ci.cramer(alpha, chisqr, r, c, n)
}
\arguments{
\item{alpha}{alpha value for 1-alpha confidence}

\item{chisqr}{Pearson chi-square test statistic of independence}

\item{r}{number of rows in contingency table}

\item{c}{number of columns in contingency table}

\item{n}{sample size}
}
\value{
Returns a 1-row matrix. The columns are:
\itemize{
\item Estimate - estimate of Cramer's V
\item SE - recovered standard error
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes a confidence interval for a population Cramer's V coefficient
of nominal association for an r x s contingency table. The confidence interval
is based on a noncentral chi-square distribution, and an approximate standard
error is recovered from the confidence interval.
}
\examples{
ci.cramer(.05, 19.21, 2, 3, 200)

# Should return:
# Estimate     SE     LL     UL
#   0.3099 0.0718 0.1601 0.4417
 

}
\references{
\insertRef{Smithson2003}{statpsych}
}
