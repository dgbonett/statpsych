% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{ci.etasqr}
\alias{ci.etasqr}
\title{Confidence interval for eta-squared}
\usage{
ci.etasqr(alpha, etasqr, df1, df2)
}
\arguments{
\item{alpha}{alpha value for 1-alpha confidence}

\item{etasqr}{estimated eta-squared}

\item{df1}{degrees of freedom for effect}

\item{df2}{error degrees of freedom}
}
\value{
Returns a 1-row matrix. The columns are:
\itemize{
\item Eta-squared - estimate of eta-squared
\item adj Eta-squared - bias adjusted eta-squared estimate
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes a confidence interval for a population eta-squared, partial
eta-squared, or generalized eta-squared in a fixed-factor between-subjects
design. An approximate bias adjusted estimate is also computed.
}
\examples{
ci.etasqr(.05, .241, 3, 116)

# Should return:
#      Eta-squared adj Eta-squared        LL        UL
# [1,]       0.241       0.2213707 0.1040229 0.3493431
 

}