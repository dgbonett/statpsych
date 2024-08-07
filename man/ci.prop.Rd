% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych3.R
\name{ci.prop}
\alias{ci.prop}
\alias{ci.prop1}
\title{Confidence intervals for a proportion}
\usage{
ci.prop(alpha, f, n)

ci.prop1(alpha, f, n)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{f}{number of participants who have the attribute}

\item{n}{sample size}
}
\value{
Returns a 2-row matrix. The columns of row 1 are:
\itemize{
\item Estimate - adjusted estimate of proportion
\item SE - adjusted standard error
\item LL - lower limit of the adjusted Wald confidence interval
\item UL - upper limit of the adjusted Wald confidence interval
}

The columns of row 2 are:
\itemize{
\item Estimate - ML estimate of proportion
\item SE - standard error
\item LL - lower limit of the Wilson confidence interval
\item UL - upper limit of the Wilson confidence interval
}
}
\description{
Computes adjusted Wald and Wilson confidence intervals for a population
proportion. The Wilson confidence interval uses a continuity correction.

ci.prop1 is deprecated and will soon be removed from statpsych; please switch to ci.prop
}
\examples{
ci.prop(.05, 12, 100)

# Should return:
#                  Estimate         SE         LL        UL
# Adjusted Wald   0.1346154 0.03346842 0.06901848 0.2002123
# Wilson with cc  0.1200000 0.03249615 0.06625153 0.2039772


}
\references{
\insertRef{Agresti1998}{statpsych}
}
