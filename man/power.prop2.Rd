% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych3.R
\name{power.prop2}
\alias{power.prop2}
\title{Approximates the power of a 2-group proportion test for planned sample
sizes}
\usage{
power.prop2(alpha, n1, n2, p1, p2, es)
}
\arguments{
\item{alpha}{alpha level for hypothesis test}

\item{n1}{planned sample size for group 1}

\item{n2}{planned sample size for group 2}

\item{p1}{planning value of proportion for group 1}

\item{p2}{planning value of proportion for group 2}

\item{es}{planning value of proportion difference}
}
\value{
Returns the approximate power of the test
}
\description{
Computes the approximate power for a test of equal population proportions
in a 2-group design for the planned sample sizes. This function requires
planning values for both proportions. Set the proportion planning values
to .5 for a conservatively low power estimate. The planning value for the
proportion difference could be set to the difference of the two proportion
planning values or it could be set to a minimally interesting effect size.
}
\examples{
power.prop2(.05, 60, 40, .5, .5, .2)

# Should return:
#     Power
# 0.4998959


}
