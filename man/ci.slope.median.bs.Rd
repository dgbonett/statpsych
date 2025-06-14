% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych2.R
\name{ci.slope.median.bs}
\alias{ci.slope.median.bs}
\title{Confidence interval for the slope of medians in a one-factor experimental
design with a quantitative between-subjects factor}
\usage{
ci.slope.median.bs(alpha, m, se, x)
}
\arguments{
\item{alpha}{alpha level for 1-alpha confidence}

\item{m}{vector of sample median}

\item{se}{vector of standard errors}

\item{x}{vector of quantitative factor values}
}
\value{
Returns a 1-row matrix. The columns are:
\itemize{
\item Estimate - estimated slope
\item SE - standard error
\item z - z test statistic
\item p - two-sided p-value
\item LL - lower limit of the confidence interval
\item UL - upper limit of the confidence interval
}
}
\description{
Computes a distrbution-free test and confidence interval for the slope
of medians in a one-factor experimental design with a quantitative
between-subjects factor using sample group medians and standard errors
as input. The sample median and standard error for each group can be
computed using the \link[statpsych]{ci.median} function.
}
\examples{
m <- c(33.5, 37.9, 38.0, 44.1)
se <- c(0.84, 0.94, 1.65, 2.98)
x <- c(5, 10, 20, 30)
ci.slope.median.bs(.05, m, se, x)

# Should return:
#   Estimate        SE        z           p        LL        UL
#  0.3664407 0.1163593 3.149216 0.001637091 0.1383806 0.5945008


}
