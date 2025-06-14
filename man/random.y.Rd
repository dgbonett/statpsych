% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{random.y}
\alias{random.y}
\title{Generate random sample of scores}
\usage{
random.y(n, m, sd, min, max, dec)
}
\arguments{
\item{n}{sample size}

\item{m}{population mean of scores}

\item{sd}{population standard deviation of scores}

\item{min}{minimum allowable value}

\item{max}{maximum allowable value}

\item{dec}{number of decimal points}
}
\value{
Returns a vector of randomly generated scores.
}
\description{
Generates a random sample of scores from a normal distribution with a
specified population mean and standard deviation. This function is useful
for generating hypothetical data for classroom demonstrations.
}
\examples{
random.y(10, 3.6, 2.8, 1, 7, 0) 

# Should return random numbers such as:
# [1] 2 7 7 1 6 3 1 3 2 1
 

}
