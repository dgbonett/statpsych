% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych1.R
\name{size.equiv.mean.ps}
\alias{size.equiv.mean.ps}
\title{Sample size for a paired-samples mean equivalence test}
\usage{
size.equiv.mean.ps(alpha, pow, var, es, cor, h)
}
\arguments{
\item{alpha}{alpha level for hypothesis test}

\item{pow}{desired power}

\item{var}{planning value of average variance of the two measurements}

\item{es}{planning value of mean difference}

\item{cor}{planning value of the correlation between measurements}

\item{h}{upper limit for range of practical equivalence}
}
\value{
Returns the required sample size
}
\description{
Computes the sample size required to perform an equivalence test for the
difference in population means with desired power in a paired-samples
design. The value of h specifies a range of practical equivalence, -h to h,
for the difference in population means. The planning value for the absolute
mean difference must be less than h. Equivalence tests often require a
very large sample size. Equivalence tests usually use 2 x alpha rather than
alpha (e.g., use alpha = .10 rather alpha = .05). Set the Pearson correlation
value to the smallest value within a plausible range, and set the variance
planning value to the largest value within a plausible range for a
conservatively large sample size.
}
\examples{
size.equiv.mean.ps(.10, .85, 15, .5, .7, 1.5)

# Should return:
# Sample size
#          68
 

}
