% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statpsych2.R
\name{size.test.lc.ancova}
\alias{size.test.lc.ancova}
\title{Sample size for a mean linear contrast test in an ANCOVA}
\usage{
size.test.lc.ancova(alpha, pow, evar, es, s, d, v)
}
\arguments{
\item{alpha}{alpha level for hypothesis test}

\item{pow}{desired power}

\item{evar}{planning value of within-group (error) variance}

\item{es}{planning value of linear contrast}

\item{s}{number of covariates}

\item{d}{largest standardized mean difference for all covariates}

\item{v}{vector of between-subjects contrast coefficients}
}
\value{
Returns the required sample size for each group
}
\description{
Computes the sample size for each group (assuming equal sample sizes) required
to test a linear contrast of population means in an ANCOVA model with desired
power. In a nonexperimental design, the sample size is affected by the magnitude
of covariate mean differences across groups. The covariate mean differences can
be approximated by specifying the largest standardized covariate mean difference
across all pairwise comparisons and for all covariates. In an experiment, this
standardized mean difference is set to 0. Set the error variance planning
value to the largest value within a plausible range for a conservatively
large sample size.
}
\examples{
v <- c(.5, .5, -1)
size.test.lc.ancova(.05, .9, 1.37, .7, 1, 0, v)

# Should return:
# Sample size per group
#                    46
 

}
