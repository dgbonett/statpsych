# Sample size for a Mann-Whitney test

Computes the sample size in each group (assuming equal sample sizes)
required for the Mann-Whitney test with desired power. A planning value
of the Mann-Whitney parameter is required. In a 2-group experiment, this
parameter is the proportion of members in the population with scores
that would be larger under treatment 1 than treatment 2. In a 2-group
nonexperiment where participants are sampled from two subpopulations of
sizes N1 and N2, the parameter is the proportion of all N1 x N2 pairs in
which a member from subpopulation 1 has a larger score than a member
from subpopulation 2.

For more details, see Section 2.14 of Bonett (2021, Volume 1)

## Usage

``` r
size.test.mann(alpha, pow, p)
```

## Arguments

- alpha:

  alpha level for hypothesis test

- pow:

  desired power

- p:

  planning value of Mann-Whitney parameter

## Value

Returns the required sample size for each group

## References

Noether GE (1987). “Sample size determination for some common
nonparametric tests.” *Journal of the American Statistical Association*,
**82**(398), 645–647. ISSN 0162-1459,
[doi:10.1080/01621459.1987.10478478](https://doi.org/10.1080/01621459.1987.10478478)
. Bonett DG (2021url). *Statistical Methods for Psychologists*.

## Examples

``` r
size.test.mann(.05, .90, .3)
#>  Sample size per group
#>                     44

# Should return:
# Sample size per group
#                    44
 
```
