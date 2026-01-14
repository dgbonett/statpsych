# Sample size for a second-stage confidence interval

Computes the second-stage sample size required to obtain desired
confidence interval precision. This function can use either the total
sample size for all groups in the first stage sample or a single group
sample size in the first stage sample. If the total first-stage sample
size is given, then the function computes the total sample size required
in the second-stage sample. If a single group first-stage sample size is
given, then the function computes the single-group sample size required
in the second-stage sample. The second-stage sample is combined with the
first-stage sample to obtain the desired confidence interval width.

For more details, see Section 1.30 of Bonett (2021, Volume 1)

## Usage

``` r
size.ci.second(n0, w0, w)
```

## Arguments

- n0:

  first-stage sample size

- w0:

  confidence interval width in first-stage sample

- w:

  desired confidence interval width

## Value

Returns the required sample size for the second-stage sample

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.second(25, 4.38, 2.0)
#>  Second-stage sample size
#>                        95

# Should return:
# Second-stage sample size
#                       95
 
```
