# Sample size for an indirect effect confidence interval

Computes the approximate sample size required to estimate a population
standardized indirect effect in a simple mediation model. The direct
effect of the independent (exogenous) variable on the response variable,
controlling for the mediator variable, is assumed to be negligible.

For more details, see Section 1.19 of Bonett (2021, Volume 4)

## Usage

``` r
size.ci.indirect(alpha, cor1, cor2, w)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor1:

  planning value of correlation between the independent and mediator
  variables

- cor2:

  planning value of correlation between the mediator and response
  variables

- w:

  desired confidence interval width

## Value

Returns the required sample size

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.indirect(.05, .4, .5, .2)
#>  Sample size
#>          106

# Should return:
# Sample size
#         106
 
```
