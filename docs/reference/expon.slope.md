# Confidence interval for an exponentiated slope

Computes confidence intervals for exp(B) - 1 (as a percent) and exp(B)
where B is a population slope coefficient in a binary logit, ordinal
logit, or log-Poisson model. This function is useful with software that
does not have an option to compute exp(B) and exp(B) - 1.

For more details, see Section 4.1 of Bonett (2021, Volume 3)

## Usage

``` r
expon.slope(alpha, b, se)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- b:

  estimated slope coefficient

- se:

  slope standard error

## Value

Returns a 2-row matrix. The first row gives the results for exp(B), and
the the second row gives the results for exp(B) - 1 (as a percent). The
columns are:

- Estimate - estimate of exp(B) or exp(B) - 1

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
expon.slope(.05, .502, .0396)
#>                   Estimate        LL       UL
#> exp(B)            1.652022  1.528651  1.78535
#> 100[exp(B) - 1]% 65.202201 52.865066 78.53502

# Should return:
#                   Estimate        LL       UL
# exp(B)            1.652022  1.528651  1.78535
# 100[exp(B) - 1]% 65.202201 52.865066 78.53502

```
