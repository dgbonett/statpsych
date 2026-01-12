# Prediction interval for a difference of scores in a 2-group experiment

For a 2-group experimental design, this function computes a prediction
interval for how the response variable score for one randomly selected
person from the study population would differ under the two treatment
conditions. Both equal variance and unequal variance prediction
intervals are computed.

For more details, see Section 2.6 of Bonett (2021, Volume 1)

## Usage

``` r
pi.score2(alpha, m1, m2, sd1, sd2, n1, n2)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estamated mean for group 1

- m2:

  estimated mean for group 1

- sd1:

  estimated standard deviation for group 1

- sd2:

  estimated standard deviation for group 2

- n1:

  sample size for group 1

- n2:

  sample size for group 2

## Value

Returns a 2-row matrix. The columns are:

- Predicted - predicted difference in scores

- df - degrees of freedom

- LL - lower limit of the prediction interval

- UL - upper limit of the prediction interval

## References

Hahn GJ (1977). “A prediction interval on the difference between two
future sample means and its application to a claim of product
superiority.” *Technometrics*, **19**(2), 131–134. ISSN 0040-1706,
[doi:10.1080/00401706.1977.10489520](https://doi.org/10.1080/00401706.1977.10489520)
.

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
pi.score2(.05, 19.4, 11.3, 2.70, 2.10, 40, 40)
#>                              Predicted    df       LL       UL
#> Equal Variances Assumed:           8.1 78.00 1.205659 14.99434
#> Equal Variances Not Assumed:       8.1 73.54 1.199067 15.00093

# Should return:
#                              Predicted    df       LL       UL
# Equal Variances Assumed:           8.1 78.00 1.205659 14.99434
# Equal Variances Not Assumed:       8.1 73.54 1.199073 15.00093
 
```
