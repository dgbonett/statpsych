# Confidence interval for a semipartial correlation

Computes a Fisher confidence interval for a population semipartial
correlation. This function requires an (unadjusted) estimate of the
squared multiple correlation in the full model that contains the
predictor variable of interest plus all control variables. This function
computes a modified Aloe-Becker confidence interval that uses n - 3
rather than n in the standard error and also uses a Fisher
transformation of the semipartial correlation.

For more details, see Section 2.7 of Bonett (2021, Volume 2)

## Usage

``` r
ci.spcor(alpha, cor, r2, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  estimated semipartial correlation

- r2:

  estimated squared multiple correlation in full model

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated semipartial correlation (from input)

- SE - standard error

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

Aloe AM, Becker BJ (2012). “An effect size for regression predictors in
meta-analysis.” *Journal of Educational and Behavioral Statistics*,
**37**(2), 278–297. ISSN 1076-9986,
[doi:10.3102/1076998610396901](https://doi.org/10.3102/1076998610396901)
.

## Examples

``` r
ci.spcor(.05, -.355, .230, 236)
#>  Estimate      SE      LL      UL
#>    -0.355 0.05426 -0.4565 -0.2444

# Should return:
# Estimate      SE      LL      UL
#   -0.355 0.05426 -0.4565 -0.2444
 
```
