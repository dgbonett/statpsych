# Confidence interval for a reliability coefficient

Computes a confidence interval for a population reliability coefficient
such as Cronbach's alpha or McDonald's omega using an estimate of the
reliability and its standard error. The standard error can be a robust
standard error or bootstrap standard error obtained from an SEM program.
Use
[ci.cronbach](https://dgbonett.github.io/statpsych/reference/ci.cronbach.md)
for Cronbach's alpha if parallel measurements can be assumed.

For more details, see Section 2.8 of Bonett (2021, Volume 4)

## Usage

``` r
ci.reliability(alpha, rel, se, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- rel:

  estimated reliability

- se:

  standard error of reliability

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated reliability (from input)

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
ci.reliability(.05, .88, .0147, 100)
#>  Estimate    LL     UL
#>      0.88 0.849 0.9066

# Should return:
# Estimate     LL     UL
#     0.88  0.849 0.9066
 
```
