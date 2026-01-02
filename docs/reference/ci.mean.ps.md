# Confidence interval for a paired-samples mean difference

Computes a confidence interval for a population paired-samples mean
difference using the estimated means, estimated standard deviations,
estimated correlation, and sample size. Also computes a paired-samples
t-test. Use the t.test function for raw data input.

For more details, see Section 4.2 of Bonett (2021, Volume 1)

## Usage

``` r
ci.mean.ps(alpha, m1, m2, sd1, sd2, cor, n)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- m1:

  estimated mean for measurement 1

- m2:

  estimated mean for measurement 2

- sd1:

  estimated standard deviation for measurement 1

- sd2:

  estimated standard deviation for measurement 2

- cor:

  estimated correlation between measurements

- n:

  sample size

## Value

Returns a 1-row matrix. The columns are:

- Estimate - estimated mean difference

- SE - standard error

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

- LL - lower limit of the confidence interval

- UL - upper limit of the confidence interval

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa. Bonett DG (2021url). *Statistical Methods
for Psychologists*.

## Examples

``` r
ci.mean.ps(.05, 58.2, 51.4, 7.43, 8.92, .537, 30)
#>  Estimate       SE      t df     p       LL       UL
#>       6.8 1.455922 4.6706 29 6e-05 3.822304 9.777696

# Should return:
# Estimate       SE      t df      p       LL       UL
#      6.8 1.455922 4.6706 29  6e-05 3.822304 9.777696

```
