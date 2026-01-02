# Hypothesis test for a mean

Computes a one-sample t-test for a population mean using the estimated
mean, estimated standard deviation, sample size, and null hypothesis
value. Use the t.test function for raw data input. A confidence interval
for a population mean is a recommended supplement to the t-test (see
[ci.mean](https://dgbonett.github.io/statpsych/reference/ci.mean.md)).

For more details, see Section 1.11 of Bonett (2021, Volume 1)

## Usage

``` r
test.mean(m, sd, n, h)
```

## Arguments

- m:

  estimated mean

- sd:

  estimated standard deviation

- n:

  sample size

- h:

  null hypothesis value of mean

## Value

Returns a 1-row matrix. The columns are:

- t - t test statistic

- df - degrees of freedom

- p - two-sided p-value

## References

Snedecor GW, Cochran WG (1989). *Statistical Methods*, 8th edition. ISU
University Pres, Ames, Iowa. Bonett DG (2021url). *Statistical Methods
for Psychologists*.

## Examples

``` r
test.mean(7.9, 3.05, 100, 7)
#>       t df       p
#>  2.9508 99 0.00396

# Should return:
#        t df       p
#   2.9508 99 0.00396
 
```
