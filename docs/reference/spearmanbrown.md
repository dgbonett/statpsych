# Computes the reliability of a scale with r2 measurements given the reliability of a scale with r1 measurements

Computes the reliability of a scale that is the sum or average of r2
parallel measurements given the reliability of a scale that is the sum
or average of r1 parallel measurements. The "measurements" can be items,
forms, raters, or occasions.

For more details, see Section 4.19 of Bonett (2021, Volume 1)

## Usage

``` r
spearmanbrown(rel, r1, r2)
```

## Arguments

- rel:

  reliability of the sum or average of r1 measurements

- r1:

  number of measurements in the original scale

- r2:

  number of measurements in the new scale

## Value

Returns the reliability of the sum or average of r2 measurements

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
spearmanbrown(.6, 10, 20)
#>  Reliability of r2 measurements
#>                            0.75

# Should return:
# Reliability of r2 measurements
#                            .75

```
