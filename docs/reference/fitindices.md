# SEM fit indices

Computes the normed fit index (NFI), adjusted normed fit index (adj
NFI), comparative fit index (CFI), Tucker-Lewis fit index (TLI), and
root mean square error of approximation index (RMSEA). Of the first four
indices, the adj NFI index is recommended because it has smaller
sampling variability than CFI and TLI and less negative bias than NFI.

For more details, see Section 2.14 of Bonett (2021, Volume 4)

## Usage

``` r
fitindices(chi1, df1, chi2, df2, n)
```

## Arguments

- chi1:

  chi-square test statistic for full model

- df1:

  degrees of freedom for full model

- chi2:

  chi-square test statistic for reduced model

- df2:

  degrees of freedom for reduced model

- n:

  sample size

## Value

Returns NFI, adj NFI, CFI, TLI, and RMSEA

## References

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
fitindices(14.21, 10, 258.43, 20, 300)
#>    NFI adj NFI  CFI    TLI  RMSEA
#>  0.945  0.9837 0.98 0.9647 0.0375

# Should return:
#    NFI   adj NFI   CFI    TLI  RMSEA
#  0.945    0.9837  0.98 0.9647 0.0375
 
```
