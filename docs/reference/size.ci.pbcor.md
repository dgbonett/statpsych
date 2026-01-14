# Sample size for a point-biserial correlation confidence interval

Computes the sample size required to estimate a population
point-biserial correlation with desired confidence interval precision in
a two-group nonexperimental design with simple random sampling. A
two-group nonexperimental design implies two subpopulations (e.g., all
boys and all girls in a school district). This function requires a
planning value for the proportion of population members who belong to
one of the two subpopulations. Set the correlation planning value to the
smallest absolute value within a plausible range for a conservatively
large sample size.

For more details, see Section 1.24 of Bonett (2021, Volume 2)

## Usage

``` r
size.ci.pbcor(alpha, cor, w, p)
```

## Arguments

- alpha:

  alpha level for 1-alpha confidence

- cor:

  planning value of point-biserial correlation

- w:

  desired confidence interval width

- p:

  proportion of members in one of the two subpopulations

## Value

Returns the required sample size

## References

Bonett DG (2020). “Point-biserial correlation: Interval estimation,
hypothesis testing, meta-analysis, and sample size determination.”
*British Journal of Mathematical and Statistical Psychology*,
**73**(S1), 113–144. ISSN 0007-1102,
[doi:10.1111/bmsp.12189](https://doi.org/10.1111/bmsp.12189) .

Bonett DG (2021). *Statistical Methods for Psychologists
https://dgbonett.sites.ucsc.edu/*.

## Examples

``` r
size.ci.pbcor(.05, .40, .25, .73)
#>  Sample size
#>          168

# Should return:
# Sample size
#         168
 
```
