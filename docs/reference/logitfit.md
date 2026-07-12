# The confusion matrix is computed using the observed 0 or 1 response variable scores, the predicted probabilities from the logistic model, and a specified cutpoint.

The confusion matrix is computed using the observed 0 or 1 response
variable scores, the predicted probabilities from the logistic model,
and a specified cutpoint.

## Usage

``` r
logitfit(y, p, cut)
```

## Arguments

- y:

  vector of observed 0 or 1 response scores

- p:

  vector of predicted probabilities from model

- cut:

  cutpoint (defines the predicted 0 or 1 scores)

## Value

Returns a 1-row matrix. The columns are:

- C - percent correctly classified (aka accuracy)

- TP - percent true positives (aka sensitivity, recall)

- FP - percent false positives (1 - TN)

- TN - percent true negatives (aka specificity)

- FN - percent false negatives (1 - TP)

- PPV - percent positive predicted values (aka precision)

- NPV - percent negative predicted values

- F1 - F1 score (2 x PPV x TP)/(PPV + TP)

## Examples

``` r
y <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
x1 <- c(1,1,3,2,6,8,4,5,6,2,4,3,1,5,3,9,8,9,8,6,6,7,5,3,8,6,5,7,8,9,7,8)
x2 <- c(0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0)
out <- glm(y ~ x1 + x2, family = binomial(link = "logit"))
p <- predict(out, type = "response")
logitfit(y, p, .3)
#>      C    TP    FP    TN   FN PPV     NPV      F1
#>  81.25 93.75 31.25 68.75 6.25  75 91.6667 83.3333

# Should return:
#      C    TP    FP    TN   FN PPV     NPV      F1 
#  81.25 93.75 31.25 68.75 6.25  75 91.6667 83.3333 

logitfit(y, p, .4)
#>     C    TP    FP    TN   FN     PPV     NPV      F1
#>  87.5 93.75 18.75 81.25 6.25 83.3333 92.8571 88.2353
# Should return:
#     C    TP    FP    TN   FN     PPV     NPV      F1
#  87.5 93.75 18.75 81.25 6.25 83.3333 92.8571 88.2353

```
