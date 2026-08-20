# Cross-validation of NPLS by classical K-fold CV.

This function runs ACMTF-R with cross-validation. A deterministic K–fold
partition is used: the subjects are split in order into `cvFolds`
groups. For each fold the training set consists of the other folds and
the test set is the current fold.

## Usage

``` r
ncrossreg(X, y, maxNumComponents = 5, maxIter = 120, cvFolds = dim(X)[1])
```

## Arguments

- X:

  Centered tensor of independent data

- y:

  Centered dependent variable

- maxNumComponents:

  Maximum number of components to investigate (default 5).

- maxIter:

  Maximum number of iterations (default 100).

- cvFolds:

  Number of folds to use in the cross-validation. For example, if
  `cvFolds` is 5, then the subjects are deterministically partitioned
  into 5 groups (each CV iteration uses 4/5 for training and 1/5 for
  testing). Default: equal to the number of subjects (i.e.
  jack-knifing).

## Value

A list with two elements: - **varExp**: a tibble with the
variance–explained (for X and Y) per number of components. - **RMSE**: a
tibble with the RMSE (computed over the unified CV prediction vector)
per number of components.

## Examples

``` r
set.seed(123)
X <- array(rnorm(25 * 5 * 4), dim = c(25, 5, 4))
y <- rnorm(25)  # Random response variable
result = ncrossreg(X, y, cvFolds=2, maxNumComponents=2)
```
