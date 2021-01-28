
<!-- README.md is generated from README.Rmd. Please edit that file -->

# grass

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/jackmwolf/grass.svg?branch=master)](https://travis-ci.com/jackmwolf/grass)
<!-- badges: end -->

`grass` (Genetic Regression Approximation through Summary Statistics) is
an in-development R package to describe various regression models using
only pre-computed summary statistics (PCSS) from genome-wide association
studies (GWASs) and PCSS repositories such as
[GeneAtlas](http://geneatlas.roslin.ed.ac.uk/). This eliminates the
logistic, privacy, and access concerns that accompany the use of
individual patient-level data (IPD).

The following figure highlights the information typically needed to
perform regression analysis on a set of *m* phenotypes with *p*
covariates when IPD is available, and the PCSS that are commonly needed
to approximate this same model in `grass`.

![Data needed for analysis using IPD compared to that when using
PCSS](./man/figures/IPDvsPCSS.png)

Currently, `grass` supports the linear modeling of complex phenotypes
defined via functions of other phenotypes. Supported functions include:

-   linear combinations
    (e.g. *ϕ*<sub>1</sub>*y*<sub>1</sub> + *ϕ*<sub>2</sub>*y*<sub>2</sub>)
-   products (e.g. *y*<sub>1</sub> ∘ *y*<sub>2</sub>)
-   logical combinations (e.g. *y*<sub>1</sub> ∧ *y*<sub>2</sub> or
    *y*<sub>1</sub> ∨ *y*<sub>2</sub>)

## Installation

grass is not currently available on CRAN.

You can install the in-development version of grass from
[GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("jackmwolf/grass")

## Examples

We will walk through two examples using grass to model combinations of
phenotypes using PCSS and then compare our results to those found using
IPD.

    library(grass)

### Principal Component Analysis

Let’s model the first principal component score of three phenotypes
using PCSS.

First, we’ll load in some data. We have a SNP’s minor allele counts
(`g`), a continuous covariate (`x`), and three continuous phenotypes
(`y1`, `y2`, and `y3`).

    dat <- grass::cont_data[c("g", "x", "y1", "y2", "y3")]
    head(dat)
    #>   g           x        y1         y2         y3
    #> 1 0  0.08613585  1.248199 -2.4208709  0.8490036
    #> 2 0 -1.19721387  2.277541 -0.8876128 -3.4783177
    #> 3 0 -0.89131854  2.011928 -2.4584807  0.8973664
    #> 4 1  1.19449979 -1.535056 -1.2870269  3.1705261
    #> 5 1  0.22346186  1.431919 -1.9866995  0.5214276
    #> 6 0  0.55236041 -2.282309 -0.2572408 -0.9526495

First, we need our assumed summary statistics: means, the full
covariance matrix, and our sample size.

    means <- colMeans(dat)
    covs  <- cov(dat)
    n     <- nrow(dat)

Then, we can calculate the linear model by using `model_prcomp()`. Our
`formula` will list all phenotypes as one sum, joined together by `+`
opperators and we indicate that we want the first principal component
score by setting `comp = 1`.

    model_pcss <- model_prcomp(
      y1 + y2 + y3 ~ g + x, comp = 1, n = n, means = means, covs = covs
      )
    model_pcss
    #> Model approximated using Pre-Computed Summary Statistics.
    #> 
    #> Call:
    #> model_prcomp(formula = y1 + y2 + y3 ~ g + x, comp = 1, n = n, 
    #>     means = means, covs = covs)
    #> 
    #> Coefficients:
    #>             Estimate Std. Error t value Pr(>|t|)    
    #> (Intercept)  1.04723    0.06925  15.122  < 2e-16 ***
    #> g            0.57050    0.07806   7.308 5.54e-13 ***
    #> x           -1.94340    0.04899 -39.669  < 2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #> Residual standard error: 1.578 on 997 degrees of freedom
    #> Multiple R-squared:  0.6198, Adjusted R-squared:  0.619 
    #> F-statistic: 812.6 on 2 and 996 DF,  p-value: < 2.2e-16

Here’s the same model using individual patient data.

    pc_1 <- prcomp(x = dat[c("y1", "y2", "y3")])$x[, "PC1"]

    mod_ipd <- lm(pc_1 ~ 1 + g + x, data = dat)
    summary(mod_ipd)
    #> 
    #> Call:
    #> lm(formula = pc_1 ~ 1 + g + x, data = dat)
    #> 
    #> Residuals:
    #>    Min     1Q Median     3Q    Max 
    #> -5.117 -1.075 -0.011  1.076  6.439 
    #> 
    #> Coefficients:
    #>             Estimate Std. Error t value Pr(>|t|)    
    #> (Intercept)  0.33001    0.06925   4.765 2.16e-06 ***
    #> g           -0.57050    0.07806  -7.308 5.54e-13 ***
    #> x            1.94340    0.04899  39.669  < 2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #> Residual standard error: 1.578 on 997 degrees of freedom
    #> Multiple R-squared:  0.6198, Adjusted R-squared:  0.619 
    #> F-statistic: 812.6 on 2 and 997 DF,  p-value: < 2.2e-16

In this case, our coefficient estimates for `g` and `x` are off by a
factor of -1; this is because we picked the opposite vector of principal
component weights to `prcomp`. This distinction in sign is arbitrary
(see the note in `?prcomp`).

### Logical Combination

In this example we will approximate a linear model where our response is
the logical combination “*y*<sub>1</sub> or *y*<sub>1</sub>”
(*y*<sub>1</sub> ∨ *y*<sub>2</sub>).

First we need data with binary phenotypes.

    dat <- grass::bin_data[c("g", "x", "y1", "y2")]
    head(dat)
    #>   g          x y1 y2
    #> 1 0 -0.9161478  1  0
    #> 2 0  1.2496985  0  1
    #> 3 1 -1.2708514  0  0
    #> 4 2  0.0832760  0  1
    #> 5 0  0.4686342  0  1
    #> 6 2  0.4620154  0  1

Once again we will organized our assumed PCSS.

    means <- colMeans(dat)
    covs  <- cov(dat)
    n     <- nrow(dat)

We also need to describe the distributions of both of our predictors
through objects of class `predictor`. (See `?new_predictor`.) grass has
shortcut functions to create `predictor` objects for common types of
variables, which we will use to create a list of `predictor`s.

    predictors <- list(
      g = new_predictor_snp(maf = means["g"] / 2),
      x = new_predictor_normal(mean = means["x"], sd = sqrt(covs["x", "x"]))
    )
    class(predictors[[1]])
    #> [1] "predictor"

Then we can approximate the linear model using `model_or()`.

    model_or(y1 | y2 ~ g + x, n = n, means = means, covs = covs, predictors = predictors)
    #> Model approximated using Pre-Computed Summary Statistics.
    #> 
    #> Call:
    #> model_or(formula = y1 | y2 ~ g + x, n = n, means = means, covs = covs, 
    #>     predictors = predictors)
    #> 
    #> Coefficients:
    #>             Estimate Std. Error t value Pr(>|t|)    
    #> (Intercept)  0.66927    0.01878  35.629  < 2e-16 ***
    #> g           -0.09104    0.02200  -4.137 3.81e-05 ***
    #> x            0.19003    0.01412  13.455  < 2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #> Residual standard error: 0.4441 on 997 degrees of freedom
    #> Multiple R-squared:  0.172,  Adjusted R-squared:  0.1703 
    #> F-statistic: 103.5 on 2 and 996 DF,  p-value: < 2.2e-16

And here’s the result we would get using IPD:

    model_ipd <- lm(y1 | y2 ~ g + x, data = dat)
    summary(model_ipd)
    #> 
    #> Call:
    #> lm(formula = y1 | y2 ~ g + x, data = dat)
    #> 
    #> Residuals:
    #>     Min      1Q  Median      3Q     Max 
    #> -1.0857 -0.4393  0.1538  0.3792  0.8597 
    #> 
    #> Coefficients:
    #>             Estimate Std. Error t value Pr(>|t|)    
    #> (Intercept)  0.67337    0.01888   35.67  < 2e-16 ***
    #> g           -0.09862    0.02211   -4.46 9.14e-06 ***
    #> x            0.18290    0.01419   12.88  < 2e-16 ***
    #> ---
    #> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    #> 
    #> Residual standard error: 0.4463 on 997 degrees of freedom
    #> Multiple R-squared:  0.1636, Adjusted R-squared:  0.1619 
    #> F-statistic:  97.5 on 2 and 997 DF,  p-value: < 2.2e-16

## Future Work

-   Add functions to check that all required PCSS are provided in
    `means` and `covs` (and `predictors` and `responses`, if applicable)

-   Return output from `calculate_lm` as an `S3` object of class
    `lm_pcss`, not `summary.lm_pcss`. Create a `print` method for this
    class and then define `summary` method for this class (and use
    existing `print.summary.lm_pcss`).

-   Support function notation for linear combinations of phenotypes
    (e.g. `y1 - y2 + 0.5 * y3 ~ 1 + g + x`) instead of requiring a
    seperate vector of weights

-   Support functions using `.` and `-` in the dependent variable
    (e.g. `y1 ~ .`, `y1 ~ . -x`)

-   Write a vignette

## References

Following are the key references for the functions in this package

-   Wolf, J.M., Barnard, M., Xueting, X., Ryder, N., Westra, J., and
    Tintle, N. (2020). Computationally efficient, exact,
    covariate-adjusted genetic principal component analysis by
    leveraging individual marker summary statistics from large biobanks.
    *Pacific Symposium on Biocomputing*, 25, 719-730.
    <https://doi.org/10.1142/9789811215636_0063>.

-   Gasdaska A., Friend D., Chen R., Westra J., Zawistowski M.,
    Lindsey W. and Tintle N. (2019) Leveraging summary statistics to
    make inferences about complex phenotypes in large biobanks. *Pacific
    Symposium on Biocomputing*, 24, 391-402.
    <https://doi.org/10.1142/9789813279827_0036>.
