
<!-- README.md is generated from README.Rmd. Please edit that file -->

grass
=====

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/jackmwolf/grass_alpha.svg?branch=master)](https://travis-ci.com/jackmwolf/grass_alpha)
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

Installation
------------

grass is not currently available on CRAN.

You can install the in-development version of grass from
[GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("jackmwolf/grass")

Examples
--------

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

First, we need our assumed summary statistics. We are careful so that
the order of `means` and `covs` are the same, and that their last
elements are the phenotypes of interest.

    means <- colMeans(dat)
    covs  <- cov(dat)
    n     <- nrow(dat)

In addition, need our weights. These are the the first principal
component vector of the phenotype covariance matrix, and they are in the
same order as the final elements of `means` and `covs`.

    SigmaY <- covs[c("y1", "y2", "y3"), c("y1", "y2", "y3")]
    phi <- eigen(SigmaY)$vectors[, 1]

Then, we can calculate the linear model by using `calculate_lm_combo()`.

    model_pcss <- calculate_lm_combo(
      means = means, covs = covs, n = n, phi = phi, add_intercept = TRUE
      )
    model_pcss
    #> $beta
    #> (Intercept)           g           x 
    #>    1.047230    0.570501   -1.943404 
    #> 
    #> $sd_beta
    #> (Intercept)           g           x 
    #>  0.06925072  0.07806084  0.04899036 
    #> 
    #> $t_stat
    #> (Intercept)           g           x 
    #>   15.122303    7.308415  -39.669117 
    #> 
    #> $p_val
    #>   (Intercept)             g             x 
    #>  1.142929e-46  5.535041e-13 2.836904e-207 
    #> 
    #> $sigma2
    #> [1] 2.490841

Here’s the same model using individual patient data.

    pc_1 <- prcomp(x = dat[c("y1", "y2", "y3")])$x[, "PC1"]

    mod_ipd <- lm(pc_1 ~ 1 + g + x, data = dat)
    summary(mod_ipd)$coef
    #>               Estimate Std. Error   t value      Pr(>|t|)
    #> (Intercept)  0.3300109 0.06925072  4.765452  2.164147e-06
    #> g           -0.5705010 0.07806084 -7.308415  5.535041e-13
    #> x            1.9434042 0.04899036 39.669117 2.836904e-207
    summary(mod_ipd)$sigma^2
    #> [1] 2.490841

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

Once again we will organize our PCSS such that `means` and `covs` have
the same order with the phenotypes at the end.

    means <- colMeans(dat)
    covs  <- cov(dat)
    n     <- nrow(dat)

We also need to describe the distributions of both of our predictors
through objects of class `predictor`. (See `?new_predictor`.) grass has
shortcut functions to create `predictor` objects for common types of
variables, which we will use to create a list of `predictor`s. (Note
that the order `predictors` matches the order of the predictors in
`means` and `covs`.)

    predictors <- list(
      g = new_predictor_snp(maf = means["g"] / 2),
      x = new_predictor_normal(mean = means["x"], sd = sqrt(covs["x", "x"]))
    )
    class(predictors[[1]])
    #> [1] "predictor"

Then we can approximate the linear model using `approx_or()`.

    model_pcss <- approx_or(
      means = means, covs = covs, n = n, predictors = predictors, add_intercept = TRUE
    )
    model_pcss
    #> $beta
    #> (Intercept)           g           x 
    #>  0.66927326 -0.09103735  0.19003213 
    #> 
    #> $sd_beta
    #> (Intercept)           g           x 
    #>  0.01878447  0.02200364  0.01412366 
    #> 
    #> $t_stat
    #> (Intercept)           g           x 
    #>   35.629066   -4.137377   13.454881 
    #> 
    #> $p_val
    #>   (Intercept)             g             x 
    #> 5.487501e-180  3.809339e-05  4.832955e-38 
    #> 
    #> $sigma2
    #> [1] 0.1972141

And here’s the result we would get using IPD:

    model_ipd <- lm(y1 | y2 ~ 1 + g + x, data = dat)
    summary(model_ipd)$coef
    #>                Estimate Std. Error   t value      Pr(>|t|)
    #> (Intercept)  0.67337349 0.01887924 35.667413 3.007202e-180
    #> g           -0.09862393 0.02211464 -4.459667  9.141119e-06
    #> x            0.18289582 0.01419491 12.884609  3.007696e-35
    summary(model_ipd)$sigma^2
    #> [1] 0.1992089

Future Work
-----------

-   Incorporate support for function notation (E.g.
    `y1 * y2 ~ 1 + g + x`) instead of depending on the order of input
    PCSS

-   Print model output in a more similar format to `summary.lm`

References
----------

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
