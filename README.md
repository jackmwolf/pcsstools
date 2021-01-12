
<!-- README.md is generated from README.Rmd. Please edit that file -->

grass
=====

<!-- badges: start -->
<!-- badges: end -->

grass (Genetic Regression Approximation through Summary Statistics) is
an in-development R package to describe various regression models using
only genome-wide association study (GWAS) summary statistics.

Currently, grass supports the linear modeling of complex phenotypes
defined via functions of other phenotypes. Supported functions include:

-   linear combinations (E.g.
    *ϕ*<sub>1</sub>*y*<sub>1</sub> + *ϕ*<sub>2</sub>*y*<sub>2</sub>)
-   products (E.g. *y*<sub>1</sub> ∘ *y*<sub>2</sub>)
-   logical combinations (E.g. *y*<sub>1</sub> ∧ *y*<sub>2</sub> or
    *y*<sub>1</sub> ∨ *y*<sub>2</sub>)

Installation
------------

grass is not currently available on CRAN.

You can install the in-development version of grass from
[GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("jackmwolf/grass")

Example
-------

Let’s model the first principal component of
*y*<sub>1</sub>, *y*<sub>2</sub> and *y*<sub>3</sub> using summary
statistics.

    library(grass)
    dat <- grass::cont_data[c("g", "x", "y1", "y2", "y3")]
    head(dat)
    #>   g           x        y1         y2         y3
    #> 1 0  0.08613585  1.248199 -2.4208709  0.8490036
    #> 2 0 -1.19721387  2.277541 -0.8876128 -3.4783177
    #> 3 0 -0.89131854  2.011928 -2.4584807  0.8973664
    #> 4 1  1.19449979 -1.535056 -1.2870269  3.1705261
    #> 5 1  0.22346186  1.431919 -1.9866995  0.5214276
    #> 6 0  0.55236041 -2.282309 -0.2572408 -0.9526495

First, we need our assumed summary statistics.

    means <- colMeans(dat)
    covs  <- cov(dat)
    n     <- nrow(dat)

In addition, we need our weights, the first principal component vector
of the responses’ covariance matrix.

    SigmaY <- covs[c("y1", "y2", "y3"), c("y1", "y2", "y3")]
    phi <- eigen(SigmaY)$vectors[, 1]

Then, we can calculate the linear model by using `calculate_lm_combo`

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

The difference in test statistic estimates can be attributed to the
non-uniqueness of optimal principal component vectors by multiplication
by  − 1.

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
