
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

References
----------

Following are the key references for the functions in this package

-   Wolf, J.M., Barnard, M., Xueting, X., Ryder, N., Westra, J., and
    Tintle, N. (2020). Computationally efficient, exact,
    covariate-adjusted genetic principal component analysis by
    leveraging individual marker summary statistics from large biobanks.
    *Pacific Symposium on Biocomputing*, 25, 719-730.
-   Gasdaska A., Friend D., Chen R., Westra J., Zawistowski M.,
    Lindsey W. and Tintle N. (2019) Leveraging summary statistics to
    make inferences about complex phenotypes in large biobanks. *Pacific
    Symposium on Biocomputing*, 24, 391-402.
