
<!-- README.md is generated from README.Rmd. Please edit that file -->

grass
=====

<!-- badges: start -->
<!-- badges: end -->

grass (Genetic Regression Approximation through Summary Statistics) is
an in-development R package to describe various regression models using
only summary statistics.

Currently, grass supports the linear modeling of complex phenotypes
defined via functions of other phenotypes. Supported functions include:

-   linear combinations (E.g.
    *ϕ*<sub>1</sub>*y*<sub>1</sub> + *ϕ*<sub>2</sub>*y*<sub>2</sub>)
-   products (E.g. *y*<sub>1</sub> ∘ *y*<sub>2</sub>)
-   logical combinations (E.g. *y*<sub>1</sub> ∧ *y*<sub>2</sub> or
    *y*<sub>1</sub> ∨ *y*<sub>2</sub>)

Installation
------------

You can install the in-development version of grass from
[GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("jackmwolf/grass")

Example
-------
