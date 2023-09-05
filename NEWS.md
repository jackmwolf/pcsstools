# pcsstools (development version)

# pcsstools 0.1.1
This version includes minor changes implemented for CRAN approval. Namely, return values were added to several functions which did not have a documented value in the initial submission.

# pcsstools 0.1.0

This is a new package that defines functions to describe regression models using only pre-computed summary statistics (PCSS) (i.e. means, variances, and covariances) in place of individual participant data.
Possible models include linear models for linear combinations, products, and logical combinations of phenotypes.
These methods were originally presented in 

* Wolf, J.M., Westra, J., and Tintle, N. (2021). Using summary statistics to 
  evaluate  the genetic architecture of multiplicative combinations of initially
  analyzed phenotypes with a flexible choice of covariates.
  *bioRxiv*.
  [https://doi.org/10.1101/2021.03.08.433979](https://doi.org/10.1101/2021.03.08.433979).


* Wolf, J.M., Barnard, M., Xueting, X., Ryder, N., Westra, J., and Tintle, N. 
  (2020). Computationally efficient, exact, covariate-adjusted genetic principal
  component analysis by leveraging individual marker summary statistics from 
  large biobanks. *Pacific Symposium on Biocomputing*, 25, 719-730. 
  [https://doi.org/10.1142/9789811215636_0063](https://doi.org/10.1142/9789811215636_0063).
  
* Gasdaska A., Friend D., Chen R., Westra J., Zawistowski M., Lindsey W. and 
  Tintle N. (2019) Leveraging summary statistics to make inferences about 
  complex phenotypes in large biobanks. *Pacific Symposium on Biocomputing*, 24, 
  391-402.
  [https://doi.org/10.1142/9789813279827_0036](https://doi.org/10.1142/9789813279827_0036).

## Key functions

* `pcsslm()` approximates a linear model of a combination of variables using PCSS

* `anova.pcsslm()` approximates an analysis of variance table for one or more linear models fit using PCSS

