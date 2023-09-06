
## Test environments
- local R installation, R 4.2.1
- Windows Server 2022, R-devel, 64 bit
- Fedora Linux, R-devel, clang, gfortran
- Ubuntu Linux 20.04.1 LTS, R-release, GCC
- win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

- Found the following (possibly) invalid URLs:
  URL: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6907735/
  From: man/calculate_lm_combo.Rd
        man/calculate_lm.Rd
        man/model_combo.Rd
        man/model_prcomp.Rd
        man/model_singular.Rd
        man/pcsslm.Rd
  Status: 403
  Message: Forbidden
  
- I have checked this URL. It seems to be correct.

## revdepcheck results

There are currently no downstream dependencies for this package.