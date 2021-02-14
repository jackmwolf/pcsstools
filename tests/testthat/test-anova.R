context("anova.pcsslm")

test_that("anova.pcsslmlist matches anova.lmlist", {
  ex_data <- pcsstools_example[c("g1", "x1", "x2", "x3", "y1", "y2")]
  
  pcss <- list(
    means = colMeans(ex_data),
    covs = cov(ex_data),
    n = nrow(ex_data)
  )
  
  pcss_mod_full <- pcsslm(y1 + y2 ~ g1 + x1 + x2 + x3, pcss = pcss, phi = c(1, 1))
  pcss_mod_reduced <- update(pcss_mod_full, . ~ . - x1 - x2 - x3)
  
  pcss_anova <- anova(pcss_mod_reduced, pcss_mod_full)
  
  ipd_mod_full <- lm(y1 + y2 ~ g1 + x1 + x2 + x3, data = ex_data)
  ipd_mod_reduced <- update(ipd_mod_full, . ~ . - x1 - x2 - x3)
  
  ipd_anova <- anova(ipd_mod_reduced, ipd_mod_full)
  
  expect_equal(pcss_anova$Res.Df, ipd_anova$Res.Df)
  expect_equal(pcss_anova$RSS, ipd_anova$RSS)
  expect_equal(pcss_anova$Df, ipd_anova$Df)
  expect_equal(pcss_anova$`Sum of Sq`, ipd_anova$`Sum of Sq`)
  expect_equal(pcss_anova$F, ipd_anova$F)
  expect_equal(pcss_anova$`Pr(>F)`, ipd_anova$`Pr(>F)`)
  
})
