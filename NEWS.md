statpsych version 1.8.0 (Release date: 2025/06/06)
===========

Changes:

* New functions:
    * test.mono.median.bs -- Hypothesis test of a monotonic trend in medians for an ordered between-subjects factor
    * ci.slope.median.bs -- Computes confidence interval for the slope of medians in a one-factor experimental design with a quantitative between-subjects factor
    * size.ci.cv -- Computes sample size for a coefficient of variation confidence interval in a 1-group design
    * size.ci.median -- Computes sample size for a median confidence interval in a 1-group design
    * size.ci.median2 -- Computes sample size for a median difference confidence interval in a 2-group design
    * size.ci.lc.median.bs -- Computes sample size for a linear contrast of medians confidence interval in a between-subjects design
    * size.ci.gen -- Computes sample size for a confidence interval for any type of parameter in a 1-group design
    * size.ci.gen2 -- Computes sample size for a confidence interval for the difference in any type of parameter in a 2-group design
    * size.test.gen -- Computes sample size for a test of any type of parameter in a 1-group design
    * size.test.gen2 -- Computes sample size for a test of the difference of any type of parameter in a 2-group design
* Modifications:
    * ci.prop now adds an exact confidence interval
* Deleted functions:
    * ci.mean1 has been replaced with ci.mean
    * ci.prop1 has been replaced with ci.prop
    * ci.median1 has been replaced with ci.median
    * ci.stdmean1 has been replaced with ci.stdmean
    * pi.var.upper has been replaced with pi.var
    

statpsych version 1.7.0 (Release date: 2024/12/17)
===========

Changes:

* New functions:
    * size.ci.oddsratio -- Computes sample size for an odds ratio confidence interval
    * size.ci.yule -- Computes sample size for a Yule's Q confidence interval
    * size.ci.phi -- Computes sample size for a phi correlation confidence interval
    * size.ci.biphi -- Computes sample size for a biserial-phi correlation confidence interval
    * size.ci.ancova2 -- Computes sample size for a 2-group ANCOVA confidence interval
    * size.ci.slope.gen -- Computes sample size for a slope coefficient confidence interval in a general model
    * size.test.ancova2 -- Computes sample size for a 2-group ANCOVA hypothesis test
    * size.test.slope.gen -- Computes sample size for a slope coefficient hypothesis test in a general model
    * signal -- Computes parameter estimates in a Yes/No signal detection study
    * exp.slope -- Computes confidence intervals for exp(B) and 100(exp(B) - 1)%
    * ci.bayes.cor -- Computes Bayesian credible interval for a Pearson or partial correlation with a skeptical prior
    * ci.bayes.spcor -- Computes Bayesian credible interval for a semipartial correlation with a skeptical prior
    * pi.var -- Computes one-sided or two-sided prediction limits for an estimated variance in a future study (will replace pi.var.upper)
* Modifications
    * size.ci.prop2 can now solve for equal or unequal sample sizes.  Requires a new argumnet, R, specifying the ratio of sample sizes and now returns a 2-column matrix.
    * size.ci.ratio.prop2 can now solve for equal or unequal sample sizes.  Requires a new argument, R, specifying the ratio of sample sizes and now returns a 2-column matrix.
    * size.test.cor2 can now solve for equal or unequal sample sizes.  Requires a new argument, R, specifying the ratio of sample sizes, and returns a 2-column matrix
    * pi.cor now has options for one-sided and two-sided prediction limits.  It requires a new argument, type
    * pi.prop now has options for one-sided and two-sided prediction limits.  It requires a new argument, type
    * the definition of the subscripts in the ci.2x2.mean.mixed, ci.2x2.median.mixed, and ci.2x2.stdmean.mixed functions have been changed so that the first subscript now specifies factor A and the second subscript specified factor B. 
* Error Corrections:
    * ci.2x2.stdmean.mixed -- corrected an error in the standard error computation
    * size.test.lc.ancova -- corrected a minor error in the sample size formula


statpsych version 1.6.0 (Release date: 2024/07/08)
===========

Changes:

* New functions:
    * ci.mean.fpc -- Computes confidence interval for a mean with a finite population correction
    * ci.prop.fpc -- Computes confidence interval for a proportion with a finite population correction
    * ci.poisson -- Computes confidence interval for a Poisson rate
    * ci.ratio.poisson2 -- Computes confidence interval for a ratio of Poisson rates in a 2-group design
    * ci.bscor -- Computes confidence interval for a biserial correlation
    * pi.cor -- Computes prediction interval for an estimated correlation
    * pi.prop -- Computes prediction interval for an estimated proportion
    * test.cor -- Hypothesis test for a Pearson or partial correlation (for zero or non-zero null hypotheses)
    * test.spear -- Hypothesis test for a Spearman correlation (for zero or non-zero null hypotheses)
    * test.cor2 -- Hypothesis test for a 2-group Pearson or partial correlation difference
    * test.spear2 -- Hypothesis test for a 2-group Spearman correlation difference
    * test.mean -- Hypothesis test for a mean using summary information
    * size.ci.cor2 -- Computes sample size for a 2-group Pearson correlation difference confidence interval
    * size.ci.spear2 -- Computes sample size for a 2-group Spearman correlation difference confidence interval
    * size.ci.tetra -- Computes sample size for a tetrachoric correlation confidence interval
    * size.ci.mean.prior -- Computes sample size for a mean confidence interval using a planning value from a prior study
    * size.ci.prop.prior -- Computes sample size for a proportion confidence interval using a planning value from a prior study
    * size.ci.cor.prior -- Computes sample size for a correlation confidence interval using a planning value from a prior study
    * adj.se -- Computes adjusted standard errors for slope coefficients in an exploratory analysis
    * fitindices -- Computes four SEM fit indices
* Modifications:
    * ci.var.upper now computes an exact upper limit rather than an approximate upper limit
    * power computations are now more accurate for very small effect sizes in the power.cor, power.cor2, power.lc.meanc.bs, power.mean, power.mean2, power.mean.ps, power.prop, power.pro2, and power.prop.ps functions
    * size.test.prop and size.test.prop2 now assume the test statistic will use a continuity correction
    * one-group function names that end with a "1" have been renamed and now exclude the "1" (for naming consistency and to avoid confusion with lower case L).
    * ci.mape2 has been renamed ci.ratio.mape2, and ci.cod2 has been renamed ci.ratio.cod2
    * The ci.phi function now uses a Fisher transformation for improved coverage probability performance


statpsych version 1.5.0 (Release date: 2023/12/11)
===========

Changes:

* New functions:
    * ci.cv1 -- Computes confidence interval for a coefficient of variation
    * ci.ratio.cv2 -- Computes confidence interval for a ratio of coefficients of variation
    * ci.pv -- Computes confidence intervals for positive and negative predictive values with
          retrospective sampling
    * ci.2x2.stdmean.ws -- Computes confidence intervals of standardized effects in a 2x2 
          within-subjects design
    * ci.2x2.stdmean.mixed -- Computes confidence intervals of standardized effects in a 2x2 
          mixed design
    * ci.2x2.median.ws -- Computes confidence intervals of effects in a 2x2 within-subjects 
          design for medians
    * ci.2x2.median.mixed -- Computes confidence intervals of effects in a 2x2 mixed design  
          for medians
    * spearmanbrown -- Computes the reliability of a scale with r2 measurements given the 
          reliability of a scale with r1 measurements
    * size.ci.spear -- Computes the sample size requirement for a Spearman correlation confidence interval
    * size.ci.pbcor -- Computes the sample size requirement for a point-biserial correlation confidence interval
    * size.ci.mape1 -- Computes the sample size requirement for a mean absolute prediction error confidence interval
* Error Corrections:
    * ci.cramer -- corrected an error in the CI computation 
    * ci.lc.stdmean.ws -- corrected an error in the standard error computation
* Modifications:
    * both biased and bias adjusted estimates are now reported in ci.stdmean1, ci.stdmean2, 
           ci.stdmean.ps, ci.stdmean.strat, and ci.2x2.stdmean.bs
    * ci.mape has been renamed ci.mape1


statpsych version 1.4.0 (Release date: 2023/06/26)
===========

Changes:

* New functions:
    * power.prop1 -- Computes power for 1-sample test of proportion for a planned sample size
    * power.prop2 -- Computes power for 2-sample test of proportion for planned sample sizes
    * power.prop.ps -- Computes power for paired-samples test of proportion for a planned 
      sample size
    * power.mean1 -- Computes power for 1-sample t-test for a planned sample size
    * power.mean2 -- Computes power for 2-sample t-test for planned sample sizes
    * power.mean.ps -- Computes power for paired-samples t-test for a planned sample size
    * power.lc.mean.bs -- Computes power of a test for a linear contrast of means for planned 
      sample sizes in a between-subjects design
    * power.cor1 -- Computes power for 1-sample test of correlation for a planned sample size
    * power.cor2 -- Computes power for 2-sample test of correlations for planned sample sizes
    * ci.cqv1 -- Computes confidence interval for a population coefficient of qualitative 
      variation 
    * ci.prop1.inv -- Computes confidence interval for a population proportion using inverse  
      sampling
    * ci.prop2.inv -- Computes confidence interval for a difference in population proportions 
      using inverse sampling
    * ci.agree.3rater -- Computes confidence intervals for a 3-rater design with dichotomous
      ratings
    * ci.ratio.sd2 -- Computes robust confidence interval for ratio of standard deviations in a 
      2-group design
    * size.test.cor2 -- Computes sample size for a test of equal Pearson or partial correlation in a 
      2-group design
    * size.test.cronbach2 -- Computes sample size to test equality of Cronbach reliability    
      coefficients in a 2-group design
    * size.ci.cronbach2 -- Computes sample size for a 2-group Cronbach reliability difference 
      confidence interval
    * size.ci.etasqr -- Computes sample size for an eta-squared confidence interval
    * size.ci.indirect -- Computes sample size for an indirect effect confidence interval
    * ci.mape2 -- Computes confidence interval for a ratio of mean absolute prediction errors in
      a 2-group design
    * ci.rel2 -- Computes confidence interval for a 2-group reliability difference
    * ci.cronbach2 -- Computes confidence interval for a difference in Cronbach reliabilities in  
      a 2-group design
    * ci.2x2.stdmean.bs -- Computes confidence intervals of standardized effects in a 2x2 
      between-subjects design for means
    * ci.2x2.median.bs -- Computes confidence intervals of effects in a 2x2 between-subjects 
      design for medians
    * pi.var.upper -- Computes upper prediction limit for an estimated variance
    * ci.bayes.normal -- Computes Bayesian credible interval for any parameter estimator with 
      a normal sampling distributuion using a Normal prior distribution
    * ci.bayes.prop1 -- Computes Bayesian credible interval for a single proportion using a Beta
      prior distribution
*  Modifications:
    * Corrected Example output in ci.reliability and ci.prop.ps
    * SE added to output in:  ci.cronbach, ci.oddsratio, ci.yule, ci.etasqr, ci.rsqr, ci.spear2, 
      ci.cor2, ci.cor.dep, ci.cod1, ci.mad1, ci.mape, ci.agree2, ci.pbcor, and ci.tetra
    * Improved accuracy in size.ci.rsqr
    * Three generalized Yule coefficients added to ci.yule
    * The ci.prop.ps, ci.ratio.prop.ps, and ci.2x2.prop.mixed functions now define proportions 
      for the y = 1 category rather than the y = 0 category.


statpsych v1.3.0 (Release date: 2023/01/01)
==============

Changes:

* New functions:
    * ci.theil -- Theil-Sen estimate and confidence interval for slope
    * sim.ci.median2 -- Simulates confidence interval coverage probability for a median difference in a two-group design
    * sim.ci.median.ps -- Simulates confidence interval coverage probability for a median difference in a paired design
    * sim.ci.stdmean2 -- Simulates confidence interval coverage probability for a standardized mean difference in a two-group design
    * pi.score.ps -- Prediction interval for difference of scores in a 2-level within-subjects experiment
* Updated outputs:
    * ci.cod1 -- first column is 'Estimate', no longer 'COD'
    * ci.cod2 -- first column is 'Estimate', no longer 'COD1'
    * ci.cramer -- first column is 'Estimate', no longer 'Cramer's V'
    * ci.lc.stdmean.bs -- now returns 3 rows, adding sample size for group 1 standardizer
    * ci.lc.stdmean.ws -- now returns two rows, one for each standardizer
    * ci.mad1 -- first column is 'Estimate', no longer 'MAD'
    * ci.mape -- first column is 'Estimate', no longer 'MAPE'
    * size.ci.lc.stdmean.bs -- now returns two rows, one for each standardizer
    * size.ci.lc.stdmean.ws -- now returns two rows, one for each standardizer
    * size.ci.stdmean2 -- now returns two rows, one for each standardizer
    * size.ci.stdmean.ps -- now returns two rows, one for each standardizer
    * ci.mann -- now returns a confidence interval for P(y1 > y2) rather than P(y1 < y2).
    

statpsych v1.2.0 (Release date: 2022/08/15)
==============

Changes:

* New functions:
    * ci.cramer - Confidence interval for Cramer's V
    * ci.2x2.mean.bs - Confidence intervals for effects in a 2x2 between-subjects design for means
    * ci.2x2.mean.ws - Confidence intervals for effects in a 2x2 within-subjects design for means
    * ci.2x2.mean.mixed - Confidence intervals for effects in a 2x2 mixed design for means
    * ci.2x2.prop.bs - Confidence intervals for effects in a 2x2 between-subjects design for proportions
    * ci.2x2.prop.mixed - Confidence intervals for effects in a 2x2 mixed design for proportions
    * sim.ci.mean1 – Simulation of confidence interval for a mean
    * sim.ci.mean2 – Simulation of confidence interval for mean difference in a two-group design
    * sim.ci.mean.ps – Simulation of confidence interval for mean difference in a paired-samples design
    * sim.ci.median1 – Simulation of confidence interval for a single median
    * sim.ci.cor – Simulation of confidence interval for a Pearson correlation
    * sim.ci.spear – Simulation of confidence interval for a Spearman correlation
* Modifications:
    * The ci.prop.ps function now outputs an adjusted point estimate of the proportion difference, as stated in the documentation, rather than an unadjusted estimate
    * The ci.cor, ci.cor2, and ci.cor.dep functions now uses a bias adjustment to reduce the bias of the Fisher transformed correlations
    * The ci.median1 function now uses the same standard error formula as the ci.median2, ci.ratio.median2, and ci.median.ps functions
* Error Correction:
    * ci.indirect -- Corrected an error in the standard error computation
    

statpsych v1.1.0 (Release date: 2022/06/30)
==============

Changes:

* New functions:
    * ci.agree2 - Confidence interval for G-index difference in a 2-group design
    * ci.cod2 - Confidence interval for a ratio of dispersion coefficients in a 2-group
    * ci.etasqr - Confidence interval for eta-squared
    * ci.lc.gen.bs - Confidence interval for a linear contrast of parameters in a between-subjects design
    * ci.lc.glm - Confidence interval for a linear contrast of general linear model parameters
    * ci.reliability - Confidence interval for a reliability coefficient
    * ci.rsqr - Confidence interval for squared multiple correlation
    * ci.sign1 - Confidence interval for the parameter of the one-sample sign test
    * ci.slope.mean.bs - Confidence interval for the slope of means in a single-factor design with a quantitative between-subjects factor
    * test.kurtosis - Computes Monte Carlo p-value for test of excess kurtosis
    * test.skew - Computes Monte Carlo p-value for test of skewness
    * test.mono.mean.bs - Test of a monotonic trend in means for an ordered between-subjects factor
    * test.mono.prop.bs - Test of monotonic trend in proportions for an ordered between-subjects
    * etasqr.gen.2way - Computes generalized eta-squared estimates in a two-factor design
* Updated documentation for consistency
* Changed arguments for some functions for consistency
    * size.test.cronbach now takes (alpha, pow, rel, r, h) rather than (alpha, pow, rel, a, h)
    * ci.cronbach now takes (alpha, rel, r, n) rather than (alpha, rel, a, n)
* Changed some of the column names in returned matrixes for consistency:
    * ci.median.ps, the last column is now "COV" rather than "cov"

statpsych 1.0.0 (Release date: 2021/09/09)
==============

* Initial release
