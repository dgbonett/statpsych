statpsych v1.3.0 (Release date: 2023/01/06)
==============

Changes:

* New functions:
    * ci.theil -- Theil-Sen estimate and confidence interval for slope
    * sim.ci.median2 -- Simulates confidence interval coverage probability for a median difference in a two-group design
    * sim.ci.median.ps -- Simulates confidence interval coverage probability for a median difference in a paired design
    * sim.ci.stdmean2 -- Simulates confidence interval coverage probability for a standardized mean difference in a 2-group design
    * pi.score.ps -- Prediction interval for difference of scores in a 2-level within-subjects experiment
* Updated outputs:
    * ci.cod1 -- first column is 'Estimate', no longer 'COD'
    * ci.cramer -- first column is 'Estimate', no longer 'Cramer's V'
    * ci.mad1 -- first column is 'Estimate', no longer 'MAD'
    * ci.mape -- first column is 'Estimate', no longer 'MAPE'
    * ci.lc.stdmean.bs -- now returns 3 rows, adding results for group 1 standardizer
    * ci.lc.stdmean.ws -- now returns two rows, adding results for level 1 standardizer
    * size.ci.lc.stdmean.bs -- now returns two rows, adding result for  group 1 standardizer
    * size.ci.lc.stdmean.ws -- now returns two rows, adding result for level 1 standardizer
    * size.ci.stdmean2 -- now returns two rows, adding result for group 1 standardizer
    * size.ci.stdmean.ps -- now returns two rows, adding result for level 1 standardizer
    * ci.mann -- now returns a confidence interval for P(y1 > y2) rather than P(y1 < y2).
* Error correction:
    * ci.lc.std.mean.ws -- corrected an error in the standard error
    

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
    * sim.ci.mean1 – Simulation of confidence interval for a single mean
    * sim.ci.mean2 – Simulation of confidence interval for mean difference in a two-group design
    * sim.ci.mean.ps – Simulation of confidence interval for mean difference in a paired-samples design
    * sim.ci.median1 – Simulation of confidence interval for a single median
    * sim.ci.cor – Simulation of confidence interval for a Pearson correlation
    * sim.ci.spear – Simulation of confidence interval for a Spearman correlation
* Modifications:
    * The ci.prop.ps function now outputs an adjusted point estimate of the proportion difference, as stated in the documentation, rather than an unadjusted estimate
    * The ci.cor, ci.cor2, and ci.cor.dep functions now uses a bias adjustment to reduce the bias of the Fisher transformed correlations
    * The ci.median1 function now uses the same standard error formula as the ci.median2, ci.ratio.median2, and ci.median.ps functions
* Error correction:
    * Corrected an error for the standard error computation in the ci.indirect function
    

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
