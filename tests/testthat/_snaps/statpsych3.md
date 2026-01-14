# ci.prop returns valid matrix

    Code
      res
    Output
                      Estimate         SE         LL        UL
      Adjusted Wald  0.1346154 0.03346842 0.06901848 0.2002123
      Wilson with cc 0.1200000 0.03249615 0.06625153 0.2039772
      Exact          0.1200000 0.03249615 0.06356890 0.2002357

# ci.pairs.mult  returns valid matrix

    Code
      res
    Output
              Estimate         SE          LL         UL
       1 2  0.14285714 0.04731825  0.05011508 0.23559920
       1 3  0.10963455 0.04875715  0.01407230 0.20519680
       2 3 -0.03322259 0.04403313 -0.11952594 0.05308076

# ci.prop2 returns valid matrix

    Code
      res
    Output
         Estimate         SE          LL        UL
       0.09210526 0.04476077 0.004375769 0.1798348

# ci.ratio.prop2 returns valid matrix

    Code
      res
    Output
       Estimate       LL       UL
       1.658824 1.017253 2.705025

# ci.lc.prop.bs returns valid matrix

    Code
      res
    Output
        Estimate         SE      z       p         LL        UL
       0.2119565 0.07602892 2.7878 0.00531 0.06294259 0.3609705

# ci.pairs.prop.bs returns valid matrix

    Code
      res
    Output
             Estimate         SE       z       p          LL          UL
       1 2 -0.2475248 0.04482323 -5.5222 0.00000 -0.35483065 -0.14021885
       1 3 -0.1039604 0.04833562 -2.1508 0.03149 -0.21967489  0.01175409
       2 3  0.1435644 0.04358401  3.2940 0.00099  0.03922511  0.24790360

# ci.slope.prop.bs returns valid matrix

    Code
      res
    Output
          Estimate          SE      z       p          LL         UL
       0.007542293 0.002016793 3.7397 0.00018 0.003589452 0.01149513

# ci.prop.ps returns valid matrix

    Code
      res
    Output
       Estimate         SE         LL         UL
          -0.44 0.09448809 -0.6251933 -0.2548067

# ci.ratio.prop.ps returns valid matrix

    Code
      res
    Output
       Estimate        LL       UL
         0.3125 0.1725141 0.566077

# ci.condslope.log returns valid matrix

    Code
      res
    Output
                        Estimate exp(Estimate)      z       p       LL        UL
      At low moderator    0.9328      2.541616 2.2698 0.02322 1.135802  5.687444
      At high moderator   1.7644      5.838068 2.9065 0.00365 1.776421 19.186357

# ci.oddsratio returns valid matrix

    Code
      res
    Output
       Estimate        SE       LL       UL
       2.044451 0.6154578 1.133267 3.688254

# ci.yule returns valid matrix

    Code
      res
    Output
          Estimate     SE    LL    UL
      Q:     0.343 0.1328 0.062 0.573
      Y:     0.177 0.0729 0.031 0.315
      H:     0.262 0.1051 0.047 0.454
      Y*:    0.131 0.0546 0.023 0.236

# ci.phi returns valid matrix

    Code
      res
    Output
       Estimate     SE    LL    UL
          0.123 0.0548 0.015 0.229

# ci.biphi returns valid matrix

    Code
      res
    Output
       Estimate     SE    LL    UL
          0.415 0.0755 0.251 0.546

# ci.tetra returns valid matrix

    Code
      res
    Output
       Estimate    SE   LL    UL
          0.514 0.093 0.31 0.675

# ci.kappa returns valid matrix

    Code
      res
    Output
                   Estimate     SE    LL    UL
      IC kappa:       0.674 0.0748 0.527 0.821
      Cohen kappa:    0.676 0.0734 0.532 0.820

# ci.agree returns valid matrix

    Code
      res
    Output
       Estimate      SE     LL     UL
         0.7333 0.05333 0.6133 0.8226

# ci.popsize returns valid matrix

    Code
      res
    Output
       Estimate    SE   LL   UL
           2908 49.49 2818 3012

# test.prop returns valid matrix

    Code
      res
    Output
       Estimate      z       p
        0.01188 2.5156 0.01188

# test.prop2 returns valid matrix

    Code
      res
    Output
       Estimate      z       p
           -0.3 2.8997 0.00373

# test.prop.bs returns valid matrix

    Code
      res
    Output
       Chi-square df       p
          17.4107  2 0.00017

# test.prop.ps returns valid matrix

    Code
      res
    Output
       Estimate      z     p
           0.07 2.1083 0.035

# size.ci.prop returns valid numeric

    Code
      res
    Output
       Sample size
                93

# size.ci.prop2 returns valid numeric

    Code
      res
    Output
        n1  n2
       383 192

# size.ci.ratio.prop2 returns valid numeric

    Code
      res
    Output
        n1  n2
       704 352

# size.ci.lc.prop.bs returns valid numeric

    Code
      res
    Output
       Sample size per group
                          87

# size.ci.prop.ps returns valid numeric

    Code
      res
    Output
       Sample size
               118

# size.ci.ratio.prop.ps returns valid numeric

    Code
      res
    Output
       Sample size
                67

# size.ci.agree returns valid numeric

    Code
      res
    Output
       Sample size
               139

# size.test.prop returns valid numeric

    Code
      res
    Output
       Sample size
                65

# size.test.prop2 returns valid numeric

    Code
      res
    Output
       Sample size per group
                         109

# size.test.lc.prop.bs returns valid numeric

    Code
      res
    Output
       Sample size per group
                         105

# size.equiv.prop2 returns valid numeric

    Code
      res
    Output
       Sample size per group
                         288

# size.supinf.prop2 returns valid numeric

    Code
      res
    Output
       Sample size per group
                         408

# size.test.prop.ps returns valid numeric

    Code
      res
    Output
       Sample size
               177

# size.equiv.prop.ps returns valid numeric

    Code
      res
    Output
       Sample size
               173

# size.supinf.prop.ps returns valid numeric

    Code
      res
    Output
       Sample size
               227

# iqv returns valid matrix

    Code
      res
    Output
         Simpson    Berger   Shannon
       0.7367908 0.5045045 0.7353931

# test.mono.prop.bs returns valid matrix

    Code
      res
    Output
            Estimate         SE         LL        UL
       1 2 0.1764706 0.06803446 0.01359747 0.3393437
       2 3 0.1862745 0.06726135 0.02525219 0.3472968
       3 4 0.1960784 0.05493010 0.06457688 0.3275800

# ci.agree2 returns valid matrix

    Code
      res
    Output
              Estimate      SE     LL     UL
      G1        0.8667 0.02880 0.6975 0.9481
      G2        0.5000 0.05590 0.2523 0.6852
      G1 - G2   0.3667 0.06289 0.1117 0.6089

# power.prop returns valid matrix

    Code
      res
    Output
           Power
       0.7156166

# power.prop2 returns valid matrix

    Code
      res
    Output
           Power
       0.4998959

# power.prop.ps returns valid matrix

    Code
      res
    Output
           Power
       0.6877704

# ci.prop.inv returns valid matrix

    Code
      res
    Output
         Estimate         SE         LL        UL
       0.07462687 0.03145284 0.02467471 0.1479676

# ci.prop2.inv returns valid matrix

    Code
      res
    Output
       Estimate         SE         LL        UL
       0.161385 0.05997618 0.05288277 0.2879851

# ci.agree.3rater returns valid matrix

    Code
      res
    Output
                    Estimate      LL      UL
      G(1,2)          0.5667  0.4660  0.6524
      G(1,3)          0.5000  0.3956  0.5912
      G(2,3)          0.8667  0.7970  0.9135
      G(1,2)-G(1,3)   0.0667  0.0058  0.1266
      G(1,2)-G(2,3)  -0.3000 -0.4068 -0.1892
      G(2,3)-G(1,3)  -0.3667 -0.4622 -0.2663
      G(3)            0.6444  0.5738  0.7069

# ci.bayes.prop returns valid matrix

    Code
      res
    Output
       Posterior mean Posterior SD        LL        UL
            0.1723577   0.03419454 0.1111747 0.2436185

# ci.pv returns valid matrix

    Code
      res
    Output
             Estimate        LL        UL
      PPV:  0.7640449 0.5838940 0.8819671
      NPV:  0.9779978 0.9623406 0.9872318

# ci.prop.fpc returns valid matrix

    Code
      res
    Output
        Estimate        SE         LL        UL
       0.1346154 0.0290208 0.07773565 0.1914951

# ci.poisson returns valid matrix

    Code
      res
    Output
       Estimate        SE       LL      UL
       4.380952 0.9684952 2.777148 6.57358

# ci.ratio.poisson2 returns valid matrix

    Code
      res
    Output
       Estimate       LL       UL
           5.13 1.939576 13.71481

# pi.prop returns valid matrix

    Code
      res
    Output
              LL       UL
       0.1390955 0.337095

# size.ci.tetra returns valid matrix

    Code
      res
    Output
       Sample size
               296

# size.ci.prop.prior returns valid matrix

    Code
      res
    Output
       Sample size
               384

# size.ci.oddsratio example

    Code
      res
    Output
       Sample size
               356

# size.ci.yule example

    Code
      res
    Output
       Sample size
               354

# size.ci.phi example

    Code
      res
    Output
       Sample size
               418

# ci.diversity example

    Code
      res
    Output
              Estimate      SE     LL     UL
      Berger    0.5598 0.01587 0.5287 0.5909
      Simpson   0.7722 0.01229 0.7481 0.7963
      Shannon   0.7292 0.01224 0.7052 0.7532

# ci.lc.prop.scheffe example

    Code
      res
    Output
        Estimate         SE      z       p         LL        UL
       0.2119565 0.07602892 2.7878 0.02053 0.02585698 0.3980561

