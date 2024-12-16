# ci.prop returns valid matrix

    Code
      res
    Output
                      Estimate         SE         LL        UL
      Adjusted Wald  0.1346154 0.03346842 0.06901848 0.2002123
      Wilson with cc 0.1200000 0.03249615 0.06625153 0.2039772

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
        Estimate         SE        z           p         LL        UL
       0.2119565 0.07602892 2.787841 0.005306059 0.06294259 0.3609705

# ci.pairs.prop.bs returns valid matrix

    Code
      res
    Output
             Estimate         SE         z            p          LL          UL
       1 2 -0.2475248 0.04482323 -5.522243 3.346989e-08 -0.35483065 -0.14021885
       1 3 -0.1039604 0.04833562 -2.150803 3.149174e-02 -0.21967489  0.01175409
       2 3  0.1435644 0.04358401  3.293968 9.878366e-04  0.03922511  0.24790360

# ci.slope.prop.bs returns valid matrix

    Code
      res
    Output
          Estimate          SE        z           p          LL         UL
       0.007542293 0.002016793 3.739746 0.000184206 0.003589452 0.01149513

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
                        Estimate exp(Estimate)        z           p       LL
      At low moderator    0.9328      2.541616 2.269824 0.023218266 1.135802
      At high moderator   1.7644      5.838068 2.906507 0.003654887 1.776421
                               UL
      At low moderator   5.687444
      At high moderator 19.186357

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
           Estimate         SE         LL        UL
      Q:  0.3430670 0.13280379 0.06247099 0.5734020
      Y:  0.1769015 0.07290438 0.03126603 0.3151817
      H:  0.2619244 0.10514465 0.04687994 0.4537659
      Y*: 0.1311480 0.05457236 0.02307188 0.2361941

# ci.phi returns valid matrix

    Code
      res
    Output
        Estimate         SE         LL        UL
       0.1229976 0.05477117 0.01462398 0.2285149

# ci.biphi returns valid matrix

    Code
      res
    Output
        Estimate         SE        LL       UL
       0.4145733 0.07551281 0.2508866 0.546141

# ci.tetra returns valid matrix

    Code
      res
    Output
        Estimate         SE        LL        UL
       0.5135167 0.09301703 0.3102345 0.6748546

# ci.kappa returns valid matrix

    Code
      res
    Output
                    Estimate         SE        LL        UL
      IC kappa:    0.6736597 0.07479965 0.5270551 0.8202643
      Cohen kappa: 0.6756757 0.07344761 0.5317210 0.8196303

# ci.agree returns valid matrix

    Code
      res
    Output
        Estimate         SE        LL        UL
       0.7333333 0.05333333 0.6132949 0.8226025

# ci.popsize returns valid matrix

    Code
      res
    Output
       Estimate       SE   LL   UL
           2908 49.49071 2818 3012

# test.prop returns valid matrix

    Code
      res
    Output
       Estimate        z          p
           0.45 2.515576 0.01188379

# test.prop2 returns valid matrix

    Code
      res
    Output
       Estimate        z           p
           -0.3 2.899726 0.003734895

# test.prop.bs returns valid matrix

    Code
      res
    Output
       Chi-square df            p
         17.41071  2 0.0001656958

# test.prop.ps returns valid matrix

    Code
      res
    Output
       Estimate        z          p
           0.07 2.108346 0.03500109

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
               Estimate         SE        LL        UL
      G1      0.8666667 0.02880329 0.6974555 0.9481141
      G2      0.5000000 0.05590170 0.2523379 0.6851621
      G1 - G2 0.3666667 0.06288585 0.1117076 0.6088621

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
                       Estimate          LL         UL
      G(1,2)         0.56666667  0.46601839  0.6524027
      G(1,3)         0.50000000  0.39564646  0.5911956
      G(2,3)         0.86666667  0.79701213  0.9135142
      G(1,2)-G(1,3)  0.06666667  0.00580397  0.1266464
      G(1,2)-G(2,3) -0.30000000 -0.40683919 -0.1891873
      G(2,3)-G(1,3) -0.36666667 -0.46222023 -0.2662566
      G(3)           0.64444444  0.57382971  0.7068720

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
               318

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
               416

