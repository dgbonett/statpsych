# ci.mean returns valid matrix

    Code
      res
    Output
       Estimate        SE       LL       UL
           24.5 0.5771157 23.33267 25.66733

# ci.stdmean returns valid matrix

    Code
      res
    Output
       Estimate adj Estimate      SE     LL     UL
         1.2329        1.209 0.21244 0.8165 1.6493

# ci.mean2 returns valid matrix

    Code
      res
    Output
                                   Estimate        SE      t    df p       LL
      Equal Variances Assumed:          5.1 0.7151214 7.1317 48.00 0 3.662152
      Equal Variances Not Assumed:      5.1 0.6846568 7.4490 46.17 0 3.721994
                                         UL
      Equal Variances Assumed:     6.537848
      Equal Variances Not Assumed: 6.478006

# ci.lc.mean.bs returns valid matrix

    Code
      res
    Output
                                   Estimate       SE      t    df       p        LL
      Equal Variances Assumed:        -5.35 1.300136 -4.115 36.00 0.00022 -7.986797
      Equal Variances Not Assumed:    -5.35 1.300136 -4.115 33.52 0.00024 -7.993588
                                          UL
      Equal Variances Assumed:     -2.713203
      Equal Variances Not Assumed: -2.706412

# ci.tukey returns valid matrix

    Code
      res
    Output
           Estimate       SE        t    df     p         LL         UL
       1 2    -4.71 1.142459  -4.1227 36.20 6e-04  -7.501852  -1.918148
       1 3   -13.43 1.104078 -12.1640 36.96 0e+00 -16.125710 -10.734290
       2 3    -8.72 1.228730  -7.0968 37.88 0e+00 -11.717042  -5.722958

# ci.ratio.mean2  returns valid matrix

    Code
      res
    Output
       Mean1    Mean2 Mean1/Mean2        LL      UL
        41.5 36.38462    1.140592 0.9837277 1.32247

# ci.stdmean2  returns valid matrix

    Code
      res
    Output
                               Estimate adj Estimate      SE     LL     UL
      Unweighted standardizer:   1.1745       1.1592 0.28440 0.6171 1.7319
      Weighted standardizer:     1.1745       1.1592 0.28028 0.6251 1.7238
      Group 1 standardizer:      1.1475       1.1176 0.29756 0.5643 1.7307
      Group 2 standardizer:      1.2034       1.1720 0.31205 0.5918 1.8151

# ci.stdmean.strat returns valid matrix

    Code
      res
    Output
                             Estimate adj Estimate      SE      LL     UL
      Weighted standardizer:  -0.0554      -0.0553 0.10023 -0.2518 0.1411
      Group 1 standardizer:   -0.0571      -0.0569 0.10369 -0.2604 0.1461
      Group 2 standardizer:   -0.0536      -0.0569 0.09721 -0.2441 0.1369

# ci.lc.stdmean.bs returns valid matrix

    Code
      res
    Output
                               Estimate adj Estimate      SE      LL      UL
      Unweighted standardizer:  -1.3013      -1.2740 0.36928 -2.0250 -0.5775
      Weighted standardizer:    -1.3013      -1.2740 0.35145 -1.9901 -0.6124
      Group 1 standardizer:     -1.3932      -1.2738 0.48498 -2.3438 -0.4427

# ci.mean.ps returns valid matrix

    Code
      res
    Output
       Estimate       SE      t df     p       LL       UL
            6.8 1.455922 4.6706 29 6e-05 3.822304 9.777696

# ci.ratio.mean.ps returns valid matrix

    Code
      res
    Output
        Mean1 Mean2 Mean1/Mean2      LL       UL
       3.4875 3.075    1.134146 1.09417 1.175583

# ci.stdmean.ps returns valid matrix

    Code
      res
    Output
                                  Estimate adj Estimate      SE     LL     UL
      Unweighted standardizer:      0.5550       0.5433 0.16099 0.2395 0.8706
      Measurement 1 standardizer:   0.5425       0.5254 0.16155 0.2259 0.8591
      Measurement 2 standardizer:   0.5685       0.5505 0.16930 0.2367 0.9003

# ci.lc.stdmean.ws returns valid matrix

    Code
      res
    Output
                               Estimate adj Estimate      SE      LL      UL
      Unweighted standardizer:  -1.3013      -1.2666 0.31479 -1.9182 -0.6843
      Level 1 standardizer:     -1.3932      -1.3375 0.36618 -2.1109 -0.6755

# ci.mad returns valid matrix

    Code
      res
    Output
       Estimate       SE       LL       UL
           12.5 2.876103 7.962667 19.62282

# ci.ratio.mad2 returns valid matrix

    Code
      res
    Output
           MAD1     MAD2 MAD1/MAD2        LL       UL
       5.111111 5.888889 0.8679245 0.4520879 1.666253

# ci.ratio.mad.ps returns valid matrix

    Code
      res
    Output
           MAD1 MAD2 MAD1/MAD2       LL       UL
       12.71429  7.5  1.695238 1.109176 2.590961

# ci.cv returns valid matrix

    Code
      res
    Output
        Estimate         SE        LL        UL
       0.1489796 0.01817373 0.1214381 0.1926778

# ci.ratio.cv2 returns valid matrix

    Code
      res
    Output
       Estimate       LL       UL
       1.389188 1.041478 1.854101

# ci.cod returns valid matrix

    Code
      res
    Output
        Estimate        SE        LL       UL
       0.5921053 0.1814708 0.3813259 1.092679

# ci.median returns valid matrix

    Code
      res
    Output
       Estimate       SE LL UL
             20 4.270922 10 30

# ci.median2 returns valid matrix

    Code
      res
    Output
       Median1 Median2 Median1-Median2       SE        LL          UL
          34.5      43            -8.5 4.316291 -16.95977 -0.04022524

# ci.ratio.median2 returns valid matrix

    Code
      res
    Output
       Median1 Median2 Median1/Median2       LL       UL
            43      37        1.162162 0.927667 1.455933

# ci.lc.median.bs returns valid matrix

    Code
      res
    Output
       Estimate       SE       LL       UL
          35.77 11.67507 12.88727 58.65273

# ci.median.ps returns valid matrix

    Code
      res
    Output
       Median1 Median2 Median1-Median2       SE        LL        UL      SE1      SE2
            13      30             -17 3.362289 -23.58996 -10.41004 3.085608 4.509735
            COV
       9.276849

# ci.ratio.median.ps returns valid matrix

    Code
      res
    Output
       Median1 Median2 Median1/Median2        LL        UL
            13      30       0.4333333 0.3094838 0.6067451

# ci.mann returns valid matrix

    Code
      res
    Output
       Estimate        SE LL        UL
          0.205 0.1401834  0 0.4797544

# ci.random.anova returns valid matrix

    Code
      res
    Output
                      Estimate         LL         UL
      Grand mean     59.200000 49.9363896 68.4636104
      Within SD:      9.166782  8.0509046 10.4373219
      Between SD:     8.585948  8.3239359  8.8562078
      Omega-squared:  0.467317  0.2284142  0.8480383

# ci.cronbach returns valid matrix

    Code
      res
    Output
       Estimate      SE     LL     UL
           0.85 0.02457 0.7971 0.8932

# size.ci.mean returns valid number

    Code
      res
    Output
       Sample size
                43

# size.ci.mean2 returns valid matrix

    Code
      res
    Output
       n1 n2
       47 47

# size.ci.stdmean2 returns valid matrix

    Code
      res
    Output
                                  n1  n2
      Unweighted standardizer:   132 132
      Single group standardizer: 141 141

# size.ci.ratio.mean2 returns valid matrix

    Code
      res
    Output
       n1  n2
       53 106

# size.ci.lc.mean.bs returns valid number

    Code
      res
    Output
       Sample size per group
                          34

# size.ci.stdmean.ps returns valid number

    Code
      res
    Output
       Sample size
                19

# size.ci.ratio.mean2 returns valid number

    Code
      res
    Output
                                 Sample size
      Unweighted standardizer:            46
      Single group standardizer:          52

# size.ci.ratio.mean.ps returns valid number

    Code
      res
    Output
       Sample size
                21

# size.ci.lc.stdmean.ws returns valid matrix

    Code
      res
    Output
       Sample size
                11

# size.ci.lc.mean.ws returns valid matrix

    Code
      res
    Output
                                 Sample size
      Unweighted standardizer:            26
      Single level standardizer:          35

# size.ci.cronbach returns valid number

    Code
      res
    Output
       Sample size
                89

# size.ci.second returns valid number

    Code
      res
    Output
       Second-stage sample size
                             70

# size.test.mean returns valid number

    Code
      res
    Output
       Sample size
                20

# size.test.mean2 returns valid matrix

    Code
      res
    Output
       n1 n2
       27 27

# size.test.lc.mean.bs returns valid matrix

    Code
      res
    Output
       Sample size per group
                          47

# size.equiv.mean2 returns valid matrix

    Code
      res
    Output
       Sample size per group
                          50

# size.supinf.mean2 returns valid matrix

    Code
      res
    Output
       Sample size per group
                         143

# size.test.mean.ps returns valid number

    Code
      res
    Output
       Sample size
                22

# size.test.lc.mean.ws returns valid matrix

    Code
      res
    Output
       Sample size
                29

# size.equiv.mean.ps returns valid number

    Code
      res
    Output
       Sample size
                68

# size.supinf.mean.ps returns valid number

    Code
      res
    Output
       Sample size
                38

# size.test.mann returns valid number

    Code
      res
    Output
       Sample size per group
                          44

# size.test.sign returns valid number

    Code
      res
    Output
       Sample size
                67

# size.test.sign.ps returns valid number

    Code
      res
    Output
       Sample size
                42

# size.test.cronbach returns valid number

    Code
      res
    Output
       Sample size
               139

# pi.score returns valid matrix

    Code
      res
    Output
       Predicted df       LL       UL
            24.5 39 17.02546 31.97454

# pi.score2 returns valid matrix

    Code
      res
    Output
                                   Predicted    df       LL       UL
      Equal Variances Assumed:         11.22 83.00 4.650454 17.78955
      Equal Variances Not Assumed:     11.22 72.34 4.569926 17.87007

# ci.var.upper returns valid number

    Code
      res
    Output
             UL
       17.23264

# etasqr.adj returns valid number

    Code
      res
    Output
       adj Eta-squared
                0.2824

# test.anova.bs returns valid matrix

    Code
      res
    Output
            F dfA dfE       p Eta-squared adj Eta-squared
       5.9196   2  57 0.00461       0.172          0.1429

# etasqr.gen.2way returns valid matrix

    Code
      res
    Output
                                                 A         B        AB
      A treatment, B classification:      0.300000 0.5435540 0.1811847
      A classification, B treatment:      0.484252 0.3804878 0.2047244
      A classification, B classification: 0.300000 0.3804878 0.1268293

# ci.ratio.cod2 returns valid matrix

    Code
      res
    Output
            COD1      COD2 COD1/COD2       LL       UL
       0.1333333 0.1232558  1.081761 0.494964 2.282254

# ci.etasqr returns valid matrix

    Code
      res
    Output
       Eta-squared adj Eta-squared      SE    LL     UL
             0.241          0.2214 0.06258 0.104 0.3493

# ci.reliability returns valid vector

    Code
      res
    Output
       Estimate    LL     UL
           0.88 0.849 0.9066

# ci.sign returns valid matrix

    Code
      res
    Output
       Estimate        SE        LL        UL
       0.826087 0.0790342 0.6711828 0.9809911

# ci.slope.mean.bs returns valid matrix

    Code
      res
    Output
                                    Estimate         SE      t    df     p        LL
      Equal Variances Assumed:     0.3664407 0.06770529 5.4123 36.00 0e+00 0.2291280
      Equal Variances Not Assumed: 0.3664407 0.07336289 4.9949 18.66 8e-05 0.2127008
                                          UL
      Equal Variances Assumed:     0.5037534
      Equal Variances Not Assumed: 0.5201806

# test.mono.mean.bs returns valid matrix

    Code
      res
    Output
           Estimate       SE        LL         UL
       1 2   -11.71 4.139530 -22.07803 -1.3419744
       2 3   -11.72 4.399497 -22.74731 -0.6926939
       3 4   -16.92 4.730817 -28.76921 -5.0707936

# power.mean returns valid matrix

    Code
      res
    Output
           Power
       0.8021669

# power.mean2 returns valid matrix

    Code
      res
    Output
           Power
       0.8398417

# power.mean.ps returns valid matrix

    Code
      res
    Output
           Power
       0.9074354

# power.lc.bs returns valid matrix

    Code
      res
    Output
           Power
       0.7221171

# ci.cqv returns valid matrix

    Code
      res
    Output
       Estimate        SE        LL        UL
            0.5 0.1552485 0.2617885 0.8841821

# ci.ratio.sd2 returns valid matrix

    Code
      res
    Output
            SD1      SD2   SD1/SD2       LL       UL
       5.711587 6.450667 0.8854257 0.486279 1.728396

# size.ci.etasqr returns valid matrix

    Code
      res
    Output
       Sample size per group
                          63

# ci.2x2.stdmean.bs returns valid matrix

    Code
      res
    Output
               Estimate adj Estimate      SE      LL      UL
      AB:       -1.4498      -1.4194 0.68852 -2.7992 -0.1003
      A:         0.4690       0.4592 0.33795 -0.1933  1.1314
      B:        -0.7533      -0.7375 0.34512 -1.4297 -0.0769
      A at b1:  -0.2558      -0.2505 0.46402 -1.1653  0.6536
      A at b2:   1.1939       1.1689 0.50014  0.2137  2.1742
      B at a1:  -1.4782      -1.4472 0.49284 -2.4441 -0.5122
      B at a2:  -0.0284      -0.0278 0.48204 -0.9732  0.9163

# ci.2x2.median.bs returns valid matrix

    Code
      res
    Output
               Estimate       SE         LL         UL
      AB:          -5.0 3.389735 -11.643758 1.64375833
      A:            1.5 1.694867  -1.821879 4.82187916
      B:           -2.0 1.694867  -5.321879 1.32187916
      A at b1:     -1.0 2.152661  -5.219138 3.21913797
      A at b2:      4.0 2.618464  -1.132095 9.13209504
      B at a1:     -4.5 2.311542  -9.030539 0.03053939
      B at a2:      0.5 2.479330  -4.359397 5.35939682

# ci.2x2.stdmean.ws returns valid matrix

    Code
      res
    Output
               Estimate adj Estimate      SE      LL     UL
      AB:        0.1725       0.1645 0.13655 -0.0951 0.4401
      A:         0.1092       0.1042 0.05753 -0.0035 0.2220
      B:         0.0747       0.0713 0.05921 -0.0413 0.1908
      A at b1:   0.1955       0.1864 0.08461  0.0297 0.3613
      A at b2:   0.0230       0.0219 0.09372 -0.1607 0.2067
      B at a1:   0.1610       0.1535 0.09457 -0.0244 0.3463
      B at a2:  -0.0115      -0.0110 0.08596 -0.1800 0.1570

# ci.2x2.stdmean.mixed returns valid matrix

    Code
      res
    Output
               Estimate adj Estimate      SE      LL     UL
      AB:       -1.9515      -1.8014 1.02687 -3.9642 0.0611
      A:         1.9091       1.8203 0.51904  0.8918 2.9264
      B:         1.0606       0.9790 0.37117  0.3331 1.7881
      A at b1:   0.9333       0.8348 0.60281 -0.2481 2.1148
      A at b2:   2.8849       2.5803 0.83966  1.2392 4.5306
      B at a1:   0.0848       0.0783 0.54692 -0.9871 1.1568
      B at a2:   2.0364       1.8797 0.71133  0.6422 3.4306

# ci.2x2.median.mixed returns valid matrix

    Code
      res
    Output
               Estimate       SE         LL        UL
      AB:         -3.50 2.681476 -8.7555973  1.755597
      A:           1.75 1.340738 -0.8777986  4.377799
      B:           4.25 1.028850  2.2334918  6.266508
      A at b1:     2.50 1.680868 -0.7944405  5.794441
      A at b2:     6.00 2.089258  1.9051295 10.094870
      B at a1:     0.00 1.313181 -2.5737873  2.573787
      B at a2:     3.50 1.996942 -0.4139342  7.413934

# ci.2x2.median.w returns valid matrix

    Code
      res
    Output
               Estimate        SE         LL       UL
      AB:          2.50 21.050122 -38.757482 43.75748
      A:          24.75  9.603490   5.927505 43.57250
      B:          18.25  9.101881   0.410641 36.08936
      A at b1:    26.00 11.813742   2.845491 49.15451
      A at b2:    23.50 16.323093  -8.492675 55.49267
      B at a1:    19.50 15.710347 -11.291715 50.29171
      B at a2:    17.00 11.850202  -6.225970 40.22597

# pi.var example

    Code
      res
    Output
            UL
       23.9724

# ci.bayes.normal returns valid matrix

    Code
      res
    Output
       Posterior mean Posterior SD       LL       UL
              24.9226    0.5543895 23.83602 26.00919

# spearmanbrown returns valid matrix

    Code
      res
    Output
       Reliability of r2 measurements
                                 0.75

# ci.mean.fpc returns valid matrix

    Code
      res
    Output
       Estimate        SE       LL       UL
           24.5 0.5381631 23.41146 25.58854

# test.mean returns valid number

    Code
      res
    Output
            t df       p
       2.5991 39 0.01313

# size.ci.mean.prior returns valid matrix

    Code
      res
    Output
       Sample size
                44

# size.ci.cv example

    Code
      res
    Output
       Sample size
                60

# size.ci.median example

    Code
      res
    Output
       Sample size
                64

# size.ci.median2 example

    Code
      res
    Output
       n1 n2
       72 72

# size.ci.lc.median.bs example

    Code
      res
    Output
       Sample size per group
                          51

# size.ci.sd example

    Code
      res
    Output
       Sample size
                49

# ci.mean.gen(.05, y) example

    Code
      res
    Output
                  Estimate       SE       LL       UL
      Square-root 35.79395 8.141122 18.64498 58.48619
      Geometric   33.26410 7.633328 19.47124 56.82741
      Harmonic    29.07073 7.978149 19.35236 58.39602

# ci.sd(.05, 4.65, 50) example

    Code
      res
    Output
       Estimate       LL      UL
           4.65 3.884303 5.79452

# ci.lc.mean.scheffe example

    Code
      res
    Output
       Estimate       SE       t       p        LL        UL
          -5.35 1.275231 -4.1953 0.00228 -9.089451 -1.610549

