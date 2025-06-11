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
       Estimate adj Estimate        SE        LL       UL
       1.232877     1.209015 0.2124335 0.8165146 1.649239

# ci.mean2 returns valid matrix

    Code
      res
    Output
                                   Estimate        SE        t       df            p
      Equal Variances Assumed:          5.1 0.7151214 7.131656 48.00000 4.621279e-09
      Equal Variances Not Assumed:      5.1 0.6846568 7.448987 46.17476 1.898214e-09
                                         LL       UL
      Equal Variances Assumed:     3.662152 6.537848
      Equal Variances Not Assumed: 3.721998 6.478002

# ci.lc.mean.bs returns valid matrix

    Code
      res
    Output
                                   Estimate       SE         t       df            p
      Equal Variances Assumed:        -5.35 1.300136 -4.114955 36.00000 0.0002152581
      Equal Variances Not Assumed:    -5.35 1.300136 -4.114955 33.52169 0.0002372436
                                          LL        UL
      Equal Variances Assumed:     -7.986797 -2.713203
      Equal Variances Not Assumed: -7.993583 -2.706417

# ci.tukey returns valid matrix

    Code
      res
    Output
           Estimate       SE          t       df            p         LL         UL
       1 2    -4.71 1.142459  -4.122686 36.20303 5.989090e-04  -7.501842  -1.918158
       1 3   -13.43 1.104078 -12.163998 36.95915 0.000000e+00 -16.125713 -10.734287
       2 3    -8.72 1.228730  -7.096758 37.87646 5.510492e-08 -11.717053  -5.722947

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
                               Estimate adj Estimate        SE        LL       UL
      Unweighted standardizer: 1.174493     1.159240 0.2844012 0.6170771 1.731909
      Weighted standardizer:   1.174493     1.159240 0.2802826 0.6251494 1.723837
      Group 1 standardizer:    1.147541     1.117605 0.2975582 0.5643375 1.730744
      Group 2 standardizer:    1.203438     1.172044 0.3120525 0.5918268 1.815050

# ci.stdmean.strat returns valid matrix

    Code
      res
    Output
                                Estimate adj Estimate         SE         LL        UL
      Weighted standardizer: -0.05538872  -0.05528428 0.10023259 -0.2518410 0.1410636
      Group 1 standardizer:  -0.05714286  -0.05692722 0.10368609 -0.2603639 0.1460782
      Group 2 standardizer:  -0.05357143  -0.05692722 0.09720571 -0.2440911 0.1369483

# ci.lc.stdmean.bs returns valid matrix

    Code
      res
    Output
                                Estimate adj Estimate        SE        LL         UL
      Unweighted standardizer: -1.301263    -1.273964 0.3692800 -2.025039 -0.5774878
      Weighted standardizer:   -1.301263    -1.273964 0.3514511 -1.990095 -0.6124317
      Group 1 standardizer:    -1.393229    -1.273810 0.4849842 -2.343781 -0.4426775

# ci.mean.ps returns valid matrix

    Code
      res
    Output
       Estimate       SE        t df           p       LL       UL
            6.8 1.455922 4.670578 29 6.33208e-05 3.822304 9.777696

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
                                   Estimate adj Estimate        SE        LL
      Unweighted standardizer:    0.5550319    0.5433457 0.1609934 0.2394905
      Measurement 1 standardizer: 0.5424837    0.5253526 0.1615500 0.2258515
      Measurement 2 standardizer: 0.5684932    0.5505407 0.1692955 0.2366800
                                         UL
      Unweighted standardizer:    0.8705732
      Measurement 1 standardizer: 0.8591158
      Measurement 2 standardizer: 0.9003063

# ci.lc.stdmean.ws returns valid matrix

    Code
      res
    Output
                                Estimate adj Estimate        SE        LL         UL
      Unweighted standardizer: -1.301263    -1.266557 0.3147937 -1.918248 -0.6842788
      Level 1 standardizer:    -1.393229    -1.337500 0.3661824 -2.110934 -0.6755248

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
       Estimate         SE        LL        UL
           0.85 0.02456518 0.7971254 0.8931436

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
                                   Predicted       df       LL       UL
      Equal Variances Assumed:         11.22 83.00000 4.650454 17.78955
      Equal Variances Not Assumed:     11.22 72.34319 4.603642 17.83636

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
              0.282381

# test.anova.bs returns valid matrix

    Code
      res
    Output
              F dfA dfE           p Eta-squared adj Eta-squared
       5.919585   2  57 0.004614428   0.1719831       0.1429298

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
       Eta-squared adj Eta-squared         SE        LL        UL
             0.241       0.2213707 0.06258283 0.1040229 0.3493431

# ci.reliability returns valid vector

    Code
      res
    Output
       Estimate        LL        UL
           0.88 0.8489612 0.9065575

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
                                    Estimate         SE        t       df
      Equal Variances Assumed:     0.3664407 0.06770529 5.412290 36.00000
      Equal Variances Not Assumed: 0.3664407 0.07336289 4.994905 18.65826
                                              p        LL        UL
      Equal Variances Assumed:     4.242080e-06 0.2291280 0.5037534
      Equal Variances Not Assumed: 8.468223e-05 0.2126998 0.5201815

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
                  Estimate adj Estimate        SE         LL         UL
      AB:      -1.44976487   -1.4193502 0.6885238 -2.7992468 -0.1002829
      A:        0.46904158    0.4592015 0.3379520 -0.1933321  1.1314153
      B:       -0.75330920   -0.7375055 0.3451209 -1.4297338 -0.0768846
      A at b1: -0.25584086   -0.2504736 0.4640186 -1.1653006  0.6536189
      A at b2:  1.19392401    1.1688767 0.5001423  0.2136630  2.1741850
      B at a1: -1.47819163   -1.4471806 0.4928386 -2.4441376 -0.5122457
      B at a2: -0.02842676   -0.0278304 0.4820369 -0.9732017  0.9163482

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
                  Estimate adj Estimate         SE           LL        UL
      AB:       0.17248839   0.16446123 0.13654635 -0.095137544 0.4401143
      A:        0.10924265   0.10415878 0.05752822 -0.003510596 0.2219959
      B:        0.07474497   0.07126653 0.05920554 -0.041295751 0.1907857
      A at b1:  0.19548684   0.18638939 0.08460680  0.029660560 0.3613131
      A at b2:  0.02299845   0.02192816 0.09371838 -0.160686202 0.2066831
      B at a1:  0.16098916   0.15349715 0.09457347 -0.024371434 0.3463498
      B at a2: -0.01149923  -0.01096408 0.08595873 -0.179975237 0.1569768

# ci.2x2.stdmean.mixed returns valid matrix

    Code
      res
    Output
                  Estimate adj Estimate        SE         LL         UL
      AB:      -1.95153666  -1.80141845 1.0268728 -3.9641704 0.06109706
      A:        1.90911195   1.82026682 0.5190413  0.8918096 2.92641425
      B:        1.06061775   0.97903177 0.3711681  0.3331416 1.78809392
      A at b1:  0.93334362   0.83480791 0.6028071 -0.2481367 2.11482389
      A at b2:  2.88488027   2.58031536 0.8396553  1.2391862 4.53057438
      B at a1:  0.08484942   0.07832254 0.5469232 -0.9871003 1.15679910
      B at a2:  2.03638608   1.87974099 0.7113341  0.6421968 3.43057538

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
              t df          p
       2.599132 39 0.01312665

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

