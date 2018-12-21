---
title: "Oomycete pathogenicity on Corn - Michigan"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---


```
##      ggplot2 RColorBrewer         grid    gridExtra         plyr 
##         TRUE         TRUE         TRUE         TRUE         TRUE 
##         lme4      lsmeans        knitr        tidyr        dplyr 
##         TRUE         TRUE         TRUE         TRUE         TRUE 
##         MASS     magrittr     reshape2   FactoMineR      cowplot 
##         TRUE         TRUE         TRUE         TRUE         TRUE 
##      stringr         here      cowplot 
##         TRUE         TRUE         TRUE
```



```r
#Reading the file
corn_data <- read.csv(file = here("data","Corn_clean_data"), sep = "\t")
corn_data$species <- trim(corn_data$species)
corn_data$temp <- factor(corn_data$temp)
corn_data <-  corn_data %>% filter(species != "litorale")


#Summarizing data using different parameters by plyr library
corn_sum <- ddply(corn_data, c("species","temp"), summarise,
              N = length(plant.seed),
              mean_wp = mean(plant.seed), 
              sd_wp = sd(plant.seed),
              se_wp = sd_wp/sqrt(N),
              mean_ph = mean(avg.height), 
              sd_ph = sd(avg.height),
              se_ph = sd_ph/sqrt(N)
              )
  
#Setting limits for error bars
wp_limits <- aes(ymax = mean_wp + se_wp, ymin=mean_wp - se_wp)
ph_limits <- aes(ymax = mean_ph + se_ph, ymin=mean_ph - se_ph)
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.


```r
(plot_wp <- ggplot(corn_sum, 
                   aes(x = reorder(species, mean_wp, median), 
                       y = mean_wp, 
                       colour = as.factor(temp))) +
  geom_point(stat = "identity", size = 3) +
  geom_errorbar(wp_limits, width=0.2) + 
  scale_color_manual(values = c("#1f78b4", "#33a02c")) + theme_gray() +
  theme(axis.text.x=element_text(angle=60, hjust = 1, vjust = 1, face="italic"), 
        axis.text.y=element_text(angle=90, hjust = 0.5), 
        plot.margin=unit(c(1,1,1,1), "mm")) +
  labs(x="Species", y = "Weight per plant (g)", colour = "Temperature (ºC)"))
```

![](Pathogenicity_corn_files/figure-html/unnamed-chunk-2-1.png)<!-- -->





```r
(plot_ph <- ggplot(corn_sum, 
                   aes(x = reorder(species, mean_ph, mean), 
                       y = mean_ph, 
                       colour = as.factor(temp))) +
  geom_point(stat = "identity", size = 3) +
  geom_errorbar(ph_limits, width=0.2) + 
  scale_color_manual(values = c("#1f78b4", "#33a02c")) + 
   theme_gray() +
  theme(axis.text.x=element_text(angle=60, hjust = 1, vjust = 1, face="italic"), 
        axis.text.y=element_text(angle=90, hjust = 0.5), 
        plot.margin=unit(c(1,1,1,1), "mm")) +
  labs(x="Species", y = "Plant height (cm)", colour = "Temperature (ºC)"))
```

![](Pathogenicity_corn_files/figure-html/unnamed-chunk-3-1.png)<!-- -->


```r
hist(corn_data$avg.height)
```

![](Pathogenicity_corn_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

```r
hist(corn_data$plant.seed)
```

![](Pathogenicity_corn_files/figure-html/unnamed-chunk-4-2.png)<!-- -->

```r
corn_data$Ht <- sqrt(corn_data$avg.height)

qqnorm(corn_data$plant.seed)
```

![](Pathogenicity_corn_files/figure-html/unnamed-chunk-4-3.png)<!-- -->

```r
qqnorm(corn_data$Ht)
```

![](Pathogenicity_corn_files/figure-html/unnamed-chunk-4-4.png)<!-- -->
#Plant height 15C


```r
set_15 <- subset(corn_data, corn_data$temp=="15")
fit.ht <- nlme::lme(Ht ~ species, random = ~1|isolate, method="ML", data=set_15)
fit1.ht <- nlme::lme(Ht ~ species, random = ~1|set/isolate, method="ML", data=set_15)
fit2.ht <- nlme::lme(Ht ~ species, random = ~1|set, method="ML", data=set_15)
anova(fit.ht, fit1.ht, fit2.ht)
```

```
##         Model df      AIC      BIC    logLik   Test   L.Ratio p-value
## fit.ht      1 24 1214.638 1330.939 -583.3189                         
## fit1.ht     2 25 1213.162 1334.310 -581.5813 1 vs 2   3.47536  0.0623
## fit2.ht     3 24 1458.978 1575.279 -705.4889 2 vs 3 247.81532  <.0001
```

```r
anova(fit1.ht)
```

```
##             numDF denDF  F-value p-value
## (Intercept)     1   877 475.3891  <.0001
## species        21    40   3.2886   6e-04
```

```r
fit.ht.ls <- lsmeans(fit1.ht, "species")
(fit.ht.ctrl <- contrast(fit.ht.ls, "trt.vs.ctrl", ref=7))
```

```
##  contrast                         estimate        SE df t.ratio p.value
##  acanthicum - control           -0.3130649 0.2941399 40  -1.064  0.9379
##  aff. dissotocum - control      -0.2942624 0.2941399 40  -1.000  0.9533
##  aff. torulosum - control       -0.2425588 0.2795203 40  -0.868  0.9765
##  amasculinum - control          -0.3648563 0.3245827 40  -1.124  0.9207
##  arrhenomanes - control         -0.2418390 0.3964177 40  -0.610  0.9963
##  attrantheridium - control      -0.2069052 0.3964177 40  -0.522  0.9984
##  heterothallicum - control      -0.2875268 0.2967988 40  -0.969  0.9598
##  inflatum - control             -0.4852999 0.2818818 40  -1.722  0.6187
##  irregulare - control           -1.4870401 0.2795203 40  -5.320  0.0001
##  lutarium - control             -0.2900736 0.2793673 40  -1.038  0.9445
##  millet control - control       -0.2483335 0.3218482 40  -0.772  0.9871
##  orthogonon - control           -0.8245058 0.2967988 40  -2.778  0.1081
##  perillum - control             -0.3267643 0.3245827 40  -1.007  0.9519
##  perplexum - control            -0.2545727 0.2967988 40  -0.858  0.9778
##  pleroticum - control           -0.2937420 0.3218482 40  -0.913  0.9699
##  rostratifingens - control      -0.1895375 0.2941399 40  -0.644  0.9950
##  sansomeana - control           -0.3983720 0.2941399 40  -1.354  0.8297
##  sylvaticum - control           -0.3698670 0.2967988 40  -1.246  0.8774
##  tardicrescens - control        -0.6405635 0.2787287 40  -2.298  0.2779
##  ultimum - control              -0.8439331 0.3218482 40  -2.622  0.1504
##  ultimum var. ultimum - control -0.8550596 0.2694088 40  -3.174  0.0429
## 
## P value adjustment: dunnettx method for 21 tests
```

```r
ht.15 <- data.frame(summary(fit.ht.ctrl)) %>% 
              separate(contrast, c('species','contrast'), sep = " - ", remove = TRUE)
ht.15$ht.sg <- ifelse(ht.15$p.value<0.05,"SG","NS")
```

#Plant weight 15C

```r
set_15 <- subset(corn_data, corn_data$temp=="15")
fit <- nlme::lme(plant.seed ~ species, random = ~1|isolate, method="ML", data=set_15)
fit1 <- nlme::lme(plant.seed ~ species, random = ~1|trial/isolate, method="ML", data=set_15)
fit2 <- nlme::lme(plant.seed ~ species, random = ~1|set, method="ML", data=set_15)
anova(fit, fit1, fit2)
```

```
##      Model df       AIC        BIC   logLik   Test  L.Ratio p-value
## fit      1 24 -180.4110  -64.10986 114.2055                        
## fit1     2 25 -242.3353 -121.18826 146.1676 1 vs 2  63.9243  <.0001
## fit2     3 24  155.9854  272.28652 -53.9927 2 vs 3 400.3207  <.0001
```

```r
anova(fit2)
```

```
##             numDF denDF   F-value p-value
## (Intercept)     1   917 1725.7359  <.0001
## species        21   917   11.9142  <.0001
```

```r
fit.ls <- lsmeans(fit2, "species")
(fit.ctrl <- contrast(fit.ls, "trt.vs.ctrl", ref=7))
```

```
##  contrast                           estimate         SE  df t.ratio
##  acanthicum - control           -0.038888889 0.06112440 917  -0.636
##  aff. dissotocum - control      -0.109777778 0.06112440 917  -1.796
##  aff. torulosum - control       -0.142333333 0.05798770 917  -2.455
##  amasculinum - control          -0.152333334 0.06695843 917  -2.275
##  arrhenomanes - control         -0.224666667 0.08200699 917  -2.740
##  attrantheridium - control      -0.170000000 0.08200699 917  -2.073
##  heterothallicum - control      -0.221333334 0.06112440 917  -3.621
##  inflatum - control             -0.032333333 0.05798770 917  -0.558
##  irregulare - control           -0.476700000 0.05798770 917  -8.221
##  lutarium - control              0.132969697 0.05885974 917   2.259
##  millet control - control       -0.166333333 0.06695843 917  -2.484
##  orthogonon - control           -0.273422222 0.06112440 917  -4.473
##  perillum - control             -0.143000000 0.06695843 917  -2.136
##  perplexum - control            -0.103777778 0.06112440 917  -1.698
##  pleroticum - control           -0.009666667 0.06695843 917  -0.144
##  rostratifingens - control      -0.132888889 0.06112440 917  -2.174
##  sansomeana - control           -0.028444445 0.06112440 917  -0.465
##  sylvaticum - control           -0.102866667 0.06112440 917  -1.683
##  tardicrescens - control        -0.233433333 0.05798770 917  -4.026
##  ultimum - control              -0.071266667 0.06695843 917  -1.064
##  ultimum var. ultimum - control -0.257680000 0.05602144 917  -4.600
##  p.value
##   0.9958
##   0.5628
##   0.1814
##   0.2635
##   0.0914
##   0.3788
##   0.0058
##   0.9980
##   <.0001
##   0.2718
##   0.1699
##   0.0002
##   0.3407
##   0.6293
##   1.0000
##   0.3184
##   0.9993
##   0.6393
##   0.0012
##   0.9412
##   0.0001
## 
## P value adjustment: dunnettx method for 21 tests
```

```r
wp.15 <- data.frame(summary(fit.ctrl)) %>% 
              separate(contrast, c('species','contrast'), sep = " - ", remove = TRUE)
wp.15$wp.sg <- ifelse(wp.15$p.value<0.05,"SG","NS")
```

#Plant height 20C


```r
set_20 <- subset(corn_data, corn_data$temp=="20")
fit.ht.20 <- nlme::lme(Ht ~ species, random = ~1|isolate, method="ML", data=set_20)
fit1.ht.20 <- nlme::lme(Ht ~ species, random = ~1|set/isolate, method="ML", data=set_20)
fit2.ht.20 <- nlme::lme(Ht ~ species, random = ~1|set, method="ML", data=set_20)
anova(fit.ht.20, fit1.ht.20, fit2.ht.20)
```

```
##            Model df      AIC      BIC    logLik   Test  L.Ratio p-value
## fit.ht.20      1 24 1114.228 1220.541 -533.1138                        
## fit1.ht.20     2 25 1105.389 1216.132 -527.6944 1 vs 2 10.83895   0.001
## fit2.ht.20     3 24 1134.998 1241.311 -543.4989 2 vs 3 31.60902  <.0001
```

```r
anova(fit2.ht.20)
```

```
##             numDF denDF  F-value p-value
## (Intercept)     1   597 353.3978  <.0001
## species        21   597   7.4696  <.0001
```

```r
fit.ht.ls <- lsmeans(fit2.ht.20, "species")
(fit.ht.ctrl <- contrast(fit.ht.ls, "trt.vs.ctrl", ref=7))
```

```
##  contrast                         estimate        SE  df t.ratio p.value
##  acanthicum - control           -0.4957090 0.1701441 597  -2.913  0.0576
##  aff. dissotocum - control      -0.5899104 0.1701441 597  -3.467  0.0102
##  aff. torulosum - control       -0.6792656 0.1618258 597  -4.198  0.0006
##  amasculinum - control          -0.4327785 0.1884409 597  -2.297  0.2532
##  arrhenomanes - control         -0.6066991 0.2298262 597  -2.640  0.1184
##  attrantheridium - control      -0.5821832 0.2298262 597  -2.533  0.1526
##  heterothallicum - control      -0.5149865 0.1724528 597  -2.986  0.0468
##  inflatum - control             -0.8227359 0.1638748 597  -5.021  <.0001
##  irregulare - control           -1.4009593 0.1618258 597  -8.657  <.0001
##  lutarium - control             -0.4379220 0.1611370 597  -2.718  0.0974
##  millet control - control       -0.4467228 0.1860650 597  -2.401  0.2045
##  orthogonon - control           -1.0912978 0.1884409 597  -5.791  <.0001
##  perillum - control             -0.7988777 0.1884409 597  -4.239  0.0005
##  perplexum - control            -0.5133867 0.1724528 597  -2.977  0.0481
##  pleroticum - control           -0.4286719 0.1860650 597  -2.304  0.2496
##  rostratifingens - control      -0.4957247 0.1701441 597  -2.914  0.0576
##  sansomeana - control           -0.6185534 0.1701441 597  -3.635  0.0056
##  sylvaticum - control           -0.6779128 0.1724528 597  -3.931  0.0018
##  tardicrescens - control        -0.8875334 0.1611370 597  -5.508  <.0001
##  ultimum - control              -1.0584773 0.1860650 597  -5.689  <.0001
##  ultimum var. ultimum - control -1.0917047 0.1557874 597  -7.008  <.0001
## 
## P value adjustment: dunnettx method for 21 tests
```

```r
ht.20 <- data.frame(summary(fit.ht.ctrl)) %>% 
              separate(contrast, c('species','contrast'), sep = " - ", remove = TRUE)
ht.20$ht.sg <- ifelse(ht.20$p.value<0.05,"SG","NS")
```

#Plant weight 20C


```r
fit.20 <- nlme::lme(plant.seed ~ species, random = ~1|isolate, method="ML", data=set_20)
fit1.20 <- nlme::lme(plant.seed ~ species, random = ~1|trial/isolate, method="ML", data=set_20)
fit2.20 <- nlme::lme(plant.seed ~ species, random = ~1|set, method="ML", data=set_20)
anova(fit.20, fit1.20, fit2.20)
```

```
##         Model df       AIC       BIC   logLik   Test  L.Ratio p-value
## fit.20      1 24 -250.2763 -143.9630 149.1381                        
## fit1.20     2 25 -264.0266 -153.2836 157.0133 1 vs 2 15.75034   1e-04
## fit2.20     3 24 -207.7103 -101.3970 127.8551 2 vs 3 58.31635  <.0001
```

```r
anova(fit2.20)
```

```
##             numDF denDF  F-value p-value
## (Intercept)     1   597 4054.055  <.0001
## species        21   597    9.888  <.0001
```

```r
fit.ls <- lsmeans(fit2.20, "species")
(fit.ctrl <- contrast(fit.ls, "trt.vs.ctrl", ref=7, adjust ="bon"))
```

```
##  contrast                         estimate         SE  df t.ratio p.value
##  acanthicum - control           -0.1690000 0.05787062 597  -2.920  0.0762
##  aff. dissotocum - control      -0.3150000 0.05787062 597  -5.443  <.0001
##  aff. torulosum - control       -0.3255000 0.05490089 597  -5.929  <.0001
##  amasculinum - control          -0.3035000 0.06339409 597  -4.788  <.0001
##  arrhenomanes - control         -0.1960000 0.07764158 597  -2.524  0.2488
##  attrantheridium - control      -0.3610000 0.07764158 597  -4.650  0.0001
##  heterothallicum - control      -0.2573333 0.05787062 597  -4.447  0.0002
##  inflatum - control             -0.3580000 0.05490089 597  -6.521  <.0001
##  irregulare - control           -0.5383500 0.05490089 597  -9.806  <.0001
##  lutarium - control             -0.1662500 0.05490089 597  -3.028  0.0539
##  millet control - control       -0.1770000 0.06339409 597  -2.792  0.1135
##  orthogonon - control           -0.4395000 0.06339409 597  -6.933  <.0001
##  perillum - control             -0.3450000 0.06339409 597  -5.442  <.0001
##  perplexum - control            -0.3070000 0.05787062 597  -5.305  <.0001
##  pleroticum - control           -0.1930000 0.06339409 597  -3.044  0.0511
##  rostratifingens - control      -0.2140000 0.05787062 597  -3.698  0.0050
##  sansomeana - control           -0.2523333 0.05787062 597  -4.360  0.0003
##  sylvaticum - control           -0.3246667 0.05787062 597  -5.610  <.0001
##  tardicrescens - control        -0.3085000 0.05490089 597  -5.619  <.0001
##  ultimum - control              -0.4355000 0.06339409 597  -6.870  <.0001
##  ultimum var. ultimum - control -0.4446000 0.05303930 597  -8.382  <.0001
## 
## P value adjustment: bonferroni method for 21 tests
```

```r
wp.20 <- data.frame(summary(fit.ctrl)) %>% 
              separate(contrast, c('species','contrast'), sep = " - ", remove = TRUE)
wp.20$wp.sg <- ifelse(wp.20$p.value<0.05,"SG","NS")
```


```r
ht.15$temp <- 15
wp.15$temp <- 15
ht.20$temp <- 20
wp.20$temp <- 20

ht <- bind_rows(ht.15, ht.20) %>% dplyr::select(species,temp,ht.sg)
ht$temp <- as.factor(ht$temp)
wp <- bind_rows(wp.15, wp.20) %>% dplyr::select(species,temp, wp.sg)
wp$temp <- as.factor(wp$temp)

corn_final <- full_join(corn_sum, ht, c("species","temp"))
corn_final <- full_join(corn_final, wp, c("species","temp"))

corn_final[is.na(corn_final)] <- "NS"
corn_final
```

```
##                 species temp  N    mean_wp     sd_wp      se_wp   mean_ph
## 1            acanthicum   15 45 0.45377778 0.2351623 0.03505593  4.081333
## 2            acanthicum   20 30 0.64800000 0.1738321 0.03173725 10.722333
## 3       aff. dissotocum   15 45 0.38288889 0.2174064 0.03240904  4.140889
## 4       aff. dissotocum   20 30 0.50200000 0.1626982 0.02970448 10.197667
## 5        aff. torulosum   15 60 0.35033333 0.1978612 0.02554377  4.281333
## 6        aff. torulosum   20 40 0.49150000 0.2238784 0.03539828 10.081750
## 7           amasculinum   15 30 0.34033333 0.1439991 0.02629052  3.344333
## 8           amasculinum   20 20 0.51350000 0.2003845 0.04480734  9.300500
## 9          arrhenomanes   15 15 0.26800000 0.1321795 0.03412861  3.772667
## 10         arrhenomanes   20 10 0.62100000 0.1585209 0.05012872  8.192000
## 11      attrantheridium   15 15 0.32266667 0.1917762 0.04951639  3.877333
## 12      attrantheridium   20 10 0.45600000 0.1270346 0.04017185  8.409000
## 13              control   15 30 0.49266667 0.1535934 0.02804218  5.186000
## 14              control   20 20 0.81700000 0.1146390 0.02563406 13.504500
## 15      heterothallicum   15 45 0.27133333 0.1783969 0.02659384  3.598444
## 16      heterothallicum   20 30 0.55966667 0.1845541 0.03369482  8.685333
## 17             inflatum   15 60 0.46033333 0.1770298 0.02285445  3.545667
## 18             inflatum   20 40 0.45900000 0.1506635 0.02382199  9.762750
## 19           irregulare   15 60 0.01596667 0.3428496 0.04426169  1.240667
## 20           irregulare   20 40 0.27865000 0.3381409 0.05346477  6.296000
## 21             lutarium   15 55 0.62563636 0.4264034 0.05749622  3.962545
## 22             lutarium   20 40 0.65075000 0.1319127 0.02085723 10.666250
## 23       millet control   15 30 0.32633333 0.1836175 0.03352382  4.272000
## 24       millet control   20 20 0.64000000 0.1596509 0.03569903 10.445500
## 25           orthogonon   15 45 0.21924444 0.3345595 0.04987318  2.868667
## 26           orthogonon   20 20 0.37750000 0.3685344 0.08240681  8.733500
## 27             perillum   15 30 0.34966667 0.2498701 0.04561983  3.496667
## 28             perillum   20 20 0.47200000 0.1998249 0.04468221  7.145000
## 29            perplexum   15 45 0.38888889 0.1764166 0.02629864  3.753778
## 30            perplexum   20 30 0.51000000 0.1473443 0.02690127  8.742333
## 31           pleroticum   15 30 0.48300000 0.1445361 0.02638856  3.957000
## 32           pleroticum   20 20 0.62400000 0.1506058 0.03367648 10.721500
## 33      rostratifingens   15 45 0.35977778 0.2138594 0.03188028  4.507111
## 34      rostratifingens   20 30 0.60300000 0.1778676 0.03247404 10.699000
## 35           sansomeana   15 45 0.46422222 0.1663961 0.02480487  3.459778
## 36           sansomeana   20 30 0.56466667 0.1596134 0.02914129  9.075333
## 37           sylvaticum   15 45 0.38980000 0.2192860 0.03268923  3.449778
## 38           sylvaticum   20 30 0.49233333 0.1629625 0.02975275  7.797333
## 39        tardicrescens   15 60 0.25923333 0.3441159 0.04442517  2.894333
## 40        tardicrescens   20 40 0.50850000 0.2088614 0.03302388  8.002250
## 41              ultimum   15 30 0.42140000 0.3657861 0.06678310  2.919667
## 42              ultimum   20 20 0.38150000 0.2260560 0.05054766  7.379000
## 43 ultimum var. ultimum   15 75 0.23498667 0.3067846 0.03542444  2.439867
## 44 ultimum var. ultimum   20 50 0.37240000 0.2048475 0.02896981  7.259800
##       sd_ph     se_ph ht.sg wp.sg
## 1  1.842081 0.2746012    NS    NS
## 2  3.476357 0.6346931    NS    NS
## 3  1.853335 0.2762788    NS    NS
## 4  3.956340 0.7223255    SG    SG
## 5  1.346519 0.1738349    NS    NS
## 6  4.215524 0.6665328    SG    SG
## 7  1.076507 0.1965423    NS    NS
## 8  3.086873 0.6902458    NS    SG
## 9  1.166230 0.3011193    NS    NS
## 10 2.300806 0.7275786    NS    NS
## 11 1.049942 0.2710938    NS    NS
## 12 3.196976 1.0109725    NS    SG
## 13 1.702081 0.3107561    NS    NS
## 14 2.836314 0.6342192    NS    NS
## 15 1.141731 0.1701992    NS    SG
## 16 1.951123 0.3562247    SG    SG
## 17 1.416573 0.1828788    NS    NS
## 18 4.166901 0.6588449    SG    SG
## 19 1.758816 0.2270622    SG    SG
## 20 4.218188 0.6669540    SG    SG
## 21 1.458801 0.1967047    NS    NS
## 22 3.667107 0.5798205    NS    NS
## 23 2.280090 0.4162855    NS    NS
## 24 2.571151 0.5749268    NS    NS
## 25 2.236757 0.3334360    NS    SG
## 26 5.543239 1.2395059    SG    SG
## 27 1.341701 0.2449599    NS    NS
## 28 2.235250 0.4998171    SG    SG
## 29 1.280366 0.1908657    NS    NS
## 30 2.394490 0.4371721    SG    SG
## 31 1.488357 0.2717357    NS    NS
## 32 3.720822 0.8320012    NS    NS
## 33 1.626239 0.2424254    NS    NS
## 34 3.280959 0.5990184    NS    SG
## 35 1.473446 0.2196484    NS    NS
## 36 3.203792 0.5849297    SG    SG
## 37 1.963142 0.2926479    NS    NS
## 38 2.128051 0.3885271    SG    SG
## 39 1.794759 0.2317024    NS    SG
## 40 3.348569 0.5294553    SG    SG
## 41 3.031610 0.5534937    NS    NS
## 42 4.578454 1.0237734    SG    SG
## 43 1.952790 0.2254888    SG    SG
## 44 3.827308 0.5412631    SG    SG
```


```r
plot_wp1 <- ggplot(corn_final, 
                   aes(x = species, 
                       y = mean_wp, 
                       colour = as.factor(temp))) +
  geom_point(aes(alpha=wp.sg), stat = "identity", size = 3) +
  geom_errorbar(wp_limits, width=0.2) + 
  scale_color_manual(values = c("#1f78b4", "#33a02c")) + 
  scale_alpha_discrete(range = c(0.5, 0.9)) + theme_gray(base_size = 16) +
  theme(axis.text.x=element_text(angle=60, hjust = 1, vjust = 1, face="italic"), 
        axis.text.y=element_text(angle=90, hjust = 0.5), 
        plot.margin=unit(c(1,1,1,1), "mm"),
        legend.position = c(0.2 , 0.83),
        legend.box = "horizontal",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_blank()) +
  labs(x="Species", y = "Weight per plant (g)", colour = "Temperature (ºC)",
       alpha= "Species vs control")
```

```
## Warning: Using alpha for a discrete variable is not advised.
```

```r
spp <- corn_sum %>% filter(temp == 20) %>% 
                    arrange(mean_wp) %>% 
                    dplyr::select(species)
(spp <- spp$species)
```

```
##  [1] "irregulare"           "ultimum var. ultimum" "orthogonon"          
##  [4] "ultimum"              "attrantheridium"      "inflatum"            
##  [7] "perillum"             "aff. torulosum"       "sylvaticum"          
## [10] "aff. dissotocum"      "tardicrescens"        "perplexum"           
## [13] "amasculinum"          "heterothallicum"      "sansomeana"          
## [16] "rostratifingens"      "arrhenomanes"         "pleroticum"          
## [19] "millet control"       "acanthicum"           "lutarium"            
## [22] "control"
```

```r
plot_wp1$data$species <- factor(plot_wp1$data$species, levels = spp)
```



```r
plot_ph1 <- ggplot(corn_final, 
                   aes(x = reorder(species, mean_ph, mean), 
                       y = mean_ph, 
                       colour = as.factor(temp))) +
  geom_point(aes(alpha=ht.sg), stat = "identity", size = 3) +
  geom_errorbar(ph_limits, width=0.2) + 
  scale_color_manual(values = c("#1f78b4", "#33a02c")) + 
  scale_alpha_discrete(range = c(0.5, 0.9)) + theme_gray(base_size = 16) +
  theme(axis.text.x=element_text(angle=60, hjust = 1, vjust = 1, face="italic"), 
        axis.text.y=element_text(angle=90, hjust = 0.5), 
        plot.margin=unit(c(1,1,1,1), "mm"),
        legend.position = c(0.2 , 0.83),
        legend.box = "horizontal",
        legend.background = element_rect(fill="gray95"),
        legend.key = element_blank()) +
  labs(x="Species", y = "Plant height (cm)", 
       colour = "Temperature (ºC)", 
       alpha= "Species vs control")
```

```
## Warning: Using alpha for a discrete variable is not advised.
```



```r
plot_grid(plot_ph1, plot_wp1, labels = c("A","B"), nrow = 2, ncol = 1)
```

<img src="Pathogenicity_corn_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

