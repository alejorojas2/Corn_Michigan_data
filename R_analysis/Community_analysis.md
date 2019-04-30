---
title: "Community analysis"
output:
  html_document:
    keep_md: yes
editor_options:
  chunk_output_type: console
---


```
##     phyloseq      ggplot2 RColorBrewer         plyr    tidyverse 
##         TRUE         TRUE         TRUE         TRUE         TRUE 
##        knitr     magrittr          ape        vegan       ggtree 
##         TRUE         TRUE         TRUE         TRUE         TRUE 
##      cowplot      ggrepel       gtable    gridExtra         biom 
##         TRUE         TRUE         TRUE         TRUE         TRUE 
##       ampvis         here 
##         TRUE         TRUE
```

## Importing amplicon data

Read the `biom` files to do the analysis of the data limiting to the fields sampled in Michigan in 2011 and 2012.


```r
#File reading
Corn_COI_file <- here("data","COI_soil_otu_table_open.biom")
CornCOI <- import_biom(Corn_COI_file)
```

```
## Warning in strsplit(conditionMessage(e), "\n"): input string 1 is invalid
## in this locale
```

```r
#Filtering samples to Corn fields
Corn_COI <- prune_samples((grepl('MICO', sample_names(CornCOI))), CornCOI)

#Rename tax ranks to actual names
colnames(tax_table(Corn_COI)) <- c("Phylum","Class","Order",
                                  "Family","Genus","Clade","Species","Kingdom")

#Remove non-assigned taxa
Corn_COI <- subset_taxa(Corn_COI, Phylum != "No blast hit")
```


## Community composition based on amplicon data


```r
#Summarizing by field
Oom_phylo <- tax_glom(Corn_COI, taxrank = "Clade")
Oom_phylo <- merge_samples(Oom_phylo, "group")

#Pruning the top 20 taxa
most_taxa <- sort(taxa_sums(Oom_phylo), TRUE)[1:15]
Oom_phylo_top20 <- prune_taxa(names(most_taxa), Oom_phylo)
ntaxa(Oom_phylo_top20)
```

```
## [1] 15
```

```r
#Transform counts for plot
#Oom_phylo_state <- transform_sample_counts(Oom_phylo_state, 
#                                           function(x) 100 * x/sum(x))
State_Year <- psmelt(Oom_phylo_top20)

#Color scale
pal <- colorRampPalette(brewer.pal(12, "Paired"))

#Reorder factor
Clade_factor <- State_Year %>% group_by(Clade) %>% dplyr::summarise(sum(Abundance))
Clade_factor <- Clade_factor[order(-Clade_factor$`sum(Abundance)`),]
Clade_factor <- Clade_factor$Clade
State_Year$Clade <- factor(State_Year$Clade, levels = Clade_factor)
levels(State_Year$Clade)
```

```
##  [1] "Pythium_Clade_F"      "Pythium_Clade_I"      "Pythium_Clade_D"     
##  [4] "Pythium_Clade_B"      "Pythium_Clade_J"      "Pythium_Clade_E"     
##  [7] "Pythium_sp_nov"       "Phytopythium"         "Phytophthora_Clade_8"
## [10] "Phytophthora_Clade_1" "Pythium_Clade_H"      "Pythium_Clade_A"     
## [13] "Saprolegnia"          "Achlya"               "Aphanomyces"
```

```r
data_state <- dplyr::select(State_Year, Sample, Abundance, Clade)
data_state <- data_state[with(data_state, order(Clade, as.numeric(Clade))),]

#Plot
(ByState <- ggplot(data = data_state, aes(Sample, Abundance, fill = Clade)) +
  geom_bar(stat = "identity", position = position_fill()) + coord_flip() +
  scale_fill_manual(values = pal(20)) + 
  theme_gray() + 
  theme(text = element_text(size = 15)))
```

<img src="Community_analysis_files/figure-html/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />


## Culture data: isolation frequency

Initial analysis of culture data and plotting frequency of isolation per species


```r
#Importing otu like matrix, tax table and sample data for culture collection
culture_otu <- read.table(here("data","otu_corn_culture.txt")) %>% as.matrix()

culture_tax <- read.table(here("data","tax_corn_culture.txt")) %>% as.matrix() 

culture_data <- read.table(here("data","data_corn_culture.txt"), sep = "\t", header = TRUE, row.names = "group")

otu_corn_ct <- otu_table(culture_otu, taxa_are_rows = TRUE)
tax_corn_ct <- tax_table(culture_tax)
data_corn_ct <- sample_data(culture_data)

#Creating phyloseq object for culture collection
corn_ct <- phyloseq(otu_corn_ct, tax_corn_ct, data_corn_ct)

#Species recovery by isolation
corn_ct_sp <- psmelt(corn_ct)
corn_ct2 <- corn_ct_sp %>% group_by(Species, Year) %>% summarise(Abd = sum(Abundance))
corn_ct2$Year <- as.factor(corn_ct2$Year)
corn_ct2$Species <- gsub("_"," ", corn_ct2$Species)

(ct_bar_plot <- ggplot(corn_ct2) + 
  geom_bar(aes(x = reorder(Species, Abd, mean), y = Abd, fill = Year), stat = "identity") + 
  labs(x ="Species", y = "Frequency") +
  coord_flip() +
  scale_fill_manual(values = c("#7fbf7b", "#3288bd")) +
  theme_gray() +
  theme(text = element_text(size = 15),
        axis.text.y = element_text(face = "italic")))
```

<img src="Community_analysis_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

## Importing culture data and merging data with amplicon phyloseq object


```r
#merging amplicon dataset by field
Oom_group <- merge_samples(Corn_COI, "group")
Oom_group <- subset_taxa(Oom_group, Clade != "NA")

#adding info to phyloseq files
sample_data(corn_ct)$type <- "Culture"
sample_data(Oom_group)$type <- "COI amplicon"
sample_names(corn_ct) <- paste(sample_names(corn_ct),"_ct",sep = "")
sample_names(Oom_group) <- paste(sample_names(Oom_group),"_COI",sep = "")

#transform abundances
corn_ct2 <- transform_sample_counts(corn_ct, function(OTU) (OTU/sum(OTU)))
Oom_group2 <- transform_sample_counts(Oom_group, function(OTU) (OTU/sum(OTU)))

#merge and transform abundances
ct_COI_corn <- merge_phyloseq(Oom_group, corn_ct)
ct_COI_corn2 <- merge_phyloseq(Oom_group2, corn_ct2)
```

### Comparing culture data with amplicon data


```r
#Creating data frame with abundance data
ct_COI_corn.psmelt <- psmelt(ct_COI_corn)

ct_COI_corn_clade <- ct_COI_corn.psmelt %>% 
  group_by(Clade, type) %>% 
  summarise(abd = sum(Abundance)) %>%
  group_by(type) %>%
  mutate(total = sum(abd)) %>%
  mutate(Rel.abundance = (abd/total))


ct_COI_corn_clade <-  ct_COI_corn_clade %>%
  dplyr::select(Clade,type,Rel.abundance) %>% 
  spread(type, Rel.abundance) %>% 
  top_n(n =15, wt = `COI amplicon`) %>%
  gather(type, Rel.abundance, 2:3)
ct_COI_corn_clade$Clade <- gsub("_"," ", ct_COI_corn_clade$Clade)

(clade.graph <- ggplot(ct_COI_corn_clade) + 
  geom_bar(aes(x = reorder(Clade, Rel.abundance, mean), y = Rel.abundance), stat = "identity") + 
  facet_grid(. ~ type, scales = "free_y") +
  ylab("Relative abundance") +
  xlab("Genus/Clade") +
  coord_flip() + 
  theme_gray(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(face = "bold", size = 20)))
```

<img src="Community_analysis_files/figure-html/unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

```r
ct_COI_corn_genus <- psmelt(ct_COI_corn) %>% 
  group_by(Species, type) %>% 
  summarise(abd = sum(Abundance)) %>%
  group_by(type) %>%
  mutate(total = sum(abd)) %>%
  mutate(Rel.abundance = (abd/total))

ct_COI_corn_genus <-  ct_COI_corn_genus %>%
  dplyr::select(Species,type,Rel.abundance) %>% 
  spread(type, Rel.abundance) %>% 
  top_n(n=12, wt = `COI amplicon`) %>%
  gather(type, Rel.abundance, 2:3)
ct_COI_corn_genus$Species <- gsub("_"," ", ct_COI_corn_genus$Species)

(species.graph <- ggplot(ct_COI_corn_genus) + 
  geom_bar(aes(x = reorder(Species, Rel.abundance, mean), 
               y = Rel.abundance), stat = "identity") + 
  facet_grid(. ~ type, scales = "free_y") +
  xlab("Top 12 Species") +
  ylab("Relative abundance") +
  coord_flip() + 
  theme_gray(base_size = 16) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 18),
        strip.text = element_text(face = "bold", size = 20)))
```

<img src="Community_analysis_files/figure-html/unnamed-chunk-6-2.png" style="display: block; margin: auto;" />


## Alpha diversity for amplicon data


```r
#Plotting richness
alpha.corn <- plot_richness(Corn_COI, "group", 
              measures = c("InvSimpson","Shannon","Chao1"), color = "Year")

alpha.corn + geom_boxplot() + xlab("Field") + theme_gray() +
  theme(axis.text.x = element_text(angle = 90))
```

![](Community_analysis_files/figure-html/unnamed-chunk-7-1.png)<!-- -->

__Richness table__


```r
#Table richness
Tb_richness <- estimate_richness(Corn_COI, 
                                 split = TRUE, 
                                 c("Observed", "Shannon", "Simpson")) %>% 
  rownames_to_column(var = "sample") %>%
  mutate(Evenness = Shannon/log(Observed))

samp_size <- colSums(otu_table(Corn_COI))
samp_size <- data.frame(samp_size) %>% rownames_to_column(var = "sample")

smp_state <- data.frame(sample_data(Corn_COI)[,71:72]) %>% 
  rownames_to_column(var = "sample")

Tb_richness_final <- left_join(Tb_richness, samp_size, by = "sample") %>%
  left_join(smp_state, by = "sample") %>%
  filter(Observed > 1) %>%
  ddply(c("group"), summarise, 
        N = length(sample),
        #OTUs = sum(samp_size),
        mean.Observed = mean(Observed),
        sd.Observed = sd(Observed, na.rm = TRUE),
        mean.Shannon = mean(Shannon),
        sd.Shannon = sd(Shannon, na.rm = TRUE),
        mean.Simpson = mean(Simpson),
        sd.Simpson = sd(Simpson, na.rm = TRUE),
        mean.Evenness = mean(Shannon/log(Observed)),
        sd.Evenness = sd(Shannon/log(Observed), na.rm = TRUE))

kable(Tb_richness_final, digits = 3, format = "markdown")
```



|group   |  N| mean.Observed| sd.Observed| mean.Shannon| sd.Shannon| mean.Simpson| sd.Simpson| mean.Evenness| sd.Evenness|
|:-------|--:|-------------:|-----------:|------------:|----------:|------------:|----------:|-------------:|-----------:|
|MICO_2  |  3|       238.667|      84.240|        3.812|      0.170|        0.871|      0.008|         0.703|       0.020|
|MICO_3  |  3|       318.333|      38.553|        4.058|      0.302|        0.906|      0.018|         0.704|       0.039|
|MICO_4  |  3|       154.000|      43.093|        3.782|      0.355|        0.908|      0.033|         0.754|       0.048|
|MICO_5  |  3|       336.000|     124.527|        4.133|      0.127|        0.949|      0.009|         0.719|       0.056|
|MICO_6  |  3|       172.333|       6.506|        4.144|      0.031|        0.965|      0.004|         0.805|       0.004|
|MICO2_1 |  3|       399.000|      44.193|        4.343|      0.107|        0.939|      0.009|         0.726|       0.015|
|MICO2_2 |  3|       236.667|      99.229|        4.077|      0.008|        0.938|      0.024|         0.756|       0.051|
|MICO2_3 |  3|       226.333|      37.434|        4.514|      0.189|        0.968|      0.013|         0.834|       0.011|
|MICO2_4 |  3|       326.333|      61.614|        4.444|      0.228|        0.954|      0.026|         0.770|       0.046|
|MICO2_5 |  3|       545.000|      25.710|        4.799|      0.106|        0.958|      0.006|         0.762|       0.014|

### ANOSIM


```r
Yr_grp <- get_variable(Corn_COI, "Year")
(Corn_yr_anosim <- anosim(phyloseq::distance(Corn_COI, "bray"), Yr_grp))
```

```
## 
## Call:
## anosim(x = phyloseq::distance(Corn_COI, "bray"), grouping = Yr_grp) 
## Dissimilarity: bray 
## 
## ANOSIM statistic R: 0.09103 
##       Significance: 0.048 
## 
## Permutation: free
## Number of permutations: 999
```

```r
field_grp <- get_variable(Corn_COI, "group")
(Corn_grp_anosim <- anosim(phyloseq::distance(Corn_COI, "bray"), field_grp))
```

```
## 
## Call:
## anosim(x = phyloseq::distance(Corn_COI, "bray"), grouping = field_grp) 
## Dissimilarity: bray 
## 
## ANOSIM statistic R: 0.9893 
##       Significance: 0.001 
## 
## Permutation: free
## Number of permutations: 999
```


## Beta diversity

### ADONIS


```r
df <- as(sample_data(Corn_COI), "data.frame")
d <- phyloseq::distance(Corn_COI, "bray")

Corn_adonis <- adonis(d ~ Year + Year/group, df)

kable(Corn_adonis$aov.tab, digits = 3, 
      caption = "__Table .__ Comparison of community structure (beta diversity)\
      using Bray-curtis distance by year.", format = "markdown")
```



|           | Df| SumsOfSqs| MeanSqs| F.Model|    R2| Pr(>F)|
|:----------|--:|---------:|-------:|-------:|-----:|------:|
|Year       |  1|     0.761|   0.761|   5.038| 0.073|  0.001|
|Year:group |  8|     6.576|   0.822|   5.440| 0.635|  0.001|
|Residuals  | 20|     3.022|   0.151|      NA| 0.292|     NA|
|Total      | 29|    10.359|      NA|      NA| 1.000|     NA|


### Ordination analysis for culture and amplicon


```r
##Colors
colors2 <- c("#77C7C6",
"#C37F3B",
"#869BCF",
"#7DD54E",
"#C67BB7",
"#CDC84A",
"#CC6569",
"#83D693",
"#7A7678",
"#698547",
"#D3C1A7",
"#6ca556")

##Year as factor
sample_data(ct_COI_corn)$Year <- as.factor(sample_data(ct_COI_corn)$Year)

##Beta diversity culture and amplicon
Oom_biom_ord <- ordinate(ct_COI_corn, "PCoA", "bray")
ord_plot <- plot_ordination(ct_COI_corn, Oom_biom_ord, color = "Year", shape = "type")
(ord_plot.COI <- ord_plot + geom_point(size = 5, alpha = 0.7) + 
  scale_colour_manual(values = colors2) + labs(title = "coxI") + 
  theme_gray(base_size = 18) + labs(color = "Year", shape = "Type")) 
```

![](Community_analysis_files/figure-html/Ordination_plot-1.png)<!-- -->


## Beta diversity for amplicon data only


```r
#Extracting the top 100 taxa
most_taxa <- sort(taxa_sums(Corn_COI), TRUE)[1:100]
Corn_COI_top100 <- prune_taxa(names(most_taxa), Corn_COI)

#Ordination analysis
Oom_biom_ord2 <- ordinate(Corn_COI_top100, "PCoA", "bray")
ord_plot2 <- plot_ordination(Corn_COI_top100, Oom_biom_ord2, color = "group", shape = "Year")
(ord_plot.COI <- ord_plot2 + geom_point(size = 5, alpha = 0.7) + 
  scale_colour_manual(values = colors2) + labs(title = "PCoA COI") +
  theme_gray(base_size = 18) + labs(color = "group", shape = "Year") +
  theme(legend.text = element_text(size = 20),
        legend.key.size = unit(0.9, "cm")))
```

![](Community_analysis_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

### Environmental/Edaphic factor analysis


```r
## Environment fit analysis
bray.pcoa <- ordinate(Corn_COI_top100, method = "PCoA", "bray")
env <- sample_data(Corn_COI_top100)
env <- env[,c(1,6,7,8,10,11,15,16,17,18,69,73)]
env[,1:12] <- sapply(env[,1:12], as.numeric)
str(env)
```

```
## 'data.frame':	30 obs. of  12 variables:
## Formal class 'sample_data' [package "phyloseq"] with 4 slots
##   ..@ .Data    :List of 12
##   .. ..$ : num  0.153 0.18 0.153 0.283 0.283 ...
##   .. ..$ : num  13.7 8.7 13.7 99.4 99.4 ...
##   .. ..$ : num  25.5 26.8 25.5 20.6 20.6 ...
##   .. ..$ : num  1.57 1.61 1.57 1 1 ...
##   .. ..$ : num  0 0 0 0 0 0 0 0 0 0 ...
##   .. ..$ : num  0 0 0 0 0 0 0 0 0 0 ...
##   .. ..$ : num  0.625 0.54 0.625 42.623 42.623 ...
##   .. ..$ : num  1115 1151 1115 1085 1085 ...
##   .. ..$ : num  728 830 728 845 845 ...
##   .. ..$ : num  822 955 822 905 905 ...
##   .. ..$ : num  30 29.5 30 21.4 21.4 ...
##   .. ..$ : num  7.342 7 7.342 0.874 0.874 ...
##   ..@ names    : chr  "AWC" "CEC7" "Clay" "Db3rdbar" ...
##   ..@ row.names: chr  "MICO2.1.2012.R1.COI" "MICO.3.2011.R2.COI" "MICO2.1.2012.R3.COI" "MICO2.5.2012.R1.COI" ...
##   ..@ .S3Class : chr "data.frame"
```

```r
Oom_env <- envfit(bray.pcoa$vectors, env, permutations = 999, na.rm = TRUE)

fit_data <- as.data.frame(scores(Oom_env, display = "vectors")) %>%
  add_rownames(var = "Env.var") %>%
  bind_cols(data.frame(Oom_env$vectors$r, Oom_env$vectors$pvals)) %>%
  #rename(R2 = Oom_env.vectors.r, P.value = Oom_env.vectors.pvals) %>%
  arrange(Oom_env.vectors.pvals)
  
## Supplementary material version

kable(fit_data, digits = 3, caption = "Significance and correlation\
of vectors fitted into PCoA ordination of oomycete communities associated with\
Corn seedlings", format = "markdown")
```



|Env.var         | Axis.1| Axis.2| Oom_env.vectors.r| Oom_env.vectors.pvals|
|:---------------|------:|------:|-----------------:|---------------------:|
|OrgMatter       | -0.408|  0.551|             0.471|                 0.001|
|WC3rdbar        |  0.811|  0.153|             0.682|                 0.001|
|pHwater         |  0.513| -0.541|             0.557|                 0.001|
|Clay            |  0.727|  0.101|             0.539|                 0.002|
|Db3rdbar        |  0.458| -0.427|             0.392|                 0.003|
|CEC7            | -0.340|  0.517|             0.383|                 0.004|
|AWC             | -0.041|  0.527|             0.280|                 0.016|
|EC              |  0.027| -0.482|             0.233|                 0.025|
|Precip_2011     | -0.450|  0.100|             0.212|                 0.038|
|Precip_2012     |  0.074|  0.386|             0.154|                 0.114|
|Precip_30yr_avg | -0.296|  0.150|             0.110|                 0.192|
|ECEC            |  0.000|  0.000|             0.000|                 1.000|

### Results ordination and environmental data


```r
## Vectors for plot
fit_reduced <- fit_data[fit_data$Oom_env.vectors.pvals < 0.05,] 

fit_plot <- as.data.frame(scores(Oom_env, display = "vectors")) %>%
  add_rownames(var = "Env.var") %>%
  inner_join(fit_reduced, by = "Env.var") %>%
  arrange(Oom_env.vectors.pvals) 

fit_plot$Env.var2 <- c("Clay (%)","Organic matter (%)", "Water content",
                       "Soil pH","Bulk density (g/cm^3)",
                       "CEC (meq/100g soil)",
                       "Avail. water capacity (cm water/cm soil)",
                       "EC (dS/m)",
                       "Precipitation (mm)")

ord_plot.data <- plot_ordination(Corn_COI_top100, Oom_biom_ord2, 
                            color = "group", shape = "Year", justDF = TRUE)

(ord.plot.env <- ggplot(data = ord_plot.data, aes(x = Axis.1, y = Axis.2)) + 
  geom_point(aes(color=group, shape=Year), size = 6, alpha = 0.9) + 
  #scale_color_brewer(type = "div", palette ="Spectral") +
  labs(color = "group", shape = "Year", x = "PCoA 1 [18.2%]", y = "PCoA 2 [10.4%]") +
  scale_colour_manual(values = colors2) +
  geom_segment(data = fit_plot, aes(x = 0, xend = Axis.1.x, y = 0, yend = Axis.2.x), 
               arrow = arrow(length = unit(0.1,"cm")), color = "black", size = 1) + 
  geom_label_repel(data = fit_plot, aes(x = Axis.1.x, y = Axis.2.x, label = Env.var2), 
            size = 5, force = 1) + #facet_wrap(~Year) +
  theme_gray(base_size = 18) + theme(legend.title = element_text(size = 22),
                                     legend.text = element_text(size = 20),
                                     legend.key.size = unit(0.8, "cm")))
```

<img src="Community_analysis_files/figure-html/Ord_envfit_plot-1.png" style="display: block; margin: auto;" />


## Relative abundance of top 20 species


```r
Corn_COI.f <- subset_taxa(Corn_COI, !is.na(Clade)) 
amp_heatmap(data = Corn_COI.f,
            group = "group",
            tax.show = 20,
            tax.aggregate = "Species",
            tax.add = "Clade")
```

<img src="Community_analysis_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

