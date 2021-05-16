---
title: "Genome-Environment-Geographic Association"
output:
  html_document:
    df_print: paged
    code_folding: "hide"
    theme: readable
author: Mathew Rees
date: '`r format(Sys.time(), "%d %B, %Y")`'
---


## Here I use genetic, ecological and geographical pairwise distance matrices to investigate pattern of diversification in Paubrasil.


```{r message=FALSE, warning=FALSE}
library(psych)    # Used to investigate correlations among predictors
library(vegan)    # Used to run RDA
library(remotes)
library(vcfR)
library(raster)
library(corrplot)
library(PerformanceAnalytics)
library(ade4)
library(ape)
library(geosphere)
library(tidyverse)
library(ggpubr)
library(rcompanion)
```


Load the data.


```{r}
# first we select only wild samples
samples <- read.csv("RDA.Chesla.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# we use our previously computed pairwise genetic distance matrix between all individuals following a custom script provided by Isaac Overcast.
genetic.dist <- as.matrix(read.csv("distances.wild.csv", row.names = 1))

# we compute the geographic distance matrix
x <- cbind(samples$Long_DD, samples$Lat_DD)
d.geo <- distm(x, fun = distGeo)
geo.dist<-as.dist(d.geo)
str(geo.dist)

# we format the distance matrix as a distance object
gen.dist <- as.dist(genetic.dist)
str(gen.dist)

# We compute the ecological distance matrix and scale all variables
x <- scale(samples[,c(4:22)]) # scale(samples[,c(5,7,8, 15,18)])
eco.dist <- dist(x)
str(eco.dist)
```


Now lets visualize these correlations to see if we have a linear relationships between our co-variants.

plot regression lines of distance matrices:
https://jkzorz.github.io/2019/07/09/scatter-plots.html


```{r}
aa = as.vector(eco.dist)
tt = as.vector(gen.dist)
gg = as.vector(geo.dist)

#new data frame with vectorized distance matrices
mat = data.frame(aa,tt,gg)
```

 
Let's check the relationship between ecological and geographic distance.


```{r}
mm = ggplot(mat, aes(y = aa, x = gg/1000)) + 
    geom_point(size = 4, alpha = 0.5) + 
    geom_smooth(method = "lm", colour = "red", alpha = 0.4) +
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~~~"))) +
    labs(x = "Physical separation (in km)", y = "Genetic distance") + 
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))
mm
```


Not quite the best fit but we can see a linear trend here.


```{r}
lmaa.gg<-lm(aa~gg)
summary(lmaa.gg)
plot(lmaa.gg)
```


Now let's plot genetic distance vs ecological distance.


```{r}
mm = ggplot(mat, aes(y = tt, x = aa)) + 
    geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21, aes(fill = gg/1000)) + 
    #geom_smooth(method = "lm", colour = "red", alpha = 0.5) + 
    stat_smooth() +
    #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~~~"))) +
    labs(x = "Ecological distance Distance", y = "Genetic distance", fill = "Physical Separation (km)") + 
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
           axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
           axis.title= element_text(face = "bold", size = 14, colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.position = "right",
           legend.text = element_text(size = 10, face = "bold"),
           legend.title = element_text(size = 11, face = "bold")) +
    scale_fill_continuous(high = "salmon", low = "navy")
    
mm
```


The data doesn't look normally distributed with a large gap in genetic distances. This is because we are comparing intra-population (bottom left of the graph) level with inter-population level (top of graph) at the same time.

Let's see with Genetic distance vs geographical distance.


```{r}
#genetic vs geographic distance
mm = ggplot(mat, aes(y = tt, x = gg/1000)) + 
    geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21, aes(fill = aa)) + 
    #geom_smooth(method = "lm", colour = "red", alpha = 0.5) + 
    stat_smooth() +
    #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~~~"))) +
    labs(x = "Physical distance (in km)", y = "Genetic distance", fill = "Ecological distance") + 
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
           axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
           axis.title= element_text(face = "bold", size = 14, colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.position = "right",
           legend.text = element_text(size = 10, face = "bold"),
           legend.title = element_text(size = 11, face = "bold")) +
    scale_fill_continuous(high = "salmon", low = "navy")
    
mm
```


Very similar patterns. Here we find that genetic distance follow a non-linear relationship with both ecological and geographical distance, with strong increase in genetic distance with the first few km but then flattens out. The linear model is actually a poor fit for this relationship.

Let's check the distribution of each variable.


```{r}
hist(aa, main="histogram of ecological distance")
hist(tt, main="histogram of genetic distance")
hist(gg, main="histogram of geographic distance")
```


Next we will remove all 0 values in geographic and ecological distance matrices (Ie. values from the same populations) and remove all genetic distances < 0.3 to remove interpopulation measures and avoid this bimodal distribution.


```{r, warning=FALSE}
aa = as.vector(eco.dist)
tt = as.vector(gen.dist)
gg = as.vector(geo.dist)

# to look at interpopulation variation, we will take out all 0 values in the ecological and geographic distance matrices

aa[aa==0]<-NA
tt[tt<0.3]<-NA  # here 0.3 is an arbitrary value given the bimodal distribution. 
                # It represents the cutoff point between variation within clusters and between clusters.
gg[gg==0]<-NA

#new data frame with vectorized distance matrices
mat = data.frame(aa,tt,gg)
hist(aa, main="histogram of ecological distance")
hist(gg/1000, main="histogram of geographic distance")
```


Lets look at what we have now, between genetic clusters. We plot genetic distance vs ecological distance.


```{r, warning=FALSE}
mm = ggplot(mat, aes(y = tt, x = aa)) + 
    geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21, aes(fill = gg/1000)) + 
    stat_smooth() + ## loess in blue
    labs(x = "Ecological Distance", y = "Genetic distance", fill = "Physical Separation") + 
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
           axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
           axis.title= element_text(face = "bold", size = 14, colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.position = "right",
           legend.text = element_text(size = 10, face = "bold"),
           legend.title = element_text(size = 11, face = "bold")) +
    scale_fill_continuous(high = "salmon", low = "navy")
    
mm
```


This looks pretty random with no obvious pattern. But wait, what happens it we flip this the other way round?


```{r, warning=FALSE}
# Fit polynomial regression line and add labels
formula <- y ~ poly(x, 5, raw = TRUE)
mm = ggplot(mat, aes(y = aa, x = tt)) + 
    geom_point(size = 4, alpha = 0.75, colour = "black",shape = 21, aes(fill = gg/1000)) + 
    geom_smooth(method = "lm", formula = formula, colour = "red", alpha = 0.5) + 
    #stat_smooth() +
    #stat_cor(aes(label =  paste(..rr.label.., ..p.label.., sep = "~~~~"))) +
    labs(x = "Genetic distance", y = "Ecological Distance", fill = "Physical Separation") + 
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
           axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
           axis.title= element_text(face = "bold", size = 14, colour = "black"), 
           panel.background = element_blank(), 
           panel.border = element_rect(fill = NA, colour = "black"),
           legend.position = "right",
           legend.text = element_text(size = 10, face = "bold"),
           legend.title = element_text(size = 11, face = "bold")) +
    scale_fill_continuous(high = "salmon", low = "navy")
    
mm
```


The pattern is much more complex, and a linear relationship does not fit this model. Rather a quintic polynomial  seems to fit better.

Let's try the same thing between genetic and geographic distances.


```{r, warning=FALSE}
#genetic vs geographic distance
mm = ggplot(mat, aes(y = tt, x = gg/1000)) + 
    geom_point(size = 3, alpha = 0.75, colour = "black",shape = 21, aes(fill = aa)) + 
    #geom_smooth(method = "lm", formula = formula, colour = "red", alpha = 0.4) +
    stat_smooth() +
    #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~~~"))) +
    labs(x = "Physical separation (in km)", y = "Genetic distance", fill = "Ecological distance") + 
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black")) +
    scale_fill_continuous(high = "salmon", low = "navy")
mm
```


Flip it round.


```{r, warning=FALSE}
#genetic vs geographic distance
formula <- y ~ poly(x, 2, raw = TRUE)
mm = ggplot(mat, aes(y = gg/1000, x = tt)) + 
    geom_point(size = 3, alpha = 0.75, colour = "black",shape = 21, aes(fill = aa)) + 
    geom_smooth(method = "lm", formula = formula, colour = "red", alpha = 0.4) +
    stat_smooth() +
    #stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~~~"))) +
    labs(x = "Genetic distance", y = "Physical separation", fill = "Ecological distance") + 
    theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
        axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
        axis.title= element_text(face = "bold", size = 14, colour = "black"), 
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black")) +
    scale_fill_continuous(high = "salmon", low = "navy")
mm
```


Now we try a few models to check what might be the best fit for our data.

First, between ecological and genetic distance


```{r}
# Between ecological and genetic distance
test <- lm(aa~tt)
test2 <- lm(aa~poly(tt, 2, raw = T))
test3 <- lm(aa~poly(tt, 3, raw = T))
test4 <- lm(aa~poly(tt, 4, raw = T))
test5 <- lm(aa~poly(tt, 5, raw = T))
test6 <- lm(aa~poly(tt, 6, raw = T))
compareLM(test, test2,test3,test4,test5,test6)
anova(test,test2,test3,test4,test5,test6)
summary(test5)
```


Then between geographic and genetic distance.


```{r}
#between geographic and genetic distance 
test <- lm(gg~tt)
test2 <- lm(gg~poly(tt, 2, raw = T))
test3 <- lm(gg~poly(tt, 3, raw = T))
test4 <- lm(gg~poly(tt, 4, raw = T))
test5 <- lm(gg~poly(tt, 5, raw = T))
test6 <- lm(gg~poly(tt, 6, raw = T))
compareLM(test, test2,test3,test4,test5,test6)
anova(test,test2,test3,test4,test5,test6)
summary(test2)
```


The end.

