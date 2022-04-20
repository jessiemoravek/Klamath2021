### for Jessie
library(MARSS)
library(xtable)

# Data: read in data (from Jessie's file)
##Read in data
```{r echo = T, results = "hide", messages = F}
KSV_meantemps <- readRDS('KSV_meantemps.rds')
daily_means_long_klamath <- readRDS('daily_means_long.rds')
covariate_klamath <- readRDS('covariate.rds')

daily_means_long_klamath <-rbind(daily_means_long_klamath, KSV = KSV_meantemps)
str(daily_means_long_klamath)

```
rownames(transformed_dat_klamath)
matplot(t(transformed_dat_klamath[1:2,]), type="l") # The 2 in AP
matplot(t(transformed_dat_klamath[15:17,]), type="l") # The 2
matplot(t(transformed_dat_klamath[23:25,]), type="l") # The 2 i
matplot(t(transformed_dat_klamath[c(3,6,9,11),]), type="l") # The 2 i
matplot(t(transformed_dat_klamath[1:10,]), type="l") 
# matplot(t(daily_means_long_klamath[1:10,]), type="l") # not zscored


# Steps:
# 1) Among sensor replicates, drop time series with (many) gaps. 
# 2) Then take mean across replicate sensors, so that we end up with 1 time series per habitat, 2+9+1 (12 x 378)



# Hypothesis 1: all states have different levels of stochastic (Q) and deterministic (C) variability 
# AICc: -3232.724  
mod1_klamath = list()
mod1_klamath$A = "zero" #no trend because we z scored
mod1_klamath$Z = "identity"
mod1_klamath$R = "zero" #all the sensors are same, so observation error should be same
mod1_klamath$Q = "diagonal and unequal" 
mod1_klamath$B = "diagonal and unequal" #assuming no species interactions
mod1_klamath$U = "zero" #no trend because we z scored 
mod1_klamath$C = "unequal" #Can set C to unequal because it is going off the Z matrix where I have already indicated how to split up the sites.
mod1_klamath$c = zscore(matrixc_klamath)
mod1_klamath.fit = MARSS(transformed_dat_klamath[1:10,], model=mod1_klamath)
MARSSparamCIs(mod1_klamath.fit)
plot(mod1_klamath.fit)


# Hypothesis 2: all states have the same levels of stochastic (Q) and deterministic (C) variability mod11_klamath = list()
# AICc: -3176.907 
mod2_klamath = list()
mod2_klamath$A = "zero" #no trend because we z scored
mod2_klamath$Z = "identity"
mod2_klamath$R = "zero" #all the sensors are same, so observation error should be same
mod2_klamath$Q = "diagonal and equal" 
mod2_klamath$B = "diagonal and equal" #assuming no species interactions
mod2_klamath$U = "zero" #no trend because we z scored 
mod2_klamath$C = "equal" #Can set C to unequal because it is going off the Z matrix where I have already indicated how to split up the sites.
mod2_klamath$c = zscore(matrixc_klamath)
mod2_klamath.fit = MARSS(transformed_dat_klamath[1:5,], model=mod2_klamath)
plot(mod2_klamath.fit)

# Next steps:
# Habitat type: Creeks vs ponds vs Klamath
# Watershed: Creeks & ponds watershed 1 / Creeks & ponds watershed 2

## How do we modify natrices? 

# 1st: group time series into categories
hypothesis = c("Klamath","Klamath","Klamath","Klamath",
               "Klamath","Klamath","Klamath","Klamath",
               "Klamath","Klamath","Klamath","Klamath")
# 2nd: build C matrix (12 x 1)
mod2_klamath$C = matrix(hypothesis)
mod2_klamath$c = zscore(matrixc_klamath)
mod2_klamath.fit = MARSS(transformed_dat_klamath[1:5,], model=mod2_klamath)
plot(mod2_klamath.fit)

# 3rd: build Q matrix (12 x 12, with "C vector" in its diagonal)
Q<-diag(12)
diag(Q)=hypothesis
Q # Get rid of quotes in zeros

# 4th: B, identical as Q
B=Q


# Remaining questions
# 1) should we cross hypotheses or not? (unequal vs equal; ...)
# 2) Read more on z-scoring (for now we z-score both covariates and variates)