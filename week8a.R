### ESPM 174A, Fall 2021
## Week 8 - Oct 11, 2021
library(MARSS)
library(MASS)
library(tidyverse)

## 8.1: MARSS with an environmental covariate (several processes)
# Create some sample data
nspp = 5 # 5 populations
years = 20 # monitored over 20 years
reps = 3 # we are measuring each population 3 times

# Z will tell MARSS which Y (observation) goes to which X (true population size)
Z = factor(rep(1:nspp, each=reps))
Z = matrix(0,nspp*reps,nspp)
for(i in 1:nspp) Z[(1:reps)+reps*(i-1),i]=1
Z # See species are grouped in 'chunks' of three (3 spatial replicates x species)

# We will assume independent and identically distributed (iid) observation error across species and replicates 
r = 0.01 # this is the process error variance
R=diag(r,nspp*reps)
R # This is the process error variance-covariance matrix

# Now we will create an environmental driver (covariate)
set.seed(1) # so that our values will not change
driver=cumsum(rnorm(years,0,.05)) # this creates a covariate w/ mean = 0 and sigma = 0.05
driver # This is our covariate data
plot.ts(driver)

# We will assume that the driver affects the different species in potentially different ways
set.seed(1) # so that our values will not change
C = rnorm(nspp,0,1)
C # These are the covariate effects (what we will normally want to estimate)

# We will add, on top of the driver, some independent and identically distributed (iid) process error
set.seed(1) # so that our values will not change
q = 0.01 # this is the process error variance
Q=diag(q,nspp)
Q # This is the process error variance-covariance matrix

# Now let's deal with the B matrix: we will assume no spp interactions for now
# That is, no interespecific interactions (off-diagonal) or density-dependence (diagonal), so random walk
B = diag(1,nspp)
B # aka 'identity' B matrix

# No drift in the random walk
U = matrix(0,nspp,1)

# Let's generate the data
N = matrix(0,nspp,years) # Here we will store the true population sizes
O = matrix(0,nspp*reps,years) # Here we will store the observations

# Set the initial conditions using lognormal; N is log(counts)
set.seed(1) # so that our values will not change
N[,1] = log(100) + rnorm(nspp,0,1)
O[,1] = Z%*%N[,1] + mvrnorm(1,matrix(0,nspp*reps,1),R)
for(t in 2:years){
  proc.err.t = C*driver[t]+mvrnorm(1,U,Q) # This is the process error at time t. It is the result of covariate effects + drift (here set to zero) and the iid process error
  N[,t]=B%*%N[,t-1,drop=FALSE] + proc.err.t
  O[,t] = Z%*%N[,t] + mvrnorm(1,matrix(0,nspp*reps,1),R)
}
matplot(t(N),type="l") # This is how the simulated data (N) looks like
# End of the data simulation code


## 8.2: Now fit a MARSS model to the generated data
dat=O # Letter "O", for Observations
mod0 = list() # This will be the list of parameters, each is a matrix

# X equation (process model)
mod0$B = "identity" # no species interactions allowed (not even density dependence)
mod0$U = "zero" # no drift
mod0$Q = "diagonal and equal" # single process error across pops and replicates
mod0$C = "unequal" # each population may be affected by the covariate differently
mod0$c = matrix(driver,1,years) # Covariate data, what we cooked earlier

# Y equation (observation model)
mod0$Z = factor(rep(1:nspp, each=reps))
Z # This Z matrix preserves the order of the data (i.e., replicates within pops)
mod0$R = "diagonal and equal" # Single observation error across pops and replicates
mod0$A = "zero"

# Fit the MARSS model
mod0.fit = MARSS(dat, model=mod0)

# And see whether we are recovering the covariate effects
plot(coef(mod0.fit)$C,C,xlab="estimated driver effect",ylab="true driver effect")
abline(0,1) # pretty good!


## 8.3: Deal with collinearity, fit a model affected by collinear covariates
# Generate collinear covariates
set.seed(1)
mod0$c = matrix(driver,1,years) # This was our covariate data
newc = mod0$c+rnorm(20,0,0.01) # New covariate data (same, with a little noise on top)
plot(mod0$c,newc) # highly correlated
mod0$c = rbind(mod0$c,newc) # New set of covariates
cor(mod0$c[1,],mod0$c[2,]) # highly correlated

# Fit MARSS model
mod0.fit = MARSS(dat, model=mod0) 
MARSSparamCIs(mod0.fit) # coefficients are inflated in opposite directions

# Now let's orthogonalize
z1 <- mod0$c[1,]
fit.z2 <- lm(mod0$c[2,] ~ mod0$c[1,])
resids_z2 <- residuals(fit.z2, type="response")
Z <- cbind(z1,resids_z2)
# standardize the variance
plot(z1,resids_z2)
Z<-t(zscore(t(Z)))
cor(Z[,1],Z[,2]) # drivers are not correlated any more

# Fit new model
mod1 = mod0
mod1$c <- t(Z)
mod1.fit = MARSS(dat, model=mod1)
MARSSparamCIs(mod1.fit) # much better


## 8.4: Dealing with seasonality via fixed factors
# Load Lake Washington Plankton data set
?lakeWAplanktonTrans # more info on the dataset
fulldat = lakeWAplanktonTrans
years = fulldat[,"Year"]>=1965 & fulldat[,"Year"]<1975
phyto = c("Diatoms", "Greens", "Bluegreens",
           "Unicells", "Other.algae")
dat = t(fulldat[years, phyto])
the.mean = apply(dat,1,mean,na.rm=TRUE)
the.sigma = sqrt(apply(dat,1,var,na.rm=TRUE))
dat = (dat-the.mean)*(1/the.sigma)
par(mfrow=c(2,1))
matplot(t(dat), type="l")
# Z-score covariates
covariates = rbind(
  Temp = fulldat[years,"Temp"],
  TP = fulldat[years,"TP"])
covariates = zscore(covariates)
matplot(t(covariates), type="l")

# Fixed seasonal factors
TT = ncol(dat) # number of time periods/samples
period = 12 # number of "seasons" (e.g., 12 months per year)
per.1st = 1 # first "season" (e.g., Jan = 1, July = 7)
c = diag(period) # create factors for seasons
for(i in 2:(ceiling(TT/period))) {c = cbind(c,diag(period))}
# give row names
rownames(c) = month = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
print(c)
dim(c)

# Specify model, fit model
C = "unconstrained"
B = "diagonal and unequal" # Each taxon has unique density-dependence
Q = "diagonal and unequal" # Independent process errors
U = "zero" # we demeaned data
Z = "identity"  # Each obs time series is associated with only one process
A = "zero" # The data are demeaned & fluctuate around a mean
R = "diagonal and equal" # Independent obs errors, similar variance
model.list = list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=c) # Fit model
seas.mod.1 = MARSS(dat, model=model.list)

# Estimated seasonal effects
seas.1 = coef(seas.mod.1,type="matrix")$C
rownames(seas.1) = phyto
colnames(seas.1) = month
matplot(t(seas.1), type="l") # OK


## 8.5: Dealing with seasonality via a Fourier series
# Create Fourier series
cos.t = cos(2 * pi * seq(TT) / period)
sin.t = sin(2 * pi * seq(TT) / period)
c.Four = rbind(cos.t,sin.t)
cor(c.Four[1,],c.Four[2,]) # not correlated!
matplot(t(c.Four), type="l")

# Specify and fit MARSS model
model.list = list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=c.Four)
seas.mod.2 = MARSS(dat, model=model.list)

# Recover seasonal effects
C.2 = coef(seas.mod.2, type="matrix")$C
# The time series of net seasonal effects
seas.2 = C.2 %*% c.Four[,1:period]
rownames(seas.2) = phyto
colnames(seas.2) = month
matplot(t(seas.2), type="l") # OK

# Plot seasonal effects
# From fixed effects model:
par(mfrow=c(2,1), mar=c(2,4,2,2)) 
matplot(t(seas.1),type="l",bty="n",xaxt="n", ylab="Fixed monthly", col=1:5)
axis(1,labels=month, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phyto, cex=0.6, col=1:5)

# From Fourier series  model:
matplot(t(seas.2),type="l",bty="n",xaxt="n",ylab="Fourier", col=1:5)
axis(1,labels=month, at=1:12,las=2,cex.axis=0.75)
legend("topright", lty=1:5, legend=phyto, cex=0.6, col=1:5)

# Compare models using AIC
data.frame(Model=c("Fixed", "Fourier"),
           AICc=round(c(seas.mod.1$AICc,
                        seas.mod.2$AICc),1))


## 8.6: Now let's fit model with seasonality AND an additional covariate (Total Phosphorus, TP)
newcovarsFour_TP<-rbind(c.Four, "TP"=covariates[2,])
head(newcovarsFour_TP)
matplot(t(newcovarsFour_TP), type="l", col=c("black","red","blue"))
newcovarsFour_TP[,1:6]
# Fit new MARSS model
C = "unconstrained" # See it was already unconstrained, so there was no need to specify C again (it would adapt automatically to the new dimensions of the covariate data)
model.list = list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=newcovarsFour_TP)
seas.mod.3 = MARSS(dat, model=model.list)
MARSSparamCIs(seas.mod.3)


# 8.7: Model diagnostics
par(mfrow=c(5,2), mai=c(0.1,0.5,0.2,0.1), omi=c(0.5,0,0,0))
  for (j in 1:5) {
    plot.ts(residuals<-MARSSresiduals(seas.mod.3, type = "tt1")$model.residuals[j, ],
            ylab = "Residual")
    abline(h = 0, lty = "dashed")
    acf(residuals,na.action = na.pass)
  }


## 8.8: Model diagnostics before including the Fourier series
newcovarsonlyTP<-newcovarsFour_TP[3,, drop=FALSE]
newcovarsonlyTP # OK
model.list = list(B=B,U=U,Q=Q,Z=Z,A=A,R=R,C=C,c=newcovarsonlyTP)
seas.mod.4 = MARSS(dat, model=model.list)
# Plot diagnostics again
par(mfrow=c(5,2), mai=c(0.1,0.5,0.2,0.1), omi=c(0.5,0,0,0))
for (j in 1:5) {
  plot.ts(residuals<-MARSSresiduals(seas.mod.4, type = "tt1")$model.residuals[j, ],
          ylab = "Residual")
  abline(h = 0, lty = "dashed")
  acf(residuals,na.action = na.pass)
} # ACF is much worse, and model did not even converge

### End