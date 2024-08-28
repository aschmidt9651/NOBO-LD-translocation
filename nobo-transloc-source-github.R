################################################################################
## Modeling survival, nest survival, and productivity as a function of source 
## population for translocated northern bobwhites
## Written by Amanda L. Schmidt
## 2023-03-30
################################################################################

################################################################################
## Survival Data - Source population
load("nobo-source-surv-data.gzip")
str(jags.data)

## Explanation of data in list jags.data

#ech - encounter history for each radio-tagged bobwhite for each week of the breeding season; 1 = alive, 0 = dead, NA = not yet part of the study, already recorded dead, or censored
#FLtrans - indicator that an individual is a Florida translocated bobwhite (1) or not (0)
#TXtrans - indicator that an individual is a Texas translocated bobwhite (1) or not (0)
#sex - indicator that an individual is male (1) or not (0)
#first - first marking occasion
#last - last recapture/recovery occasion
#nind - number of individuals in the data set

#dimensions for ech is equal to the number of individual bobwhites (nrow = 323) and the number of weeks during the breeding season (ncol = 27)
#length of FLtrans, TX trans, sex, first, and last vectors is equal to the number of individual bobwhites (nrow = 323)

## Load library
library(jagsUI)

## Model
filename = "nobo-source-surv-jag-model.JAG"
cat("model {
# Source Population Model

# Priors and constraints
    
 b0 ~ dnorm(0, 0.001) # intercept (FL resident survival)
 b1 ~ dnorm(0, 0.001) # effect of FL translocated bird on survival
 b2 ~ dnorm(0, 0.001) # effect of TX translocated bird on survival
 b3 ~ dnorm(0, 0.001) # effect sex (male) on survival
    
 # Likelihood 

 for (i in 1:nind){
   for (j in (first[i]+1):last[i]){ #known fate, detection = 1
     logit(phi[i,j]) <- b0 + b1*FLtrans[i] + b2*TXtrans[i] + b3*sex[i]
     
     mu[i,j] <- phi[i,j]*ech[i,j-1]
     ech[i,j] ~ dbern(mu[i,j]) #bernoulli
        
   
   } 
 } 

}", file = filename)

## Parameters to monitor
jags.par <- c("b0", "b1", "b2", "b3")

# Compile model in jags
m.surv <- jags(data = jags.data, parameters.to.save = jags.par, model.file = "nobo-source-surv-jag-model.JAG", 
             n.chains = 3, n.thin = 1, n.iter = 10000, n.burnin = 5, n.adapt = 1000, DIC = FALSE)

m.surv

## Visually inspect convergence
plot(m.surv)

## Posterior estimation

#extract
FLresident <- FLtrans <- TXtrans <- male <- numeric()

for(i in 1:length(m.surv$sims.list$b0)){
    FLresident[i] <- plogis(m.surv$sims.list$b0[i])
    FLtrans[i] <- plogis(m.surv$sims.list$b0[i] + m.surv$sims.list$b1[i])
    TXtrans[i] <- plogis(m.surv$sims.list$b0[i] + m.surv$sims.list$b2[i])
    male[i] <- plogis(m.surv$sims.list$b0[i] + m.surv$sims.list$b3[i])
}

surv <- data.frame(FLresident,FLtrans,TXtrans,male)

#breeding season survival estimates
mean_surv <- apply(surv, 2, mean)^26 #26 week breeding season
upper_surv <- apply(surv, 2, function(x) quantile(x, probs = c(0.975)))^26
lower_surv <- apply(surv, 2, function(x) quantile(x, probs = c(0.025)))^26

#difference in breeding season survival among source populations - are they significant?
#i.e., comparing proportion of the posterior distribution with the same sign as the mean

#FLresident v TXtrans
sum(surv$FLresident > surv$TXtrans)/length(surv$FLresident)

#FLresident v FLtrans
sum(surv$FLresident < surv$FLtrans)/length(surv$FLresident)

#FLtrans v TXtrans
sum(surv$FLtrans > surv$TXtrans)/length(surv$FLtrans)

rm(list=ls())

###############################################################################
## Nest Survival Data - Source population
load("nobo-source-nestsurv-data.gzip")
str(jags.data)

## Explanation of data in list jags.data

#ech - encounter history for each nest for each day of the breeding season; 1 = active or hatched (last encounter), 0 = failed, NA = not yet observed, already failed or hatched
#FLtrans - indicator that a nest is from a Florida translocated bobwhite (1) or not (0)
#TXtrans - indicator that a nest is from a Texas translocated bobwhite (1) or not (0)
#sex - indicator that a nest was incubated (1) by a male or (0) female
#first - occasion nest was first identified
#last -  occasion nest hatched or failed
#nnest - number of nests in the data set

#dimensions for ech is equal to the number of bobwhite nests (nrow = 188) and the number of days any nest was active during the breeding season (ncol = 200)
#length of FLtrans, TX trans, sex, first, and last vectors is equal to the number of bobwhite nests (nrow = 188)

## Load library
library(jagsUI)

## Model
filename = "nobo-source-nestsurv-jag-model.JAG"
cat("model {
# Source Population Model

# Priors and constraints
    
 b0 ~ dnorm(0, 0.001) # intercept (FLresident nest survival)
 b1 ~ dnorm(0, 0.001) # effect of FLtrans bird on nest survival
 b2 ~ dnorm(0, 0.001) # effect of TXtrans bird on nest survival
 b3 ~ dnorm(0, 0.001) # effect of sex (male) on nest survival
  
 for(i in 1:nnest){ # loop through nests
 
 # Likelihood

   for (j in (first[i]+1):last[i]){ #known fate, detection = 1
     logit(phi[i,j]) <- b0 + b1*FLtrans[i] + b2*TXtrans[i] + b3*sex[i]
     
     mu[i,j] <- phi[i,j]*ech[i,j-1]
     ech[i,j] ~ dbern(mu[i,j]) #bernoulli
        
   
   } 
 } 

}", file = filename)


## Parameters to monitor
jags.par <- c("b0", "b1", "b2", "b3")

# Compile model in jags 
m.nestsurv <- jags(data = jags.data, parameters.to.save = jags.par, model.file = "nobo-source-nestsurv-jag-model.JAG", 
           n.chains = 3, n.thin = 1, n.iter = 10000, n.burnin = 5, n.adapt = 1000, DIC = FALSE)

m.nestsurv

## Visually inspect convergence
plot(m.nestsurv)

## Posterior estimation

#extract
FLresident <- FLtrans <- TXtrans <- male <- numeric()

for(i in 1:length(m.nestsurv$sims.list$b0)){
  FLresident[i] <- plogis(m.nestsurv$sims.list$b0[i])
  FLtrans[i] <- plogis(m.nestsurv$sims.list$b0[i] + m.nestsurv$sims.list$b1[i])
  TXtrans[i] <- plogis(m.nestsurv$sims.list$b0[i] + m.nestsurv$sims.list$b2[i])
  male[i] <- plogis(m.nestsurv$sims.list$b0[i] + m.nestsurv$sims.list$b3[i])
}

nestsurv <- data.frame(FLresident,FLtrans,TXtrans,male)

#nest survival estimates
mean_nestsurv <- apply(nestsurv, 2, mean)^23 #23 day incubation 
upper_nestsurv <- apply(nestsurv, 2, function(x) quantile(x, probs = c(0.975)))^23
lower_nestsurv <- apply(nestsurv, 2, function(x) quantile(x, probs = c(0.025)))^23

#difference in nest survival among source populations - are they significant?
#i.e., comparing proportion of the posterior distribution with the same sign as the mean

#FLresident v TXtrans
sum(nestsurv$FLresident > nestsurv$TXtrans)/length(nestsurv$FLresident)

#FLresident v FLtrans
sum(nestsurv$FLresident < nestsurv$FLtrans)/length(nestsurv$FLresident)

#FLtrans v TXtrans
sum(nestsurv$FLtrans > nestsurv$TXtrans)/length(nestsurv$FLtrans)

rm(list=ls())

###############################################################################
## Productivity Data - Source population
load("nobo-source-prod-data.gzip")
str(jags.data)

## Explanation of data in list jags.data

#totalnest - total number of nests each individual initiated
#clutchsize - number of eggs in each nest
#chicks - number of chicks hatched in each nest
#FLtrans - indicator that a nest is from a Florida translocated bobwhite (1) or not (0)
#TXtrans - indicator that a nest is from a Texas translocated bobwhite (1) or not (0)
#sex - indicator that a nest was incubated (1) by a male or (0) female
#anynest - value of 1 = individual initiated at least one nest, value of 0 = individual did not initiate a nest
#nind - number of individuals in the data set

#dimensions for clutchsize and chicks is equal to the number of individual bobwhites (nrow = 323) and the maximum number of nests any bird initiated in a given year (ncol = 4)
##a value of NA indicates a) the data was not recorded for a given nest or b) an individual did not initiate that nest. The data anynest and totalnest are used in JAGS to index only over the nests an individual actually initiated

#length of totalnest, FLtrans, TX trans, sex, and anynest vectors is equal to the number of individual bobwhites (nrow = 323)

## Load library
library(jagsUI)

## Parameters to monitor
jags.par <-c("n0", "n1", "n2", "n3", 
            "c0", "c1", "c2", "c3", 
             "b0", "b1", "b2", "b3")

## Compile model in jags
m.prod <- jags(data = jags.data, parameters.to.save = jags.par, model.file = "nobo-source-prod-jag-model.JAG", 
           n.chains=3, n.thin = 1, n.iter = 30000, n.burnin = 5, n.adapt = 1000, DIC = FALSE)

m.prod

## Visually inspect convergence
plot(m.prod)

## Posterior estimation

#extract nest propensity
FLresident <- FLtrans <- TXtrans <- numeric()

for(i in 1:length(m.prod$sims.list$n0)){
  FLresident[i] <- exp(m.prod$sims.list$n0[i])
  FLtrans[i] <- exp(m.prod$sims.list$n0[i] + m.prod$sims.list$n1[i])
  TXtrans[i] <- exp(m.prod$sims.list$n0[i] + m.prod$sims.list$n2[i])
 }

nestprop <- data.frame(FLresident,FLtrans,TXtrans)

#nest propensity estimates
mean_nestprop <- apply(nestprop, 2, mean)
upper_nestprop <- apply(nestprop, 2, function(x) quantile(x, probs = c(0.975)))
lower_nestprop <- apply(nestprop, 2, function(x) quantile(x, probs = c(0.025)))

#difference in nest propensity among source populations - are they significant?
#i.e., comparing proportion of the posterior distribution with the same sign as the mean

#FLresident v TXtrans
sum(nestprop$FLresident > nestprop$TXtrans)/length(nestprop$FLresident)

#FLresident v FLtrans
sum(nestprop$FLresident < nestprop$FLtrans)/length(nestprop$FLresident)

#FLtrans v TXtrans
sum(nestprop$FLtrans > nestprop$TXtrans)/length(nestprop$FLtrans)

#extract clutch size
FLresident <- FLtrans <- TXtrans <- numeric()

for(i in 1:length(m.prod$sims.list$c0)){
  FLresident[i] <- exp(m.prod$sims.list$c0[i])
  FLtrans[i] <- exp(m.prod$sims.list$c0[i] + m.prod$sims.list$c1[i])
  TXtrans[i] <- exp(m.prod$sims.list$c0[i] + m.prod$sims.list$c2[i])
  }

clutch <- data.frame(FLresident,FLtrans,TXtrans)

#clutch size estimates
mean_clutch <- apply(clutch, 2, mean)
upper_clutch <- apply(clutch, 2, function(x) quantile(x, probs = c(0.975)))
lower_clutch <- apply(clutch, 2, function(x) quantile(x, probs = c(0.025)))

#difference in clutch size among source populations - are they significant?
#i.e., comparing proportion of the posterior distribution with the same sign as the mean

#FLresident v TXtrans
sum(clutch$FLresident > clutch$TXtrans)/length(clutch$FLresident)

#FLresident v FLtrans
sum(clutch$FLresident < clutch$FLtrans)/length(clutch$FLresident)

#FLtrans v TXtrans
sum(clutch$FLtrans > clutch$TXtrans)/length(clutch$FLtrans)

##extract nest success
FLresident <- FLtrans <- TXtrans <- numeric()

for(i in 1:length(m.prod$sims.list$b0)){
  FLresident[i] <- exp(m.prod$sims.list$b0[i])
  FLtrans[i] <- exp(m.prod$sims.list$b0[i] + m.prod$sims.list$b1[i])
  TXtrans[i] <- exp(m.prod$sims.list$b0[i] + m.prod$sims.list$b2[i])
}

hatch <- data.frame(FLresident,FLtrans,TXtrans)

#nest success estimates
mean_hatch <- apply(hatch, 2, mean)
upper_hatch <- apply(hatch, 2, function(x) quantile(x, probs = c(0.975)))
lower_hatch <- apply(hatch, 2, function(x) quantile(x, probs = c(0.025)))

#difference in nest success among source populations - are they significant?
#i.e., comparing proportion of the posterior distribution with the same sign as the mean

#FLresident v TXtrans
sum(hatch$FLresident > hatch$TXtrans)/length(hatch$FLresident)

#FLresident v FLtrans
sum(hatch$FLresident < hatch$FLtrans)/length(hatch$FLresident)

#FLtrans v TXtrans
sum(hatch$FLtrans > hatch$TXtrans)/length(hatch$FLtrans)

##extract fecundity
FLresident <- FLtrans <- TXtrans <- numeric()

for(i in 1:length(m.prod$sims.list$b0)){
  FLresident[i] <- exp(m.prod$sims.list$b0[i]) * exp(m.prod$sims.list$c0[i]) * exp(m.prod$sims.list$n0[i]) * (1/2)
  FLtrans[i] <- exp(m.prod$sims.list$b0[i] + m.prod$sims.list$b1[i]) * exp(m.prod$sims.list$c0[i] + m.prod$sims.list$c1[i]) * exp(m.prod$sims.list$n0[i] + m.prod$sims.list$n1[i]) * (1/2)
  TXtrans[i] <- exp(m.prod$sims.list$b0[i] + m.prod$sims.list$b2[i]) * exp(m.prod$sims.list$c0[i] + m.prod$sims.list$c2[i]) * exp(m.prod$sims.list$n0[i] + m.prod$sims.list$n2[i]) * (1/2)
}

fecundity <- data.frame(FLresident,FLtrans,TXtrans)

#nest success estimates
mean_fecundity <- apply(fecundity, 2, mean)
upper_fecundity <- apply(fecundity, 2, function(x) quantile(x, probs = c(0.975)))
lower_fecundity <- apply(fecundity, 2, function(x) quantile(x, probs = c(0.025)))

#difference in fecundity among source populations - are they significant?
#i.e., comparing proportion of the posterior distribution with the same sign as the mean

#FLresident v TXtrans
sum(fecundity$FLresident > fecundity$TXtrans)/length(fecundity$FLresident)

#FLresident v FLtrans
sum(fecundity$FLresident < fecundity$FLtrans)/length(fecundity$FLresident)

#FLtrans v TXtrans
sum(fecundity$FLtrans > fecundity$TXtrans)/length(fecundity$FLtrans)

rm(list=ls())

###############################################################################
