################################################################################
## Modeling survival and productivity as a function of holding time  
## for Texas translocated northern bobwhites
## Written by Amanda L. Schmidt
## 2023-03-30
################################################################################

################################################################################
## Survival Data - Holding time
load("nobo-holding-surv-data.gzip")
str(jags.data)

## Explanation of data in list jags.data

#ech - encounter history for each radio-tagged Texas translocated bobwhite for each day of the breeding season; 1 = alive, 0 = dead, NA = not yet part of the study, already recorded dead, or censored
#holdmat - holding times (in hours) standardized to mean 0 and standard deviation 1 for each radio-tagged Texas translocated bobwhite for each day of the breeding season
#first - first marking occasion
#last - last recapture/recovery occasion
#nind - number of individuals in the data set

#dimensions for ech and holdmat are equal to the number of individual bobwhites (nrow = 232) and the number of days during the breeding season (ncol = 187)
#length of first and last vectors is equal to the number of individual bobwhites (nrow = 232)

## Load library
library(jagsUI)

## Model
filename = "nobo-holding-surv-jag-model.JAG"
cat("model {
# Texas Translocation Holding Time Model

# Priors and constraints
    
 b0 ~ dnorm(0, 0.001) # intercept
 b1 ~ dnorm(0, 0.001) # effect of holding time
 
 # Likelihood 

 for (i in 1:nind){
   for (j in (first[i]+1):last[i]){ #known fate, detection = 1
     logit(phi[i,j]) <- b0 + b1*holdmat[i,j-1]
     
     mu[i,j] <- phi[i,j]*ech[i,j-1]
     ech[i,j] ~ dbern(mu[i,j]) #bernoulli
        
   
   } 
 } 

}", file = filename)

## Parameters to monitor
jags.par <- c("b0", "b1")

# Compile model in jags
m.surv <- jags(data = jags.data, parameters.to.save = jags.par, model.file = "nobo-holding-surv-jag-model.JAG", 
           n.chains = 3, n.thin = 1, n.iter = 10000, n.burnin = 5, n.adapt = 1000, DIC = FALSE)

m.surv

## Visually inspect convergence
plot(m.surv)

## Posterior estimation

#load in holding matrices
load("holdingmatrix_scale.gzip")
load("holdingmatrix.gzip")

#explanation of holding matrix data

#holdmat_scale - holding times (in hours) standardized to mean 0 and standard deviation 1 for each radio-tagged Texas translocated bobwhite for each day of the breeding season
#holdmat_new - holding times (in hours) for each radio-tagged Texas translocated bobwhite for each day of the breeding season

mean_hold_time <- mean(holdmat_new, na.rm = T) #107.07 hours
sd_hold_time <- sd(holdmat_new, na.rm = T) #91.02 hours

#extract
holding_vec <- seq(min(holdmat_scale, na.rm = T), 
                   max(holdmat_scale, na.rm = T), 0.1)

surv <- matrix(NA,nrow =length(m.surv$sims.list$b0), ncol = length(holding_vec))

for(i in 1:length(m.surv$sims.list$b0)){
  for(j in 1:length(holding_vec))
    surv[i,j] <- plogis(m.surv$sims.list$b0[i] + m.surv$sims.list$b1[i] * holding_vec[j])
}

#breeding season survival estimates for various holding times
mean_surv <- apply(surv, 2, mean)^180 #180 days breeding season 
upper_surv <- apply(surv, 2, function(x) quantile(x, probs = c(0.975)))^180
lower_surv <- apply(surv, 2, function(x) quantile(x, probs = c(0.025)))^180

#estimating specific values
# 48 hour holding time
hold.48 <- (48 - mean_hold_time)/sd_hold_time
h48 <- plogis(m.surv$sims.list$b0 + m.surv$sims.list$b1 * hold.48)

# 144 hour (6 days) holding time
hold.144 <- (144 - mean_hold_time)/sd_hold_time
h144 <- plogis(m.surv$sims.list$b0 + m.surv$sims.list$b1 * hold.144)

# 288 hour (12 days) holding time
hold.288 <- (288 - mean_hold_time)/sd_hold_time
h288 <- plogis(m.surv$sims.list$b0 + m.surv$sims.list$b1 * hold.288)

#breeding season survival estimates for three holding times
surv1 <- data.frame(h48, h144, h288)
mean_surv1 <- apply(surv1, 2, mean)^180 #180 days breeding season 
upper_surv1 <- apply(surv1, 2, function(x) quantile(x, probs = c(0.975)))^180
lower_surv1 <- apply(surv1, 2, function(x) quantile(x, probs = c(0.025)))^180

rm(list=ls())

###############################################################################
## Productivity Data - Holding time
load("nobo-holding-prod-data.gzip")
str(jags.data)

## Explanation of data in list jags.data

#totalnest - total number of nests each individual initiated
#clutchsize - number of eggs in each nest
#chicks - number of chicks hatched in each nest
#holdingtime - holding times (in hours) standardized to mean 0 and standard deviation 1 for each radio-tagged Texas translocated bobwhite
#sex - indicator that a nest was incubated (1) by a male or (0) female
#anynest - value of 1 = individual initiated at least one nest, value of 0 = individual did not initiate a nest
#nind - number of individuals in the data set

#dimensions for clutchsize and chicks is equal to the number of individual bobwhites (nrow = 232) and the maximum number of nests any bird initiated in a given year (ncol = 3)
##a value of NA indicates a) the data was not recorded for a given nest or b) an individual did not initiate that nest. The data anynest and totalnest are used in JAGS to index only over the nests an individual actually initiated

#length of totalnest, holdingtime, sex, and anynest vectors is equal to the number of individual bobwhites (nrow = 232)

## Load library
library(jagsUI)

## Parameters to monitor
jags.par <- c("n0", "n1", "n2", 
              "c0", "c1", "c2",
              "b0", "b1", "b2")

## Compile model in jags
m.prod <- jags(data = jags.data, parameters.to.save = jags.par, model.file = "nobo-holding-prod-jag-model.JAG", 
           n.chains=3, n.thin = 1, n.iter = 30000, n.burnin = 5, n.adapt = 1000, DIC = FALSE)

m.prod

## Visually inspect convergence
plot(m.prod)

## Posterior estimation

#load in holding matrices
load("holdingtime_scale.gzip")
load("holdingtime.gzip")

#explanation of holding time vector

#holdingtime_scale - holding times (in hours) standardized to mean 0 and standard deviation 1 for each radio-tagged Texas translocated bobwhite
#holdingtime - holding times (in hours) for each radio-tagged Texas translocated bobwhite

mean_hold_time <- mean(holdingtime, na.rm = T) #115.62 hours
sd_hold_time <- sd(holdingtime, na.rm = T) #95.92 hours

#extract nest prop
holding_vec <- seq(min(holdingtime_scale, na.rm = T), 
                   max(holdingtime_scale, na.rm = T), 0.1)

nestprop <- matrix(NA,nrow =length(m.prod$sims.list$n0), ncol = length(holding_vec))

for(i in 1:length(m.prod$sims.list$n0)){
  for(j in 1:length(holding_vec))
    nestprop[i,j] <- exp(m.prod$sims.list$n0[i] + m.prod$sims.list$n1[i] * holding_vec[j])
}

##nest propensity estimates for various holding times
mean_nestprop <- apply(nestprop, 2, mean) 
upper_nestprop <- apply(nestprop, 2, function(x) quantile(x, probs = c(0.975)))
lower_nestprop <- apply(nestprop, 2, function(x) quantile(x, probs = c(0.025)))

#estimating specific values
# 48 hour holding time
hold.48 <- (48 - mean_hold_time)/sd_hold_time
h48 <- exp(m.prod$sims.list$n0 + m.prod$sims.list$n1 * hold.48)

# 144 hour (6 days) holding time
hold.144 <- (144 - mean_hold_time)/sd_hold_time
h144 <- exp(m.prod$sims.list$n0 + m.prod$sims.list$n1 * hold.144)

# 288 hour (12 days) holding time
hold.288 <- (288 - mean_hold_time)/sd_hold_time
h288 <- exp(m.prod$sims.list$n0 + m.prod$sims.list$n1 * hold.288)

#nest propensity estimates for three holding times
nestprop1 <- data.frame(h48, h144, h288)
mean_nest1 <- apply(nestprop1, 2, mean)
upper_nest1 <- apply(nestprop1, 2, function(x) quantile(x, probs = c(0.975)))
lower_nest1 <- apply(nestprop1, 2, function(x) quantile(x, probs = c(0.025)))

#extract clutch size
clutch <- matrix(NA,nrow =length(m.prod$sims.list$c0), ncol = length(holding_vec))

for(i in 1:length(m.prod$sims.list$c0)){
  for(j in 1:length(holding_vec))
    clutch[i,j] <- exp(m.prod$sims.list$c0[i] + m.prod$sims.list$c1[i] * holding_vec[j])
}

##clutch size estimates for various holding times
mean_clutch <- apply(clutch, 2, mean) 
upper_clutch <- apply(clutch, 2, function(x) quantile(x, probs = c(0.975)))
lower_clutch <- apply(clutch, 2, function(x) quantile(x, probs = c(0.025)))

#estimating specific values
# 48 hour holding time
hold.48 <- (48 - mean_hold_time)/sd_hold_time
h48 <- exp(m.prod$sims.list$c0 + m.prod$sims.list$c1 * hold.48)

# 144 hour (6 days) holding time
hold.144 <- (144 - mean_hold_time)/sd_hold_time
h144 <- exp(m.prod$sims.list$c0 + m.prod$sims.list$c1 * hold.144)

# 288 hour (12 days) holding time
hold.288 <- (288 - mean_hold_time)/sd_hold_time
h288 <- exp(m.prod$sims.list$c0 + m.prod$sims.list$c1 * hold.288)

#clutch estimates for three holding times
clutch1 <- data.frame(h48, h144, h288)
mean_clutch1 <- apply(clutch1, 2, mean)
upper_clutch1 <- apply(clutch1, 2, function(x) quantile(x, probs = c(0.975)))
lower_clutch1 <- apply(clutch1, 2, function(x) quantile(x, probs = c(0.025)))

#extract nest success
hatch <- matrix(NA,nrow =length(m.prod$sims.list$b0), ncol = length(holding_vec))

for(i in 1:length(m.prod$sims.list$b0)){
  for(j in 1:length(holding_vec))
    hatch[i,j] <- exp(m.prod$sims.list$b0[i] + m.prod$sims.list$b1[i] * holding_vec[j])
}

##nest success estimates for various holding times
mean_hatch <- apply(hatch, 2, mean) 
upper_hatch <- apply(hatch, 2, function(x) quantile(x, probs = c(0.975)))
lower_hatch <- apply(hatch, 2, function(x) quantile(x, probs = c(0.025)))

#estimating specific values
# 48 hour holding time
hold.48 <- (48 - mean_hold_time)/sd_hold_time
h48 <- exp(m.prod$sims.list$b0 + m.prod$sims.list$b1 * hold.48)

# 144 hour (6 days) holding time
hold.144 <- (144 - mean_hold_time)/sd_hold_time
h144 <- exp(m.prod$sims.list$b0 + m.prod$sims.list$b1 * hold.144)

# 288 hour (12 days) holding time
hold.288 <- (288 - mean_hold_time)/sd_hold_time
h288 <- exp(m.prod$sims.list$b0 + m.prod$sims.list$b1 * hold.288)

#nest success for three holding times
hatch1 <- data.frame(h48, h144, h288)
mean_hatch1 <- apply(hatch1, 2, mean)
upper_hatch1 <- apply(hatch1, 2, function(x) quantile(x, probs = c(0.975)))
lower_hatch1 <- apply(hatch1, 2, function(x) quantile(x, probs = c(0.025)))

#extract fecundity
fecundity <- matrix(NA,nrow =length(m.prod$sims.list$b0), ncol = length(holding_vec))

for(i in 1:length(m.prod$sims.list$b0)){
  for(j in 1:length(holding_vec))
    fecundity[i,j] <- (exp(m.prod$sims.list$b0[i] + m.prod$sims.list$b1[i] * holding_vec[j])) * 
                      (exp(m.prod$sims.list$c0[i] + m.prod$sims.list$c1[i] * holding_vec[j])) * 
                      (exp(m.prod$sims.list$n0[i] + m.prod$sims.list$n1[i] * holding_vec[j])) * (1/2)
}

##fecundity estimates for various holding times
mean_fecundity <- apply(fecundity, 2, mean) 
upper_fecundity <- apply(fecundity, 2, function(x) quantile(x, probs = c(0.975)))
lower_fecundity <- apply(fecundity, 2, function(x) quantile(x, probs = c(0.025)))

#estimating specific values
# 48 hour holding time
hold.48 <- (48 - mean_hold_time)/sd_hold_time
h48 <- (exp(m.prod$sims.list$b0 + m.prod$sims.list$b1 * hold.48)) * 
       (exp(m.prod$sims.list$c0 + m.prod$sims.list$c1 * hold.48)) * 
       (exp(m.prod$sims.list$n0 + m.prod$sims.list$n1 * hold.48)) * (1/2)

# 144 hour (6 days) holding time
hold.144 <- (144 - mean_hold_time)/sd_hold_time
h144 <- (exp(m.prod$sims.list$b0 + m.prod$sims.list$b1 * hold.144)) * 
        (exp(m.prod$sims.list$c0 + m.prod$sims.list$c1 * hold.144)) * 
        (exp(m.prod$sims.list$n0 + m.prod$sims.list$n1 * hold.144)) * (1/2)

# 288 hour (12 days) holding time
hold.288 <- (288 - mean_hold_time)/sd_hold_time
h288 <- (exp(m.prod$sims.list$b0 + m.prod$sims.list$b1 * hold.288)) * 
        (exp(m.prod$sims.list$c0 + m.prod$sims.list$c1 * hold.288)) * 
        (exp(m.prod$sims.list$n0 + m.prod$sims.list$n1 * hold.288)) * (1/2) 

#fecundity estimates for three holding times
fecundity1 <- data.frame(h48, h144, h288)
mean_fecundity1 <- apply(fecundity1, 2, mean)
upper_fecundity1 <- apply(fecundity1, 2, function(x) quantile(x, probs = c(0.975)))
lower_fecundity1 <- apply(fecundity1, 2, function(x) quantile(x, probs = c(0.025)))

rm(list=ls())

###############################################################################
