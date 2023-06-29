##Libraries##

library(tidyverse)
library(nimble)
library(ggmcmc)
library(coda)
library(MCMCvis)
library(viridis)


##For fine-grained model

##Prep data for model##

##Bring in data files (FF.csv, Obs1.csv, Obs2.csv, seat.csv)
##Name the objects as FF, Obs1, Obs2, and seat

OBS <- array(c(Obs1, Obs2), dim = c(nrow(Obs1), ncol(Obs1), 2))

Obs1.total <- apply(Obs1, 1, sum)

Obs2.total <- apply(Obs2, 1, sum)

OBS.total <- cbind(Obs1.total, Obs2.total)

FF.total <- apply(FF, 1, sum)

#Set constants#

nspecies <- 31
nsites <- 321
nobs <- 2

##Nimble code##

code <- nimbleCode({
  
  #-Priors-#
  
  for(o in 1:nobs){
    
    #Observer-specific composition
    phi[1:(nspecies+1), o] ~ ddirch(phi.ones[1:(nspecies+1),o])
    
    for(i in 1:(nspecies+1)){
      
      phi.ones[i,o] <- 1
      
    }#end i
    
    #Intercept on observation error
    int.epsilon[o] ~ dnorm(0, 0.01)
    
  }#end o
  
  for(i in 1:nspecies){
    
    #Intercept for expected species-specific abundances
    lambda0[i] ~ dnorm(0, 0.01)
    
  }#end i
  
  #Effect of rear seat
  beta ~ dnorm(0, 0.01)
  
  #-Likelihood-#
  
  for(j in 1:nsites){
    
    #Front facing camera (true) composition
    FF[j,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[j])
    
    #Front facing camera (true) community-abundance
    FF.total[j] ~ dpois(lambda.total)
    
    for(o in 1:nobs){
      
      #Observer composition
      OBS[j,1:(nspecies+1),o] ~ dmulti(phi[1:(nspecies+1),o], OBS.total[j,o])
      
      #Observer community-abundance
      OBS.total[j,o] ~ dpois(lambda.total * E.epsilon[j,o])
      
      #Linear predictor of observation error
      log(E.epsilon[j,o]) <- int.epsilon[o] + beta * seat[j,o]
      
    }#end o
  }#end j
  
  for(i in 1:nspecies){
    
    #Composition as proportion of species-specific abundance to community-abundance
    pi[i] <- lambda[i]/lambda.total
    
    #Linear predictor of species-specific abundance
    log(lambda[i]) <- lambda0[i]
    
    for(o in 1:nobs){
      
      #Correction factor for species- and observer-specific observations
      correction[i,o] <- exp(int.epsilon[o]) * phi[i,o]/pi[i]
      
    }#end o
  }#end i
  
  #Constrain expected community-abundance to be sum of expected species-specific abundances
  lambda.total <- sum(lambda[1:nspecies])
  
})

##Compile data##

data <- list(FF = FF[,1:nspecies],
             FF.total = FF.total,
             OBS = OBS,
             OBS.total = OBS.total,
             seat = seat
)


constants <- list(nspecies = nspecies, nsites = nsites, nobs = 2)

##Initial values##

inits <- function(){list(pi = apply(FF[,1:nspecies]/FF.total, 2, mean, na.rm = TRUE),
                         int.epsilon = sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))),
                         beta = runif(1, -1, 1),
                         phi = apply(apply(OBS, c(1,2,3), max)/apply(apply(OBS, c(1,2,3), max), 1, sum), c(2,3), mean,
                                     na.rm = TRUE)
)}

##Parameters to save##

params <- c(
  "pi",
  "phi",
  "beta",
  "int.epsilon",
  "lambda0"
)

params2 <- c(
  "lambda",
  "lambda.total",
  "correction"
)

##MCMC settings##

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params, monitors2 = params2)

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 20000
nb <- 10000
nt <- 1
nt2 <- 10

##Run model##

f.out <- runMCMC(compiled.model$MCMC,
                 niter = ni, nburnin = nb,
                 nchains = nc, thin = nt, thin2 = nt2,
                 samplesAsCodaMCMC = TRUE)





######
##For coarse-grained model

##Bring in data files (FF_c.csv, Obs1_c.csv, Obs2_c.csv, seat.csv)
##Name the objects as FF, Obs1, Obs2, and seat


##Prep data for model##

OBS <- array(c(Obs1, Obs2), dim = c(nrow(Obs1), ncol(Obs1), 2))

Obs1.total <- apply(Obs1, 1, sum)

Obs2.total <- apply(Obs2, 1, sum)

OBS.total <- cbind(Obs1.total, Obs2.total)

FF.total <- apply(FF, 1, sum)

ID.obs <- ID.obs[order(ID.obs$group),]

##Set constants##

nspecies <- 15
nsites <- 321
nobs <- 2


#Run model code above#



##Compile data##

data <- list(FF = FF[,1:nspecies],
             FF.total = FF.total,
             OBS = OBS,
             OBS.total = OBS.total,
             seat = seat
)


constants <- list(nspecies = nspecies, nsites = nsites, nobs = 2)


##Initial values##

inits <- function(){list(pi = apply(FF[,1:nspecies]/FF.total, 2, mean, na.rm = TRUE),
                         int.epsilon = sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))),
                         beta = runif(1, -1, 1),
                         phi = apply(apply(OBS, c(1,2,3), max)/apply(apply(OBS, c(1,2,3), max), 1, sum), c(2,3), mean,
                                     na.rm = TRUE)
)}


##MCMC settings##

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params, monitors2 = params2)

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

##Run model##

c.out <- runMCMC(compiled.model$MCMC,
                 niter = ni, nburnin = nb,
                 nchains = nc, thin = nt, thin2 = nt2,
                 samplesAsCodaMCMC = TRUE)

