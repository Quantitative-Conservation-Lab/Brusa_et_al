#-----------#
#-Libraries-#
#-----------#

library(tidyverse)
library(nimble)
library(coda)

#-----------#
#-Functions-#
#-----------#

#Generate composition proportions
comp.fun <- function(nspecies)
{
  pi <- runif(nspecies, 0, 1)
  pi <- pi/sum(pi)
  return(pi)
}

#Generate miss ID rates
missID.fun <- function(nspecies, nmissID)
{
  phi.psi <- matrix(0, ncol = nspecies, nrow = nspecies)
  phi.psi[sample(x = 1:nspecies^2, size = nmissID)] <- runif(nmissID, 0, 0.25)
  diag(phi.psi) <- runif(nspecies, 0.8, 1)
  phi.psi <- phi.psi/apply(phi.psi, MARGIN = 1, sum)
  return(phi.psi)
}

#----------------#
#-Set parameters-#
#----------------#

#Number of species
nspecies <- 5

#Number of sites
nsites <- 500

#Community composition
pi <- comp.fun(nspecies)

#Community expected abundance
lambda.total <- runif(1, 10000, 20000)

#Species-specific abundance
lambda <- lambda.total * pi

alpha <- rnorm(nspecies, runif(1, 0.5, 1), 0.1)

#Difference in field-of-view
observer.offset <- runif(1,1,1.5)

alpha.OBS <- alpha * observer.offset

#Expected movement rate
E.alpha.OBS <- sum(pi * alpha.OBS)

#Expected miss ID rate
phi.psi <- missID.fun(nspecies, nmissID = 8)

#Detection probability
p <- runif(1, 0.25, 1)

#Derived product of movement and detection
E.epsilon <- E.alpha.OBS * p


#---------------#
#-Simulate data-#
#---------------#

#Front facing camera, point of view camera, latent correct ID abundance, latent miss ID abundance
FF <- N <- C <- matrix(NA, ncol = nspecies, nrow = nsites)

#Observer data
OBS <- array(NA, dim = c(nsites, nspecies, 2))

#Confusion matrix
confusion.matrix <- array(NA, dim = c(nsites, nspecies, nspecies))

for(j in 1:nsites){
  
  #Front facing camera data
  FF[j,] <- rpois(nspecies, lambda.total * pi)
  
  #Latent abundance w/correct ID
  N[j,] <- rpois(nspecies, lambda.total * pi * alpha.OBS)
  
  for(i in 1:nspecies){
    confusion.matrix[j,i,] <- rmultinom(1, N[j,i], phi.psi[i,])
    
  }
  
  #Latent abundance w/miss ID
  C[j,] <- apply(confusion.matrix[j,,], MARGIN = 2, sum)
  
  #Observer 1 data
  OBS[j,,1] <- rbinom(n = nspecies, size = C[j,], prob = p)
  
  #Observer 2 data
  OBS[j,,2] <- rbinom(n = nspecies, size = C[j,], prob = p)
  
}


#--------------------#
#-Estimation Model-#
#--------------------#

code <- nimbleCode({
  
  #-Priors-#
  
  #Composition of latent abundance (corrected for imperfect detection)
  phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  
  #Derived product of movement and detection
  E.epsilon ~ dnorm(0, 0.01)
  
  #-Likelihood-#
  
  for(j in 1:nsites){
    
    #Front facing camera composition
    FF[j,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[j])
    
    #Front facing camera total abundance
    FF.total[j] ~ dpois(lambda.total)
    
    for(o in 1:nobs){
      
      OBS[j,1:nspecies,o] ~ dmulti(phi[1:nspecies], OBS.total[j,o])
      
      OBS.total[j,o] ~ dpois(lambda.total * E.epsilon)
      
    }#end o
    
  }#end j
  
  for(i in 1:nspecies){
    
    pi[i] <- lambda[i]/lambda.total
    
    #psi.ones[i] <- 1
    
    phi.ones[i] <- 1
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
    correction[i] <- E.epsilon * phi[i]/pi[i]
    
  }#end i
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Compile data-#

data <- list(FF = FF, 
             FF.total = apply(FF, 1, sum),
             OBS = OBS,
             OBS.total = apply(OBS, c(1,3), sum)
)

constants <- list(nspecies = nspecies, nsites = nsites, nobs = 2)

#-Initial values-#

inits <- function(){list(pi = apply(FF/apply(FF, 1, sum), 2, mean),
                         E.epsilon = E.epsilon,
                         phi = apply(apply(OBS, c(1,2), max)/apply(apply(OBS, c(1,2), max), 1, sum), 2, mean)
)}

#-Parameters to save-#

params <- c(
  "pi",
  "phi",
  "lambda.total",
  "lambda",
  "E.epsilon",
  "correction"
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 20000
nb <- 10000
nt <- 1

#-Run model-#

out2 <- runMCMC(compiled.model$MCMC,
                niter = ni, nburnin = nb,
                nchains = nc, thin = nt,
                samplesAsCodaMCMC = TRUE)

