#-Packages-#

library(tidyverse)
library(nimble)

#-Functions-#

#Generate composition proportions
comp.fun <- function(nspecies)
{
  pi <- runif(nspecies, 0, 1)
  pi <- pi/sum(pi)
  return(pi)
}

#Generate miss ID rates
missID.fun <- function(nspecies)
{
  phi.psi <- matrix(runif(nspecies^2,0,0.25),ncol = nspecies, nrow = nspecies)
  diag(phi.psi) <- runif(nspecies, 0.8, 1)
  phi.psi <- phi.psi/apply(phi.psi, MARGIN = 1, sum)
  return(phi.psi)
}

#-Set Parameters-#

#Number of species
nspecies <- 10

#Number of groups
ngroups <- 1000

#Community composition
pi <- comp.fun(nspecies)

#Community epected abundance
lambda.total <- runif(1, 500, 1000)

#Species expected abundance
lambda <- lambda.total * pi

#Movement rate
mu.alpha.POV <- runif(1, 0.5, 1)
sd.alpha.POV <- 0.1

mu.alpha.OBS <- mu.alpha.POV * runif(1,1,1.5)
sd.alpha.OBS <- 0.1

alpha.POV <- rnorm(nspecies, mu.alpha.POV, sd.alpha.POV)
alpha.OBS <- rnorm(nspecies, mu.alpha.OBS, sd.alpha.OBS)

#Expected movement rate
E.alpha.POV <- sum(pi * alpha.POV)

#Expected community composition (post-contact)
psi <- alpha.POV * pi / E.alpha.POV 

#Expected miss ID rate
phi.psi <- missID.fun(nspecies)

#Law of total probability
phi <- t(phi.psi) %*% psi

#Bayes rule
psi.phi <- matrix(NA, nrow = nspecies, ncol = nspecies)
for(j in 1:nspecies){
  for(k in 1:nspecies){
    psi.phi[j,k] <- psi[j] * phi.psi[j,k] / phi[k]
  }
}

#Detection probability
p <- runif(1, 0.5, 1)

#-Simulate data-#

#Front facing camera, point of view camera, latent correct ID abundance, latent miss ID abundance
FF <- POV <- N <- C <- matrix(NA, ncol = nspecies, nrow = ngroups)

#Observer data
OBS <- array(NA, dim = c(ngroups, nspecies, 2))

for(i in 1:ngroups){
  
  #Front facing camera data
  FF[i,] <- rpois(nspecies, lambda)
  
  #Point of view camera data
  POV[i,] <- rpois(nspecies, lambda * alpha.POV)
  
  #Latent abundance w/correct ID
  N[i,] <- rpois(nspecies, lambda * alpha.OBS)
  
  confusion.matrix <- NULL
  
  for(j in 1:nspecies){
    
    confusion.matrix <- cbind(confusion.matrix, rmultinom(1, N[i,j], phi.psi[j,]))
    
  }
  
  #Latent abundance w/miss ID
  C[i,] <- apply(confusion.matrix, MARGIN = 1, sum)
  
  #Observer 1 data
  OBS[i,,1] <- rbinom(n = nspecies, size = C[i,], prob = p)
  
  #Observer 2 data
  OBS[i,,2] <- rbinom(n = nspecies, size = C[i,], prob = p)
  
}

#-Nimble Code-#

code <- nimbleCode({
  
  #-Priors-#
  
  #Detection probability
  p ~ dunif(0, 1)
  
  #Movement rate / availability for point of view camera
  mu.alpha.POV ~ dnorm(0, 0.01)
  
  #Movement rate / availability for observers
  mu.alpha.OBS ~ dnorm(0, 0.01)
  
  #Composition of point of view camera
  psi[1:nspecies] ~ ddirch(psi.ones[1:nspecies])
  
  #Composition of latent abundance (corrected for imperfect detection)
  phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  
  #-Likelihood-#
  
  for(i in 1:ngroups){
    
    #Front facing camera composition
    FF[i,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[i])
    
    #Point of view camera composition
    POV[i,1:nspecies] ~ dmulti(psi[1:nspecies], POV.total[i])
    
    #Composition of latent abundance (corrected for imperfect detection)
    C[i,1:nspecies] ~ dmulti(phi[1:nspecies], N.total[i])
    
    #Front facing camera total abundance
    FF.total[i] ~ dpois(lambda.total)
    
    #Point of view camera total abundance
    POV.total[i] ~ dpois(lambda.total * mu.alpha.POV)
    
    #Total latent abundance (corrected for imperfect detection)
    N.total[i] ~ dpois(lambda.total * mu.alpha.OBS)
    
    for(j in 1:nspecies){

      for(o in 1:nobs){

        OBS[i,j,o] ~ dbin(p, C[i,j])

      }#end o

    }#end j
    
  }#end i
  
  for(j in 1:nspecies){
    
    pi[j] <- lambda[j]/lambda.total
    
    psi.ones[j] <- 1
    alpha.POV[j] <- psi[j] * mu.alpha.POV/pi[j]
    
    phi.ones[j] <- 1
    alpha.OBS[j] <- phi[j] * mu.alpha.OBS/pi[j]
    
    log(lambda[j]) <- lambda0[j]
    lambda0[j] ~ dnorm(0, 0.01)
    
  }#end j
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Compile data-#

data <- list(FF = FF, 
             FF.total = apply(FF, 1, sum),
             POV = POV,
             POV.total = apply(POV, 1, sum),
             OBS = OBS
             )

constants <- list(nspecies = nspecies, ngroups = ngroups, nobs = 2)

#-Initial values-#

inits <- function(){list(pi = apply(FF/apply(FF, 1, sum), 2, mean),
                         p = p,
                         alpha.POV = alpha.POV,
                         alpha.OBS = alpha.OBS,
                         mu.alpha.POV = mu.alpha.POV,
                         mu.alpha.OBS = mu.alpha.OBS,
                         psi = apply(POV/apply(POV, 1, sum), 2, mean),
                         phi = apply(apply(OBS, c(1,2), max)/apply(apply(OBS, c(1,2), max), 1, sum), 2, mean),
                         C = apply(OBS, c(1,2), max),
                         N.total = apply(apply(OBS, c(1,2), max), 1, sum)
)}

#-Parameters to save-#

params <- c(
            "p",
            "alpha.POV",
            "alpha.OBS",
            "mu.alpha.POV",
            "mu.alpha.OBS",
            "pi",
            "psi",
            "phi",
            "lambda",
            "lambda.total"
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

out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)



