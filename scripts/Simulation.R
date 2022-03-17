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
nspecies <- 5

#Number of sites
nsites <- 500

#Community composition
pi <- comp.fun(nspecies)

#Community expected abundance
lambda.total <- runif(1, 10000, 20000)

#Movement rate
mu.alpha.POV <- runif(1, 0.5, 1)
sd.alpha.POV <- 0.1

observer.offset <- runif(1,1,1.5)

alpha.POV <- rnorm(nspecies, mu.alpha.POV, sd.alpha.POV)
alpha.OBS <- alpha.POV * observer.offset

#Expected movement rate
E.alpha.POV <- sum(pi * alpha.POV)

#Expected movement rate
E.alpha.OBS <- sum(pi * alpha.OBS)

#Expected community composition (post-contact)
psi <- alpha.POV * pi / E.alpha.POV 

#Expected miss ID rate
phi.psi <- missID.fun(nspecies)

#Law of total probability
phi <- t(phi.psi) %*% psi

#Bayes rule
psi.phi <- matrix(NA, nrow = nspecies, ncol = nspecies)
for(i in 1:nspecies){
  for(k in 1:nspecies){
    psi.phi[i,k] <- psi[i] * phi.psi[i,k] / phi[k]
  }
}

#Detection probability
p <- runif(1, 0.25, 1)

#Derived product of movement and detection
E.epsilon <- E.alpha.OBS * p

#Species-specific product of movement and detection
epsilon <- alpha.OBS * p

#-Simulate data-#

#Front facing camera, point of view camera, latent correct ID abundance, latent miss ID abundance
FF <- POV <- N <- C <- matrix(NA, ncol = nspecies, nrow = nsites)

#Observer data
OBS <- array(NA, dim = c(nsites, nspecies, 2))

for(j in 1:nsites){
  
  #Front facing camera data
  FF[j,] <- rpois(nspecies, lambda.total * pi)
  
  #Point of view camera data
  POV[j,] <- rpois(nspecies, lambda.total * pi * alpha.POV)
  
  #Latent abundance w/correct ID
  N[j,] <- rpois(nspecies, lambda.total * pi * alpha.OBS)
  
  confusion.matrix <- NULL
  
  for(i in 1:nspecies){
    
    confusion.matrix <- cbind(confusion.matrix, rmultinom(1, N[j,i], phi.psi[i,]))
    
  }
  
  #Latent abundance w/miss ID
  C[j,] <- apply(confusion.matrix, MARGIN = 1, sum)
  
  #Observer 1 data
  OBS[j,,1] <- rbinom(n = nspecies, size = C[j,], prob = p)
  
  #Observer 2 data
  OBS[j,,2] <- rbinom(n = nspecies, size = C[j,], prob = p)
  
}

#-Nimble Code-#

code <- nimbleCode({
  
  #-Priors-#
  
  #Composition of point of view camera
  psi[1:nspecies] ~ ddirch(psi.ones[1:nspecies])
  
  #Composition of latent abundance (corrected for imperfect detection)
  phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  
  #Derived product of movement and detection
  E.epsilon ~ dnorm(0, 0.01)
  
  #-Likelihood-#
  
  for(j in 1:nsites){
    
    #Front facing camera composition
    FF[j,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[j])
    
    #Point of view camera composition
    POV[j,1:nspecies] ~ dmulti(psi[1:nspecies], POV.total[j])
    
    #Front facing camera total abundance
    FF.total[j] ~ dpois(lambda.total)
    
    for(o in 1:nobs){
      
      OBS[j,1:nspecies,o] ~ dmulti(phi[1:nspecies], OBS.total[j,o])
      
      OBS.total[j,o] ~ dpois(lambda.total * E.epsilon)
      
    }#end o
    
  }#end j
  
  for(i in 1:nspecies){
    
    pi[i] <- lambda[i]/lambda.total
    
    psi.ones[i] <- 1
    
    phi.ones[i] <- 1
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
    epsilon[i] <- psi[i] * E.epsilon / pi[i]
    
  }#end i
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Compile data-#

data <- list(FF = FF, 
             FF.total = apply(FF, 1, sum),
             POV = POV,
             POV.total = apply(POV, 1, sum),
             OBS = OBS,
             OBS.total = apply(OBS, c(1,3), sum)
             )

constants <- list(nspecies = nspecies, nsites = nsites, nobs = 2)

#-Initial values-#

inits <- function(){list(pi = apply(FF/apply(FF, 1, sum), 2, mean),
                         E.epsilon = E.epsilon,
                         epsilon = epsilon,
                         psi = apply(POV/apply(POV, 1, sum), 2, mean),
                         phi = apply(apply(OBS, c(1,2), max)/apply(apply(OBS, c(1,2), max), 1, sum), 2, mean)
)}

#-Parameters to save-#

params <- c(
            "pi",
            "psi",
            "phi",
            "lambda",
            "lambda.total",
            "epsilon",
            "E.epsilon"
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

#-Output-#

#plot(out[1:3][,!grepl("N", attr(out$chain1, "dimnames")[[2]])])

#plot(out[1:3][,grep("alpha.OBS\\[1\\]", attr(out$chain1, "dimnames")[[2]])])

# output <- round(cbind(unlist(sapply(attr(summary(out)$statistics, "dimnames")[[1]], function(x) eval(parse(text=x)))),
#                 summary(out)$statistics[,"Mean"],
#                 summary(out)$quantile[,c(1,5)]), digits = 3)

