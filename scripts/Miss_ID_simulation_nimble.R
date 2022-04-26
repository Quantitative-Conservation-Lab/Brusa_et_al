library(tidyverse)
library(nimble)

# #Logit function
# logit <- function(pp)
# {
#   log(pp) - log(1-pp)
# }
# 
# #Inverse logit function
# expit <- function(eta)
# {
#   1/(1+exp(-eta))
# }

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
  #phi.psi <- matrix(runif(nspecies^2,0,0.25),ncol = nspecies, nrow = nspecies)
  phi.psi <- matrix(0, ncol = nspecies, nrow = nspecies)
  phi.psi[sample(x = 1:nspecies^2, size = nmissID)] <- runif(nmissID, 0, 0.25)
  diag(phi.psi) <- runif(nspecies, 0.8, 1)
  phi.psi <- phi.psi/apply(phi.psi, MARGIN = 1, sum)
  return(phi.psi)
}

#Number of species
nspecies <- 5

#Number of groups
ngroups <- 2000

#Community composition
pi <- comp.fun(nspecies)

#Community epected abundance
lambda.total <- runif(1, 300, 1000)

#Species expected abundance
lambda <- lambda.total * pi

#Availability probability
alpha <- runif(nspecies, 0.5, 1)

#Mean availability probability
#mu.alpha <- mean(alpha)
mu.alpha <- sum(pi * alpha)

#Expected community composition (post-contact)
psi <- alpha * pi / mu.alpha

#Expected miss ID rate
phi.psi <- missID.fun(nspecies, nmissID = 4)

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
p <- runif(1, 0, 1)

#Front facing and point of view camera data
POV <- FF <- matrix(NA, ncol = nspecies, nrow = ngroups)

#Latent miss ID count
C <- matrix(NA, ncol = nspecies, nrow = ngroups)

#Observer data
OBS <- array(NA, dim = c(ngroups, nspecies, 2))

for(i in 1:ngroups){
  
  FF[i,] <- rpois(nspecies, lambda) 
  
  POV[i,] <- rbinom(n = nspecies, size = FF[i,], prob = alpha)
  
  tmp <- NULL #Maybe save these confusion matrix values 
  
  for(j in 1:nspecies){
    
    tmp <- cbind(tmp, rmultinom(1, POV[i,j], phi.psi[j,]))
    
  }
  
  C[i,] <- apply(tmp, MARGIN = 1, sum)
  
  OBS[i,,1] <- rbinom(n = nspecies, size = C[i,], prob = p)
  
  OBS[i,,2] <- rbinom(n = nspecies, size = C[i,], prob = p)
  
}


code <- nimbleCode({
  
  #Priors
  
  #Detection probability
  #p ~ dunif(0, 1)
  
  #Mean availability probability
  # mu.alpha ~ dunif(0, 1)
  
  #pi[1:nspecies] ~ ddirch(pi.c[1:nspecies])
  psi[1:nspecies] ~ ddirch(psi.c[1:nspecies])
  
  #Derived product of movement and detection
  E.epsilon ~ dnorm(0, 0.01)
  
  # phi.psi[1,1] <- 1 - sum(phi.psi[2:nspecies,1])
  phi.psi[1,1] <- 1 - sum(phi.psi[1,2:nspecies])
  
  # phi.psi[2,2] <- 1 - phi.psi[1,2] - sum(phi.psi[3:nspecies,2])
  # phi.psi[(nspecies-1),(nspecies-1)] <- 1 - sum(phi.psi[1:(nspecies-2),(nspecies-1)]) - phi.psi[nspecies,(nspecies-1)]
  
  #phi.psi[nspecies,nspecies] <- 1 - sum(phi.psi[1:(nspecies-1),nspecies])
  phi.psi[nspecies,nspecies] <- 1 - sum(phi.psi[nspecies,1:(nspecies-1)])

  for(j in 2:(nspecies-1)){

    # phi.psi[j,j] <- 1 - sum(phi.psi[1:(j-1),j]) - sum(phi.psi[(j+1):nspecies,j])
    phi.psi[j,j] <- 1 - sum(phi.psi[j,1:(j-1)]) - sum(phi.psi[j,(j+1):nspecies])
    
  }

  for(k in 1:nspecies){

    for(j in (k+1):nspecies){

      phi.psi[j,k] <- exp(m.phi.psi[j,k])/(1+exp(m.phi.psi[j,k]) + exp(m.phi.psi[j,k]) + exp(m.phi.psi[j,k]))
      
      m.phi.psi[j,k] ~ dnorm(0, 0.01)

    }#end j


    for(j in 1:(k-1)){

      phi.psi[j,k] <- exp(m.phi.psi[j,k])/(1+exp(m.phi.psi[j,k]) + exp(m.phi.psi[j,k]) + exp(m.phi.psi[j,k]))

      m.phi.psi[j,k] ~ dnorm(0, 0.01)

    }#end j

  }# end k
  
  for(j in 1:nspecies){
    
    #pi.c[j] <- 1
    
    pi[j] <- lambda[j]/lambda.total
    
    psi.ones[j] <- 1
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
  }#end j
  
  lambda.total <- sum(lambda[1:nspecies])
  
  #Likelihood
  
  for(i in 1:ngroups){
    
    #Availability of birds post aircraft contact
    #POV.total[i] ~ dbin(mu.alpha, FF.total[i])
    
    #Species composition prior to aircraft contact
    FF[i,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[i])
    
    #Front facing camera total abundance
    FF.total[i] ~ dpois(lambda.total)
    
    #Species composition post aircraft contact
    POV[i,1:nspecies] ~ dmulti(psi[1:nspecies], POV.total[i])
    
    for(o in 1:nobs){
      
      OBS[i,1:nspecies,o] ~ dmulti(phi[1:nspecies], OBS.total[i,o])
      
      OBS.total[i,o] ~ dpois(lambda.total * E.epsilon)
      
    }#end o
    
    # for(j in 1:nspecies){
    #   
    #   for(o in 1:nobs){
    #     
    #     #Observation process (imperfect detection)
    #     OBS[i,j,o] ~ dbin(p * phi[j], POV.total[i])
    #     
    #   }#end o
    #   
    # }#end j
    
  }#end i
  
  for(j in 1:nspecies){
    
    for(k in 1:nspecies){
      
      psi.phi[j,k] <- psi[j] * phi.psi[j,k] / phi[k]
      
      phi.hold[j,k] <- psi[j] * phi.psi[j,k]
      
      beta[j,k] <- phi[k]/psi[j]
      
    }#end k
    
    #Species-specific availability
    #alpha[j] <- psi[j] * mu.alpha / pi[j]
    
    #Species proportion post aircraft contact
    phi[j] <- sum(phi.hold[1:nspecies,j]) #this should be in a k loop for clarity
    
  }#end j
  
})

#-Informed data-#
phi.psi.informed <- psi.phi.informed <- matrix(NA, nrow = nspecies, ncol = nspecies)

for(i in 1:nspecies){
  for(k in 1:nspecies){
    if(phi.psi[i,k] < 0.025){
      phi.psi.informed[i,k] <- 0
    }
    if(psi.phi[i,k] < 0.025){
      psi.phi.informed[i,k] <- 0
    }
  }
}

#-Compile data-#

data <- list(FF = FF, 
             FF.total = apply(FF, 1, sum),
             POV = POV, 
             #POV.total = apply(POV, 1, sum),
             OBS = OBS,
             OBS.total = apply(OBS, c(1,3), sum),
             #phi.psi = phi.psi.informed,
             psi.phi = psi.phi)

constants <- list(nspecies = nspecies, ngroups = ngroups, nobs = 2)

#-Initial values-#

inits <- function(){list(#mu.alpha = mean(apply(POV/FF, 2, mean)),
                         p = p,
                         #pi = apply(FF/apply(FF, 1, sum), 2, mean),
                         psi = apply(POV/apply(POV, 1, sum), 2, mean))}

#-Parameters to save-#

params <- c("p",
            #"pi", 
            "psi", 
            "phi",
            #"alpha", 
            #"mu.alpha", 
            "beta",
            "psi.phi",
            "phi.psi"
            
)

#-MCMC settings-#

model <- nimbleModel(code = code,
                     constants = constants,
                     data = data,
                     inits = inits())

MCMCconf <- configureMCMC(model, monitors = params)

# MCMCconf$addSampler(target = c("m.phi.psi[2, 1]", "m.phi.psi[3, 1]", "m.phi.psi[4, 1]"),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c("m.phi.psi[1, 2]", "m.phi.psi[3, 2]", "m.phi.psi[4, 2]"),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c("m.phi.psi[1, 3]", "m.phi.psi[2, 3]", "m.phi.psi[4, 3]"),
#                     type = "RW_block")
# 
# MCMCconf$addSampler(target = c("m.phi.psi[1, 4]", "m.phi.psi[2, 4]", "m.phi.psi[3, 4]"), 
#                     type = "RW_block")


MCMC <- buildMCMC(MCMCconf)

compiled.model <- compileNimble(model, MCMC)

nc <- 3
ni <- 20000
nb <- 10000
nt <- 1

#-Run model-#

out1 <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#-Convergence-#
gelman.diag(out)

plot(out)