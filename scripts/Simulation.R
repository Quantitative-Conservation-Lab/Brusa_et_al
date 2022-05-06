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
  #phi.psi <- matrix(runif(nspecies^2,0,0.25),ncol = nspecies, nrow = nspecies)
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
nsites <- 1000

#Community composition
pi <- comp.fun(nspecies)

#Community expected abundance
lambda.total <- runif(1, 10000, 20000)

#Species-specific abundance
lambda <- lambda.total * pi

#Movement rate
mu.alpha.POV <- runif(1, 0.5, 1)
sd.alpha.POV <- 0.1

#Difference in field-of-view
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
phi.psi <- missID.fun(nspecies, nmissID = 8)

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

#---------------#
#-Simulate data-#
#---------------#

#Front facing camera, point of view camera, latent correct ID abundance, latent miss ID abundance
FF <- POV <- N <- C <- matrix(NA, ncol = nspecies, nrow = nsites)

#Observer data
OBS <- array(NA, dim = c(nsites, nspecies, 2))

#Confusion matrix
confusion.matrix <- array(NA, dim = c(nsites, nspecies, nspecies))

for(j in 1:nsites){
  
  #Front facing camera data
  FF[j,] <- rpois(nspecies, lambda.total * pi)
  
  #Point of view camera data
  POV[j,] <- rpois(nspecies, lambda.total * pi * alpha.POV)
  
  #Latent abundance w/correct ID
  N[j,] <- rpois(nspecies, lambda.total * pi * alpha.OBS)
  
  for(i in 1:nspecies){
    
    # confusion.matrix[j,i,] <- cbind(confusion.matrix, rmultinom(1, N[j,i], phi.psi[i,]))
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
#-Estimation Model 1-#
#--------------------#

code <- nimbleCode({
  
  #Priors
  
  #Detection probability
  p ~ dunif(0, 1)
  
  #Mean availability probability
  E.alpha ~ dnorm(0, 0.01)
  
  psi[1:nspecies] ~ ddirch(psi.ones[1:nspecies])
  
  #Derived product of movement and detection
  E.epsilon <- p * E.alpha
  
  #Conditional probabilities
  for(i in 1:nspecies){
    
    phi.psi[i,1:nspecies] ~ ddirch(phi.psi.ones[i,1:nspecies])
    
    phi.psi.ones[i,1:nspecies] <- 1
    
    pi[i] <- lambda[i]/lambda.total
    
    psi.ones[i] <- 1
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
  }#end i
  
  lambda.total <- sum(lambda[1:nspecies])
  
  #Likelihood
  
  for(j in 1:nsites){
    
    #Front facing camera total abundance
    FF.total[j] ~ dpois(lambda.total)
    
    #Latent available abundance
    N.total[j] ~ dpois(lambda.total * E.alpha)
    
    #Species composition prior to aircraft contact
    FF[j,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[j])
    
    #Point of view camera composition
    POV[j,1:nspecies] ~ dmulti(psi[1:nspecies], POV.total[j])
    
    for(i in 1:nspecies){
      
      #Latent species composition  post aircraft contact
      N[j,i] ~ dbin(psi[i], N.total[j])
      
      for(k in 1:nspecies){
       
        #Confusion matrix
        confusion.matrix[j,i,k] ~ dbin(phi.psi[i,k], N[j,i])
         
      }#end k
      
    }#end i
    
    for(o in 1:nobs){
      
      OBS.total[j,o] ~ dpois(lambda.total * E.epsilon)
      
      for(k in 1:nspecies){
        
        OBS[j,k,o] ~ dbin(p, C[j,k])
        
      }#end k
      
    }#end o
    
    for(k in 1:nspecies){
      
      C[j,k] <- sum(confusion.matrix[j,1:nspecies,k])
      
    }#end k
    
  }#end j
  
  for(i in 1:nspecies){
    
    for(k in 1:nspecies){
      
      psi.phi[i,k] <- psi[i] * phi.psi[i,k] / phi[k]
      
      phi.hold[i,k] <- psi[i] * phi.psi[i,k]
      
    }#end k
    
    #Species-specific availability
    alpha[i] <- psi[i] * E.alpha / pi[i]
    
    #Species-specific movement/detection
    epsilon[i] <- psi[i] * E.epsilon / pi[i]
    
    #Species-specific misidentification
    missID[i] <- phi[i]/psi[i]
    
  }#end j
  
  for(k in 1:nspecies){
    
    #Species proportion post aircraft contact
    phi[k] <- sum(phi.hold[1:nspecies,k])
    
  }#end k
  
})


#-Informed data-#
# phi.psi.informed <- psi.phi.informed <- matrix(NA, nrow = nspecies, ncol = nspecies)
# 
# for(i in 1:nspecies){
#   for(k in 1:nspecies){
#     if(phi.psi[i,k] < 0.025){
#       phi.psi.informed[i,k] <- 0
#     }
#     if(psi.phi[i,k] < 0.025){
#       psi.phi.informed[i,k] <- 0
#     }
#   }
# }

#-Compile data-#

data <- list(FF = FF[1:(nsites/2),], 
             FF.total = apply(FF[1:(nsites/2),], 1, sum),
             POV = POV[1:(nsites/2),],
             POV.total = apply(POV[1:(nsites/2),], 1, sum),
             OBS = OBS[1:(nsites/2),,],
             OBS.total = apply(OBS[1:(nsites/2),,], c(1,3), sum)
)

constants <- list(nspecies = nspecies, nsites = nsites/2, nobs = 2)

#-Initial values-#
inits <- function(){list(pi = apply(FF[1:(nsites/2),]/apply(FF[1:(nsites/2),], 1, sum), 2, mean),
                         E.epsilon = E.epsilon,
                         epsilon = epsilon,
                         psi = apply(POV[1:(nsites/2),]/apply(POV[1:(nsites/2),], 1, sum), 2, mean),
                         confusion.matrix = confusion.matrix[1:(nsites/2),,],
                         C = C[1:(nsites/2),],
                         N = N[1:(nsites/2),],
                         N.total = apply(N[1:(nsites/2),], 1, sum)
)}

#-Parameters to save-#

params <- c(
  "p", 
  "psi", 
  "phi",
  "E.alpha", 
  "lambda.total",
  "lambda",
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

out1 <- runMCMC(compiled.model$MCMC,
                niter = ni, nburnin = nb,
                nchains = nc, thin = nt,
                samplesAsCodaMCMC = TRUE)

#--------------------#
#-Estimation Model 2-#
#--------------------#

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
    
    missID[i] <- phi[i]/psi[i]
    
  }#end i
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Compile data-#

data <- list(FF = FF[1:(nsites/2),], 
             FF.total = apply(FF[1:(nsites/2),], 1, sum),
             POV = POV[1:(nsites/2),],
             POV.total = apply(POV[1:(nsites/2),], 1, sum),
             OBS = OBS[1:(nsites/2),,],
             OBS.total = apply(OBS[1:(nsites/2),,], c(1,3), sum)
)

constants <- list(nspecies = nspecies, nsites = nsites/2, nobs = 2)

#-Initial values-#

inits <- function(){list(pi = apply(FF[1:(nsites/2),]/apply(FF[1:(nsites/2),], 1, sum), 2, mean),
                         E.epsilon = E.epsilon,
                         epsilon = epsilon,
                         psi = apply(POV[1:(nsites/2),]/apply(POV[1:(nsites/2),], 1, sum), 2, mean),
                         phi = apply(apply(OBS[1:(nsites/2),,], c(1,2), max)/apply(apply(OBS[1:(nsites/2),,], c(1,2), max), 1, sum), 2, mean)
)}

#-Parameters to save-#

params <- c(
            "pi",
            "psi",
            "phi",
            "lambda.total",
            "lambda",
            "E.epsilon",
            "missID"
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

#---------------#
#-Out of sample-#
#---------------#

code <- nimbleCode({
  
  E.epsilon ~ dnorm(mean.E.epsilon, sd = sd.E.epsilon)

  for(j in 1:nsites){
    
    for(o in 1:nobs){
      
      for(i in 1:nspecies){
        
        OBS[j,i,o] ~ dpois(lambda[i] * missID[i] * epsilon[i])
        
      }
      
    }#end o
    
  }#end j
  
  for(i in 1:nspecies){
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
    epsilon[i] ~ dnorm(mean.epsilon[i], sd = sd.epsilon[i])
    missID[i] ~ dnorm(mean.missID[i], sd = sd.missID[i])

  }#end i
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Informed priors-#

mean.E.epsilon <- summary(out2)[[1]]["E.epsilon","Mean"]
sd.E.epsilon <- summary(out2)[[1]]["E.epsilon","SD"]

mean.epsilon <- summary(out2)[[1]][grepl("epsilon\\[", attr(summary(out2)[[1]], "dimnames")[[1]]), "Mean"]
sd.epsilon <- summary(out2)[[1]][grepl("epsilon\\[", attr(summary(out2)[[1]], "dimnames")[[1]]), "SD"]

mean.missID <- summary(out2)[[1]][grepl("missID\\[", attr(summary(out2)[[1]], "dimnames")[[1]]), "Mean"]
sd.missID <- summary(out2)[[1]][grepl("missID\\[", attr(summary(out2)[[1]], "dimnames")[[1]]), "SD"]

#-Compile data-#

data <- list(OBS = OBS[(nsites/2+1):nsites,,],
             mean.E.epsilon = mean.E.epsilon,
             sd.E.epsilon = sd.E.epsilon,
             mean.epsilon = mean.epsilon,
             sd.epsilon = sd.epsilon,
             mean.missID = mean.missID, 
             sd.missID = sd.missID
)

constants <- list(nspecies = nspecies, nsites = nsites/2, nobs = 2)

#-Initial values-#

inits <- function(){list(lambda.total = lambda.total,
                         lambda = lambda.total * pi
)}

#-Parameters to save-#

params <- c(
  "lambda.total",
  "lambda"
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

out3 <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)

#--------#
#-Output-#
#--------#


output <- data.frame(matrix(NA, nrow = 8, ncol = 3))

output[1:3,1] <- "Model1"
output[4:6,1] <- "Model2"
output[7:8,1] <- "OOS"

output[1,2] <- (summary(out1)[[1]]["p","Mean"] - p)/p
output[2,2] <- (summary(out1)[[1]]["E.alpha","Mean"] - E.alpha.OBS)/E.alpha.OBS
output[3,2] <- (summary(out1)[[1]]["E.epsilon","Mean"] - E.epsilon)/E.epsilon
output[4,2] <- (summary(out2)[[1]]["E.epsilon","Mean"] - E.epsilon)/E.epsilon
output[5,2] <- (summary(out2)[[1]]["lambda.total","Mean"] - lambda.total)/lambda.total
output[6,2] <- mean((summary(out2)[[1]][grepl("lambda\\[", attr(summary(out2)[[1]], "dimnames")[[1]]), "Mean"] - lambda)/lambda)
output[7,2] <- (summary(out3)[[1]]["lambda.total","Mean"] - lambda.total)/lambda.total
output[8,2] <- mean((summary(out3)[[1]][grepl("lambda\\[", attr(summary(out3)[[1]], "dimnames")[[1]]), "Mean"] - lambda)/lambda)

output[1,3] <- gelman.diag(out1[c(1:3)][,"p"])$psrf[,1] < 1.1
output[2,3] <- gelman.diag(out1[c(1:3)][,"E.alpha"])$psrf[,1] < 1.1
output[3,3] <- gelman.diag(out1[c(1:3)][,"E.epsilon"])$psrf[,1] < 1.1
output[4,3] <- gelman.diag(out2[c(1:3)][,"E.epsilon"])$psrf[,1] < 1.1
output[5,3] <- gelman.diag(out2[c(1:3)][,"lambda.total"])$psrf[,1] < 1.1
output[6,3] <- all(gelman.diag(out2[c(1:3)][,grepl("lambda\\[", attr(summary(out2)[[1]], "dimnames")[[1]])])$psrf[,1] < 1.1)
output[7,3] <- gelman.diag(out3[c(1:3)][,"lambda.total"])$psrf[,1] < 1.1
output[8,3] <- all(gelman.diag(out3[c(1:3)][,grepl("lambda\\[", attr(summary(out3)[[1]], "dimnames")[[1]])])$psrf[,1] < 1.1)

#-------------#
#-Save output-#
#-------------#

ID <- length(list.files("./output/")) + 1
save(output, file = paste("./output/output", ID, ".Rds", sep=""))
