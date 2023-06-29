#---------------#
#-Out of sample-#
#---------------#

#Community expected abundance
lambda.total.oos <- runif(1, 10000, 20000)

#Community composition
pi.oos <- comp.fun(nspecies)

#Species-specific abundance
lambda.oos <- lambda.total.oos * pi.oos

#Out of sample data
N <- C <- matrix(NA, ncol = nspecies, nrow = nsites)

#Observer data
OBS <- array(NA, dim = c(nsites, nspecies, 2))

#Confusion matrix
confusion.matrix <- array(NA, dim = c(nsites, nspecies, nspecies))

for(j in 1:nsites){
  
  #Latent abundance w/correct ID
  N[j,] <- rpois(nspecies, lambda.total.oos * pi.oos * alpha.OBS)
  
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

code <- nimbleCode({
  
  for(j in 1:nsites){
    
    for(o in 1:nobs){
      
      for(i in 1:nspecies){
        
        
        OBS[j,i,o] ~ dpois(lambda[i] * correction[i])
        
      }
      
    }#end o
    
  }#end j
  
  for(i in 1:nspecies){
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
    
    correction[i] ~ dnorm(mean.correction[i], sd = sd.correction[i])
    
  }#end i
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Informed priors-#

mean.correction <- summary(out2)[[1]][grepl("correction\\[", attr(summary(out2)[[1]], "dimnames")[[1]]), "Mean"]

sd.correction <- summary(out2)[[1]][grepl("correction\\[", attr(summary(out2)[[1]], "dimnames")[[1]]), "SD"]

#-Compile data-#

data <- list(OBS = OBS,
             mean.correction = mean.correction,
             sd.correction = sd.correction
)

constants <- list(nspecies = nspecies, nsites = nsites, nobs = 2)

#-Initial values-#

inits <- function(){list(lambda.total = lambda.total.oos,
                         lambda = lambda.oos
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


