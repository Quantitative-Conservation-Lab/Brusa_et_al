library(tidyverse)
library(nimble)
library(coda)
library(ggmcmc)
library(coda)
library(here)


#Bring in data
sb_df <- read.csv(here('Data', 'seaducks.csv'))


#Mis-Identification by group, then by individual, then integrate mis-ID info into detection model
#Remove uninformative records
sb_df <- sb_df %>% filter(sb_df$SPECIES != "PEFA" & sb_df$SPECIES != "UNSD" & sb_df$SPECIES != "BAEA" 
                          & sb_df$SPECIES != "BLOY" & sb_df$SPECIES != "GBHE" & sb_df$SPECIES != "NWCR"
                          & sb_df$SPECIES != "HAPO" & sb_df$SPECIES != "HASE" & sb_df$SPECIES != "BLLA"
                          & sb_df$SPECIES != "UNMA" & sb_df$SPECIES !=  "UNPD" & sb_df$SPECIES !=  "UMSD" 
                          & sb_df$SPECIES != "UNDD" & sb_df$SPECIES != "UNDU" & sb_df$SPECIES !=  "USSD"
                          & sb_df$SPECIES != "UNSD" & sb_df$SPECIES != "UNSB" & sb_df$SPECIES != "RODO")

sb_df <- sb_df %>% filter(platform == 1)

#Set up array for observer records
sb_df$Transect <- str_trunc(sb_df$TransectGroup, 12, "right", ellipsis = "")
ID.obs <- sb_df

tran <- unique(ID.obs$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))

ID.obs$numeric_tran <- 0
for(i in dumb$tran){
  ID.obs$numeric_tran[ID.obs$Transect == i] <- match(i, dumb$tran)
}


# ID.obs$SPECIES <- ifelse(ID.obs$SPECIES == "SUSC", "SUSC", ifelse(ID.obs$SPECIES == "BUFF", "BUFF", 
#                                                                   ifelse(ID.obs$SPECIES == "WEGR", "WEGR", "other")))

#Run the model with species only in the 25th percentile or more for occurrence or counts
ID.obs <- ID.obs %>% filter(SPECIES != "ANMU" & SPECIES != "BAGO" & SPECIES != "BLSC" & SPECIES != "CAGO" & 
                              SPECIES != "COME" & SPECIES != "DCCO" & SPECIES != "EUWI" & SPECIES != "GWGU" &
                              SPECIES != "NOPI" & SPECIES != "PALO" & SPECIES != "RUDU" & SPECIES != "UNAC" & 
                              SPECIES != "UNGR" & SPECIES != "UNME" & SPECIES != "USAC")

ID.BM <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.BM,
              values_fn = sum, values_fill = 0)

#ID.BM <- as.matrix(ID.BM[,c(3, 4, 5, 2)])
ID.BM <- as.matrix(ID.BM[,-1])

ID.TC <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.TC,
              values_fn = sum, values_fill = 0)

#ID.TC <- as.matrix(ID.TC[,c(3, 4, 5, 2)])
ID.TC <- as.matrix(ID.TC[,-1])

OBS <- array(c(ID.BM, ID.TC), dim = c(nrow(ID.BM), ncol(ID.BM), 2))

BM.total <- apply(ID.BM, 1, sum)

TC.total <- apply(ID.TC, 1, sum)

OBS.total <- cbind(BM.total, TC.total)


ID.POV <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.POV,
              values_fn = sum, values_fill = 0)

#ID.POV <- as.matrix(ID.POV[,c(3, 4, 5, 2)])
ID.POV <- as.matrix(ID.POV[,-1])

POV <- ID.POV

POV.total <- apply(POV, 1, sum)


ID.FF <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.FF,
              values_fn = sum, values_fill = 0)

#ID.FF <- as.matrix(ID.FF[,c(3, 4, 5, 2)])
ID.FF <- as.matrix(ID.FF[,-1])

FF <- ID.FF

FF.total <- apply(FF, 1, sum)


nspecies <- length(unique(ID.obs$SPECIES))
nsites <- length(unique(ID.obs$numeric_tran))
nobs <- 2

#Build a confusion matrix based on proportions of observers to FF camera
#How to deal with zeros for FF? Maybe can set FF to 1
#How to deal with proportions > 1??
#Taking the proportion of the species over the transect

#Not sure if the following lines will become useful...

# proportions <- ID.obs
# proportions$FF.prop <- proportions$Count.FF/proportions$Count.FF
# 
# is.nan.data.frame <- function(x)
#   do.call(cbind, lapply(x, is.nan))
# 
# proportions[is.nan(proportions)] <- 1
# 
# proportions$POV.prop <- proportions$Count.POV/proportions$Count.FF
# 
# proportions$POV.prop[(!is.finite(proportions$POV.prop))] <- 0
# 
# proportions$BM.prop <- proportions$Count.BM/proportions$Count.FF
# 
# proportions$BM.prop[(!is.finite(proportions$BM.prop))] <- 0
# 
# proportions$TC.prop <- proportions$Count.TC/proportions$Count.FF
# 
# proportions$TC.prop[(!is.finite(proportions$TC.prop))] <- 0
# 
# prop.POV <- proportions %>%
#   pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = POV.prop,
#               values_fn = sum, values_fill = 0)
# 
# prop.POV <- as.matrix(prop.POV[,-1])
# 
# POV <- prop.POV
# 
# POV.total <- apply(POV, 1, sum)
# 
# prop.BM <- proportions %>%
#   pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = BM.prop,
#               values_fn = sum, values_fill = 0)
# 
# prop.BM <- as.matrix(prop.BM[,-1])
# 
# BM <- prop.BM
# 
# BM.total <- apply(BM, 1, sum)
# 
# prop.TC <- proportions %>%
#   pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = TC.prop,
#               values_fn = sum, values_fill = 0)
# 
# prop.TC <- as.matrix(prop.TC[,-1])
# 
# TC <- prop.TC
# 
# TC.total <- apply(TC, 1, sum)





#-Nimble Code-#

code <- nimbleCode({
  
  #-Priors-#
  
  #Composition of point of view camera
  psi[1:nspecies] ~ ddirch(psi.ones[1:nspecies])
  
  #Composition of latent abundance (corrected for imperfect detection)
  phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  
  #Derived product of movement and detection
  E.epsilon ~ dgamma(1, 1)
  
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
    
    pi[i] <- lambda[i]/(lambda.total)
    
    psi.ones[i] <- 1
    
    phi.ones[i] <- 1
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
    epsilon[i] <- psi[i] * E.epsilon / (pi[i])
    
  }#end i
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Compile data-#

data <- list(FF = FF, 
             FF.total = FF.total,
             POV = POV,
             POV.total = POV.total,
             OBS = OBS,
             OBS.total = OBS.total
)

constants <- list(nspecies = nspecies, nsites = nsites, nobs = 2)

#-Initial values-#

inits <- function(){list(pi = apply(FF/apply(FF, 1, sum), 2, mean, na.rm = TRUE),
                         E.epsilon = sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))
                         *runif(1, 0.25, 1)),
                         epsilon = rnorm(nspecies, runif(1, 0.5, 1), 0.1) * runif(1,1,1.5) * runif(1, 0.25, 1),
                         psi = apply(POV/apply(POV, 1, sum), 2, mean, na.rm = TRUE),
                         phi = apply(apply(OBS, c(1,2), max)/apply(apply(OBS, c(1,2), max), 1, sum), 2, mean, 
                                     na.rm = TRUE)
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

#Diagnostics
out.mcmc <- as.mcmc.list(out)
out.ggs <- ggs(out.mcmc)

out.diag <- ggs_diagnostics(out.ggs)
out.gelman <- gelman.diag(out.mcmc)
MCMCtrace(covs_sampsM, params = c("sigma0", "beta.NPGO", "beta.Depth", "beta.ho", "beta.upwell", 
                                  "sigma.eps.year", "beta.shore7", "beta.shore5", "sigma.eps.pair", 
                                  "beta.shore9A", "beta.shore1A", "beta.shore6A", "beta.FT", "beta.ST",
                                  "beta.shore4", "beta.shore6D", "beta.shore8A", "beta.river", "beta.sst",
                                  "beta.chl", "beta.sal", "beta.BM", "beta.bss.1", "beta.bss.2", "beta.bss.3",
                                  "r.N", "beta.offshore"))

ggs_traceplot(out.ggs, c("E.epsilon"))
ggs_traceplot(out.ggs, c("lambda.total"))
ggs_traceplot(out.ggs, c("phi"))
ggs_traceplot(out.ggs, c("psi."))
ggs_traceplot(out.ggs, c("pi."))

#Work with output
#Check dimensions - 20K iterations minus 10K for burn-in, 22 parameters
dim(out[[1]])
mcmc.params <- as.mcmc.list(out)
params.birds <- data.frame(as.matrix(mcmc.params)) 

params.birds = as.matrix(params.birds)

out.birds <- data.frame(#Params = c("E.epsilon", "epsilon.SUSC", "epsilon.BUFF", "Epsilon.WEGR", "Epsilon.other", 
                                   # "lambda.SUSC", "lambda.BUFF", "lambda.WEGR", "lambda.other", "lambda.total", 
                                   # "phi.SUSC", "phi.BUFF", "phi.WEGR", "phi.other", "pi.SUSC", "pi.BUFF", "pi.WEGR", 
                                   # "pi.other", "psi.SUSC", "psi.BUFF", "psi.WEGR", "psi.other"),
                        Mean = apply(params.birds, 2, mean),
                        lcl = apply(params.birds, 2, quantile, probs = c(.05)),
                        ucl = apply(params.birds, 2, quantile, probs = c(.95)),
                        SD = apply(params.birds, 2, sd))

out.summ <- MCMCsummary(mcmc.params)


