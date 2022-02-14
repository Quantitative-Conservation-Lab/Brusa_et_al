library(tidyverse)
library(stringr)
library(runjags)
library(ggmcmc)
library(ggplot2)
library(here)
library(rsample)


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
ID.obs <- sb_df[seq(1, nrow(sb_df), 2), ]
#Pull out a few unique transects/groups to make group1 and group2 have the same no. of transects and groups
ID.obs <- ID.obs[-c(1, 3, 7, 10, 11, 14),]
group1 <- ID.obs
tran <- unique(ID.obs$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))

ID.obs$numeric_tran <- 0
for(i in dumb$tran){
  ID.obs$numeric_tran[ID.obs$Transect == i] <- match(i, dumb$tran)
}

tran.grp <- unique(ID.obs$TransectGroup)
dumb1 <- data.frame(tran.grp = tran.grp,
                    numbers = 1:length(tran.grp))

ID.obs$numeric_tran.grp <- 0
for(i in dumb1$tran.grp){
  ID.obs$numeric_tran.grp[ID.obs$TransectGroup == i] <- match(i, dumb1$tran.grp)
}


ID.obs$SPECIES <- ifelse(ID.obs$SPECIES == "SUSC", "SUSC", ifelse(ID.obs$SPECIES == "BUFF", "BUFF", 
                                                                  ifelse(ID.obs$SPECIES == "WEGR", "WEGR", "other")))
ID.BM <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran, numeric_tran.grp), names_from = SPECIES, values_from = Count.BM,
              values_fn = sum, values_fill = 0)

ID.BM <- as.matrix(ID.BM[,c(4, 5, 6, 3)])

ID.TC <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran, numeric_tran.grp), names_from = SPECIES, values_from = Count.TC,
              values_fn = sum, values_fill = 0)

ID.TC <- as.matrix(ID.TC[,c(4, 5, 6, 3)])

ID <- array(c(ID.BM, ID.TC), dim = c(nrow(ID.BM), 4, 2))

M.BM <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran, numeric_tran.grp), names_from = SPECIES, values_from = Count.BM,
              values_fn = sum, values_fill = 0)

M.BM <- M.BM[,3:6]

total.species <- vector()
for(i in 1:nrow(M.BM)){
  total.species[i] <- sum(M.BM[i,])
  M.BM$total.species[i] <- total.species[i]
}

M.BM <- M.BM$total.species


M.TC <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran, numeric_tran.grp), names_from = SPECIES, values_from = Count.TC,
              values_fn = sum, values_fill = 0)

M.TC <- M.TC[,3:6]

total.species <- vector()
for(i in 1:nrow(M.TC)){
  total.species[i] <- sum(M.TC[i,])
  M.TC$total.species[i] <- total.species[i]
}

M.TC <- M.TC$total.species

M.obs <- cbind(M.BM, M.TC)


ID.POV <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran, numeric_tran.grp), names_from = SPECIES, values_from = Count.POV,
              values_fn = sum, values_fill = 0)

ID.POV <- as.matrix(ID.POV[,c(4, 5, 6, 3)])

POV <- as.matrix(ID.POV)

M.POV <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran, numeric_tran.grp), names_from = SPECIES, values_from = Count.POV,
              values_fn = sum, values_fill = 0)

M.POV <- M.POV[,3:6]

total.species <- vector()
for(i in 1:nrow(M.POV)){
  total.species[i] <- sum(M.POV[i,])
  M.POV$total.species[i] <- total.species[i]
}

M <- M.POV$total.species


group1$all <- paste(group1$TransectGroup, group1$SPECIES, group1$Count.FF, group1$Count.BM, group1$Count.POV,
                    group1$Count.TC)
sb_df$all <- paste(sb_df$TransectGroup, sb_df$SPECIES, sb_df$Count.FF, sb_df$Count.BM, sb_df$Count.POV,
                   sb_df$Count.TC)

group2 <- subset(sb_df, !(all %in% group1$all))
Counts <- group2

#Convert groups and transects to numbers
tran <- unique(Counts$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))
Counts$numeric_tran <- 0
for(i in dumb$tran){
  Counts$numeric_tran[Counts$Transect == i] <- match(i, dumb$tran)
}

tran.grp <- unique(Counts$TransectGroup)
dumb1 <- data.frame(tran.grp = tran.grp,
                    numbers = 1:length(tran.grp))

Counts$numeric_tran.grp <- 0
for(i in dumb1$tran.grp){
  Counts$numeric_tran.grp[Counts$TransectGroup == i] <- match(i, dumb1$tran.grp)
}

Counts$SPECIES <- ifelse(Counts$SPECIES == "SUSC", "SUSC", ifelse(Counts$SPECIES == "BUFF", "BUFF", 
                                                                  ifelse(Counts$SPECIES == "WEGR", "WEGR", "other")))

Counts = Counts %>% group_by(numeric_tran, numeric_tran.grp, SPECIES) %>% 
  summarise(Count.FF = sum(Count.FF), Count.BM = sum(Count.BM), Count.TC = sum(Count.TC),
            Count.POV = sum(Count.POV))

Counts.FF <- Counts %>%
  pivot_wider(id_cols = c(numeric_tran.grp), names_from = SPECIES, values_from = Count.FF,
              values_fn = sum, values_fill = 0)

Counts.FF <- Counts.FF[,c(3, 4, 5, 2)]
FF <- as.matrix(Counts.FF)

Counts.BM <- Counts %>%
  pivot_wider(id_cols = c(numeric_tran.grp), names_from = SPECIES, values_from = Count.BM,
              values_fn = sum, values_fill = 0)

Counts.BM <- Counts.BM[,c(3, 4, 5, 2)]
BM <- as.matrix(Counts.BM)

Counts.TC <- Counts %>%
  pivot_wider(id_cols = c(numeric_tran.grp), names_from = SPECIES, values_from = Count.TC,
              values_fn = sum, values_fill = 0)

Counts.TC <- Counts.TC[,c(3, 4, 5, 2)]
TC <- as.matrix(Counts.TC)

obs <- array(c(BM, TC), dim = c(nrow(BM), 4, 2))


Counts.tr <- Counts %>% group_by(numeric_tran.grp) %>% slice(1L)
tr <- Counts.tr$numeric_tran

Counts.grp <- Counts %>% group_by(numeric_tran.grp) %>% slice(1L)
n.groups <- length(unique(Counts$numeric_tran.grp))
n.transects <- length(unique(Counts$numeric_tran))
n.species <- length(unique(Counts$SPECIES))
n.observers <- 2

cat("model {

#########################SPECIES MIS-ID COMPONENT#########################

#priors on the confusion probabilities - mlogit link
for(o in 1:n.observers){
  pr.obsJ.spK[1,1,o] <- 1-pr.obsJ.spK[2,1,o]-pr.obsJ.spK[3,1,o]-pr.obsJ.spK[4,1,o]
  pr.obsJ.spK[2,1,o] <- exp(m.pr.obsJ.spK[2,1,o])/(1+exp(m.pr.obsJ.spK[2,1,o]) + exp(m.pr.obsJ.spK[3,1,o]) + 
  exp(m.pr.obsJ.spK[4,1,o]))
  pr.obsJ.spK[3,1,o] <- exp(m.pr.obsJ.spK[3,1,o])/(1+exp(m.pr.obsJ.spK[2,1,o]) + exp(m.pr.obsJ.spK[3,1,o]) + 
  exp(m.pr.obsJ.spK[4,1,o]))
  pr.obsJ.spK[4,1,o] <- exp(m.pr.obsJ.spK[4,1,o])/(1+exp(m.pr.obsJ.spK[2,1,o]) + exp(m.pr.obsJ.spK[3,1,o]) + 
  exp(m.pr.obsJ.spK[4,1,o]))

  pr.obsJ.spK[1,2,o] <- exp(m.pr.obsJ.spK[1,2,o])/(1+exp(m.pr.obsJ.spK[1,2,o]) + exp(m.pr.obsJ.spK[3,2,o]) + 
  exp(m.pr.obsJ.spK[4,2,o]))
  pr.obsJ.spK[2,2,o] <- 1-pr.obsJ.spK[1,2,o]-pr.obsJ.spK[3,2,o]-pr.obsJ.spK[4,2,o]
  pr.obsJ.spK[3,2,o] <- exp(m.pr.obsJ.spK[3,2,o])/(1+exp(m.pr.obsJ.spK[1,2,o]) + exp(m.pr.obsJ.spK[3,2,o]) + 
  exp(m.pr.obsJ.spK[4,2,o]))
  pr.obsJ.spK[4,2,o] <- exp(m.pr.obsJ.spK[4,2,o])/(1+exp(m.pr.obsJ.spK[1,2,o]) + exp(m.pr.obsJ.spK[3,2,o]) + 
  exp(m.pr.obsJ.spK[4,2,o]))

  pr.obsJ.spK[1,3,o] <- exp(m.pr.obsJ.spK[1,3,o])/(1+exp(m.pr.obsJ.spK[1,3,o]) + exp(m.pr.obsJ.spK[2,3,o]) + 
  exp(m.pr.obsJ.spK[4,3,o]))
  pr.obsJ.spK[2,3,o] <- exp(m.pr.obsJ.spK[2,3,o])/(1+exp(m.pr.obsJ.spK[1,3,o]) + exp(m.pr.obsJ.spK[2,3,o]) + 
  exp(m.pr.obsJ.spK[4,3,o]))
  pr.obsJ.spK[3,3,o] <- 1-pr.obsJ.spK[1,3,o]-pr.obsJ.spK[2,3,o]-pr.obsJ.spK[4,3,o]
  pr.obsJ.spK[4,3,o] <- exp(m.pr.obsJ.spK[4,3,o])/(1+exp(m.pr.obsJ.spK[1,3,o]) + exp(m.pr.obsJ.spK[2,3,o]) + 
  exp(m.pr.obsJ.spK[4,3,o]))

  pr.obsJ.spK[1,4,o] <- exp(m.pr.obsJ.spK[1,4,o])/(1+exp(m.pr.obsJ.spK[1,4,o]) + exp(m.pr.obsJ.spK[2,4,o]) + 
  exp(m.pr.obsJ.spK[3,4,o]))
  pr.obsJ.spK[2,4,o] <- exp(m.pr.obsJ.spK[2,4,o])/(1+exp(m.pr.obsJ.spK[1,4,o]) + exp(m.pr.obsJ.spK[2,4,o]) + 
  exp(m.pr.obsJ.spK[3,4,o]))
  pr.obsJ.spK[3,4,o] <- exp(m.pr.obsJ.spK[3,4,o])/(1+exp(m.pr.obsJ.spK[1,4,o]) + exp(m.pr.obsJ.spK[2,4,o]) + 
  exp(m.pr.obsJ.spK[3,4,o]))
  pr.obsJ.spK[4,4,o] <- 1-pr.obsJ.spK[1,4,o]-pr.obsJ.spK[2,4,o]-pr.obsJ.spK[3,4,o]

  m.pr.obsJ.spK[2,1,o] ~ dnorm(0,0.01)      
  m.pr.obsJ.spK[3,1,o] ~ dnorm(0,0.01)      
  m.pr.obsJ.spK[4,1,o] ~ dnorm(0,0.01) 
  m.pr.obsJ.spK[1,2,o] ~ dnorm(0,0.01)
  m.pr.obsJ.spK[3,2,o] ~ dnorm(0,0.01)
  m.pr.obsJ.spK[4,2,o] ~ dnorm(0,0.01)
  m.pr.obsJ.spK[1,3,o] ~ dnorm(0,0.01)
  m.pr.obsJ.spK[2,3,o] ~ dnorm(0,0.01)
  m.pr.obsJ.spK[4,3,o] ~ dnorm(0,0.01)
  m.pr.obsJ.spK[1,4,o] ~ dnorm(0,0.01)
  m.pr.obsJ.spK[2,4,o] ~ dnorm(0,0.01)
  m.pr.obsJ.spK[3,4,o] ~ dnorm(0,0.01)
}

#probability of observing species j on transect t for observer o is a function of how common species k is there and the 
#probability of observing j given true species k
#pr.spK is estimated from the POV camera, so we are just modeling pr.obsJ.spK in this part of the code  
for(t in 1:n.transects){
  for(o in 1:n.observers){
    pr.obsJ[t,1,o] <- pr.spK[t,1]*pr.obsJ.spK[1,1,o] + pr.spK[t,2]*pr.obsJ.spK[1,2,o] +
                      pr.spK[t,3]*pr.obsJ.spK[1,3,o] + pr.spK[t,4]*pr.obsJ.spK[1,4,o]

    pr.obsJ[t,2,o] <- pr.spK[t,1]*pr.obsJ.spK[2,1,o] + pr.spK[t,2]*pr.obsJ.spK[2,2,o] +
                      pr.spK[t,3]*pr.obsJ.spK[2,3,o] + pr.spK[t,4]*pr.obsJ.spK[2,4,o]

    pr.obsJ[t,3,o] <- pr.spK[t,1]*pr.obsJ.spK[3,1,o] + pr.spK[t,2]*pr.obsJ.spK[3,2,o] +
                      pr.spK[t,3]*pr.obsJ.spK[3,3,o] + pr.spK[t,4]*pr.obsJ.spK[3,4,o]

    pr.obsJ[t,4,o] <- pr.spK[t,1]*pr.obsJ.spK[4,1,o] + pr.spK[t,2]*pr.obsJ.spK[4,2,o] +
                      pr.spK[t,3]*pr.obsJ.spK[4,3,o] + pr.spK[t,4]*pr.obsJ.spK[4,4,o]

  }
}

#derived probability of species K given J was observed, Bayes rule 
#then take the mean across observers 
for(t in 1:n.transects){
  for(j in 1:n.species){
    for(k in 1:n.species){
      for(o in 1:n.observers){
        pr.spK.obsJ[t,k,j,o] <- (pr.obsJ.spK[j,k,o]*pr.spK[t,1])/(pr.obsJ[t,j,o])
      }
      pr.spK.obsJ.mean[t,k,j] <- (pr.spK.obsJ[t,k,j,1] + pr.spK.obsJ[t,k,j,2])/2 #mean per spp pair per transect
    }
  }
}


#the species ID data by observer are multinomial given the probability of observing species j and the total observed      
for(i in 1:n.groups){
  for(o in 1:n.observers){
    ID[i,,o] ~ dmultinom(pr.obsJ[tr[i],,o],M.obs[i,o])
  }
}

#########################SPECIES COMPOSITION COMPONENT(POV)#########################

#prior on species composition for transect t
for(t in 1:n.transects){
  pr.spK[t,1:4] ~ ddirch(alpha.sp)
}

  #POV camera data - multinomial sample from the total number of birds in the group, with pr.spK probabilities   
for(i in 1:n.groups){
  POV[i,] ~ dmultinom(pr.spK[tr[i],], M[i])
}
#########################AVAILABILITY COMPONENT(FF vs OBSERVERS)#########################

for(i in 1:n.groups){
  for(k in 1:n.species){
    #FF camera data - abundance of species k in group i prior to plane passing (data), modeled with a mean mu
    #might need a negative binomial here 
    FF[i,k] ~ dpois(mu[i,k])
    
    #the relationship between abundance from observers and abundance before the plane passed, by species
    mu[i,k] <- beta[k] * lambda[i,k] #N.spK instead of lambda
  }
}

#prior on beta 
for(k in 1:n.species){
  beta[k] ~ dunif(0,5) #each species gets is own beta
}  

#########################DETECTION COMPONENT#########################
  
#prior for abundance of species k in group i after plane passing (predicted from N-mix)
#might need a negative binomial here 
for(i in 1:n.groups){
  for(k in 1:n.species){
    N.obsJ[i,k] ~ dpois(lambda[i,k]) #maybe want N.spK instead of lambda, but seems circular??
    
    #prior on abundance of species j in group i#
    lambda[i,k] ~ dgamma(1,1) 
  }
}

#estimation of the abundance of putative species J - Nmixture
for(j in 1:n.species){
  for(o in 1:n.observers){
    for(i in 1:n.groups){
      obs[i,j,o] ~ dbinom(p[j,o],N.obsJ[i,j])
    }
    #prior for detection probability 
    p[j,o] ~ dbeta(1,1)
  }
}  

#########################CORRECTING FOR MIS-ID#########################
#calculate the number of species K expected for each J observed 
for(i in 1:n.groups){
  for(j in 1:n.species){
    N.spK.obsJ[i,1:4,j] ~ dmultinom(pr.spK.obsJ.mean[tr[i],1:4,j],N.obsJ[i,j])
  }
}
#sum over observed species J to get the abundance of species K 
for(i in 1:n.transects){
  for(k in 1:n.species){
    N.spK[i,k] <- sum(N.spK.obsJ[i,k,1:4])
  }
}
#take the mean of K  
for(k in 1:n.species){
  N.spK.mean[k] <- mean(N.spK[,k]) #providing a mean per transect, n.species (4) x n.transects (301)
}

#########################BASIC DETECTION COMPONENT#########################
# THIS IS THE BASIC DETECTION MODEL THAT DOESN'T ACCOUNT FOR MIS-ID 
#prior for abundance of species k in group i after plane passing (predicted from N-mix)
#might need a negative binomial here 
#for(i in 1:n.groups){
#  for(k in 1:n.species){
#    N[i,k] ~ dpois(lambda[i,k]) #N.spK instead of lambda
    
#    # prior on abundance of species j in group i
#    lambda[i,k] ~ dgamma(1,1) #N.spK instead of lambda
#  }
#}

##prior for detection probability 
#for(k in 1:n.species){
#  for(o in 1:n.observers){
#    p[k,o] ~ dbeta(1,1)
#  }
#}  
  
##observation data 
#for(i in 1:n.groups){
#  for(j in 1:n.species){
#    for(o in 1:n.observers){
#      obs[i,j,o] ~ dbin(p[j,o],N[i,j])
#    }
#  }
#}

}",file = "ducks_mult.txt")


#DATA
data <- list(ID = ID, M.obs = M.obs, obs = obs, POV = POV, M = M, FF = FF, tr = tr, n.transects = n.transects, 
             n.groups = n.groups, n.observers = n.observers, n.species = n.species, alpha.sp = c(1,1,1,1))

#INITIAL VALUES 
init.pr.spK <- matrix(NA,nrow=n.transects,ncol = n.species)
for(t in 1:n.transects){
  init.pr.spK[t,] <- c(0.25,0.25,0.25,0.25)
}
init.lambda <- matrix(NA,nrow = n.groups,ncol = n.species) #N.spK instead of lambda
for(i in 1:n.groups){
  for(j in 1:n.species){
    init.lambda[i,j] <- max(obs[i,j,])+1 #N.spK instead of lambda
  }
}

inits.mult<-function(){list(pr.spK = init.pr.spK, beta = rep(1,4), lambda = init.lambda)} #N.spK instead of lambda


#PARAMETERS
params.mult <- c("beta","p","N.spK.mean") #for diagnostics
#params.mult <- c("pr.spK","beta","p","N.spK.mean","pr.spK.obsJ.mean")

#Generate samples from the posterior distribution
mult.mod<-run.jags(model = "ducks_mult.txt",
                   monitor = params.mult,
                   data = data,
                   n.chains = 3,
                   adapt = 5000,
                   burnin = 5000,
                   sample = 5000,
                   inits = inits.mult)

#Diagnostics
mult.mod_mcmc <- as.mcmc.list(mult.mod)
mult.mcmc <- as.matrix(as.mcmc(mult.mod))
mult.mod_params <- data.frame(as.matrix(mult.mcmc)) %>%
  dplyr::select(beta.1., beta.2., beta.3., beta.4., N.spK.mean.1., N.spK.mean.2., N.spK.mean.3., N.spK.mean.4.)

mult_mcmc.slim <- as.mcmc(mult.mod_params)
mult.mcmc.slim <- as.mcmc.list(mult_mcmc.slim)

mult.mod_ggs <- ggs(mult.mod_mcmc)
ggs_geweke(mult.mod_ggs)
ggs_Rhat(mult.mod_ggs)


ggs_traceplot(mult.mod_ggs, c("beta"))
ggs_traceplot(mult.mod_ggs, c("p"))
ggs_traceplot(mult.mod_ggs, c("pr.spK.obsJ.mean"))

ggs_density(mult.mod_ggs, c("beta"))
ggs_density(mult.mod_ggs, c("p"))
ggs_density(mult.mod_ggs, c("pr.spK.obsJ.mean"))


mult.df <- data.frame(as.matrix(mult.mcmc))
mult.summ <- mult.df %>%
  group_by() %>%
  summarize_all(mean)


out.mult <- data.frame(Parameter = c("beta.1.", "beta.2.", "beta.3.", "beta.4.", "p.1.1.", "p.2.1.", "p.3.1.",
                                     "p.4.1.", "p.1.2.", "p.2.2.", "p.3.2.", "p.4.2.", "N.spK.mean.1.", "N.spK.mean.2.",
                                     "N.spK.mean.3.", "N.spK.mean.4."),
                       Mean = apply(mult.mcmc, 2, mean),
                       lcl = apply(mult.mcmc, 2, quantile, probs = c(.05)),
                       ucl = apply(mult.mcmc, 2, quantile, probs = c(.95)))

