library(tidyverse)
library(stringr)
library(runjags)
library(ggmcmc)
library(ggplot2)
library(here)


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

ID <- array(c(ID.BM, ID.TC), dim = c(1758, 4, 2))

Counts <- sb_df

#Convert groups and transects to numbers
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

obs <- array(c(BM, TC), dim = c(1758, 4, 2))

M.BM <- Counts.BM

total.species <- vector()
for(i in 1:nrow(M.BM)){
  total.species[i] <- sum(M.BM[i,])
  M.BM$total.species[i] <- total.species[i]
}

M.BM <- M.BM$total.species

M.TC <- Counts.TC

total.species <- vector()
for(i in 1:nrow(M.TC)){
  total.species[i] <- sum(M.TC[i,])
  M.TC$total.species[i] <- total.species[i]
}

M.TC <- M.TC$total.species

M.obs <- cbind(M.BM, M.TC)

Counts.POV <- Counts %>%
  pivot_wider(id_cols = c(numeric_tran.grp), names_from = SPECIES, values_from = Count.POV,
              values_fn = sum, values_fill = 0)

Counts.POV <- Counts.POV[,c(3, 4, 5, 2)]
POV <- as.matrix(Counts.POV)

M.POV <- Counts.POV
total.species <- vector()
for(i in 1:nrow(M.POV)){
  total.species[i] <- sum(M.POV[i,])
  M.POV$total.species[i] <- total.species[i]
}

M <- M.POV$total.species

Counts.tr <- Counts %>% group_by(numeric_tran.grp) %>% slice(1L)
tr <- Counts.tr$numeric_tran

Counts.grp <- Counts %>% group_by(numeric_tran.grp) %>% slice(1L)
n.groups <- length(unique(Counts$numeric_tran.grp))
n.transects <- length(unique(Counts$numeric_tran))
n.species <- length(unique(Counts$SPECIES))
n.observers <- 2




ducks.mult <- "model{
###MIS-IDENTIFICATION MODEL 
#Priors and constraints

#Using Bayes theorem to get the probability that the species was K given that J was observed (sawJ.givK)  
#you can map transect to group so we model misID at the transect scale and detection at the group scale 

#priors on the confusion probabiities - multinomial logit
pr.sawJ.spK[1,1] <- 1-pr.sawJ.spK[2,1]-pr.sawJ.spK[3,1]-pr.sawJ.spK[4,1]
pr.sawJ.spK[2,1] <- exp(m.pr.sawJ.spK[2,1])/(1+exp(m.pr.sawJ.spK[2,1]) + exp(m.pr.sawJ.spK[3,1]) + 
                    exp(m.pr.sawJ.spK[4,1]))
pr.sawJ.spK[3,1] <- exp(m.pr.sawJ.spK[3,1])/(1+exp(m.pr.sawJ.spK[2,1]) + exp(m.pr.sawJ.spK[3,1]) +  
                    exp(m.pr.sawJ.spK[4,1]))
pr.sawJ.spK[4,1] <- exp(m.pr.sawJ.spK[4,1])/(1+exp(m.pr.sawJ.spK[2,1]) + exp(m.pr.sawJ.spK[3,1]) + 
                    exp(m.pr.sawJ.spK[4,1]))

pr.sawJ.spK[1,2] <- exp(m.pr.sawJ.spK[1,2])/(1+exp(m.pr.sawJ.spK[1,2]) + exp(m.pr.sawJ.spK[3,2]) + 
                    exp(m.pr.sawJ.spK[4,2]))
pr.sawJ.spK[2,2] <- 1-pr.sawJ.spK[1,2]-pr.sawJ.spK[3,2]-pr.sawJ.spK[4,2]
pr.sawJ.spK[3,2] <- exp(m.pr.sawJ.spK[3,2])/(1+exp(m.pr.sawJ.spK[1,2]) + exp(m.pr.sawJ.spK[3,2]) + 
                    exp(m.pr.sawJ.spK[4,2]))
pr.sawJ.spK[4,2] <- exp(m.pr.sawJ.spK[4,2])/(1+exp(m.pr.sawJ.spK[1,2]) + exp(m.pr.sawJ.spK[3,2]) + 
                    exp(m.pr.sawJ.spK[4,2]))

pr.sawJ.spK[1,3] <- exp(m.pr.sawJ.spK[1,3])/(1+exp(m.pr.sawJ.spK[1,3]) + exp(m.pr.sawJ.spK[2,3]) + 
                    exp(m.pr.sawJ.spK[4,3]))
pr.sawJ.spK[2,3] <- exp(m.pr.sawJ.spK[2,3])/(1+exp(m.pr.sawJ.spK[1,3]) + exp(m.pr.sawJ.spK[2,3]) + 
                    exp(m.pr.sawJ.spK[4,3]))
pr.sawJ.spK[3,3] <- 1-pr.sawJ.spK[1,3]-pr.sawJ.spK[2,3]-pr.sawJ.spK[4,3]
pr.sawJ.spK[4,3] <- exp(m.pr.sawJ.spK[4,3])/(1+exp(m.pr.sawJ.spK[1,3]) + exp(m.pr.sawJ.spK[2,3]) + 
                    exp(m.pr.sawJ.spK[4,3]))

pr.sawJ.spK[1,4] <- exp(m.pr.sawJ.spK[1,4])/(1+exp(m.pr.sawJ.spK[1,4]) + exp(m.pr.sawJ.spK[2,4]) + 
                    exp(m.pr.sawJ.spK[3,4]))
pr.sawJ.spK[2,4] <- exp(m.pr.sawJ.spK[2,4])/(1+exp(m.pr.sawJ.spK[1,4]) + exp(m.pr.sawJ.spK[2,4]) + 
                    exp(m.pr.sawJ.spK[3,4]))
pr.sawJ.spK[3,4] <- exp(m.pr.sawJ.spK[3,4])/(1+exp(m.pr.sawJ.spK[1,4]) + exp(m.pr.sawJ.spK[2,4]) + 
                    exp(m.pr.sawJ.spK[3,4]))
pr.sawJ.spK[4,4] <- 1-pr.sawJ.spK[1,4]-pr.sawJ.spK[2,4]-pr.sawJ.spK[3,4]

m.pr.sawJ.spK[2,1] ~ dnorm(0,0.01)      
m.pr.sawJ.spK[3,1] ~ dnorm(0,0.01)      
m.pr.sawJ.spK[4,1] ~ dnorm(0,0.01) 
m.pr.sawJ.spK[1,2] ~ dnorm(0,0.01)
m.pr.sawJ.spK[3,2] ~ dnorm(0,0.01)
m.pr.sawJ.spK[4,2] ~ dnorm(0,0.01)
m.pr.sawJ.spK[1,3] ~ dnorm(0,0.01)
m.pr.sawJ.spK[2,3] ~ dnorm(0,0.01)
m.pr.sawJ.spK[4,3] ~ dnorm(0,0.01)
m.pr.sawJ.spK[1,4] ~ dnorm(0,0.01)
m.pr.sawJ.spK[2,4] ~ dnorm(0,0.01)
m.pr.sawJ.spK[3,4] ~ dnorm(0,0.01)


#j is the species seen, k is the true species 
for(i in 1:n.groups){
  for(j in 1:n.species){  
  
    #first calculate the denominator for each species - Pr(saw J)
    #probability of seeing species j given species k * probabilities of species k on transect tr[i]
  
    pr.sawJ[tr[i],j] <- pr.sawJ.spK[j,1]*pr.spK[tr[i],1] + pr.sawJ.spK[j,2]*pr.spK[tr[i],2] + pr.sawJ.spK[j,3]*pr.spK[tr[i],3] + 
                        pr.sawJ.spK[j,4]*pr.spK[tr[i],4] 

    #then get the probability of the true species k given what you saw 
    #probability of seeing species j, given that it is species k * probability of species k on transect
    #divided by probability of seeing species j
  
    for(k in 1:n.species){  
      pr.spK.sawJ[tr[i],k,j] <- pr.sawJ.spK[j,k]*pr.spK[tr[i],k] / pr.sawJ[tr[i],j]
    }
  }
}

#prior on species composition for transect t
for(t in 1:n.transects){
  pr.spK[t,] ~ ddirch(alpha[])
}
alpha <- rep(1, 4)

#this is modeled at the transect level, tr[group i]
for (i in 1:n.groups){
  for(k in 1:n.species){
    for(o in 1:n.observers){
      for(j in 1:n.species){
    
        psi[tr[i],j,k,o] <- pr.sawJ.spK[j,k]*pr.spK[tr[i],k]                  
      }

        #observer data - multinomial sample from total number of birds in the group
        #probabilities are pr(see J) but we estimate these separately for every true species k 
        ID[i,,o] ~ dmultinom(psi[tr[i],,k,o], M.obs[i,o]) 
    }
  }
}


###DETECTION AND AVAILABILITY MODELS  
###Everything here in terms of true k except original observation data 
#Need to finalize data input to be the same for each section of the model (detection/availability and ID)
for(i in 1:n.groups){
  for(k in 1:n.species){


    ## Process model
    #POV camera data - multinomial sample from the total number of birds in the group, with pr.spK probabilities   
    
    #I am not totally sure if it is going to work to model pr.spK like this here
    #if not, we can just take the mean over transects elsewhere
    #and use the mean over transects in the misID model 
    POV[i,j] ~ dmultinom(pr.spK[tr[i],j], M[i])

    #FF camera data - abundance of species j in group i prior to plane passing (data), modeled with a mean mu
    FF[i,j] ~ dpois(mu[i,j])
    
    # prior on abundance of species j in group i
    lambda[i,j] ~ dgamma(0.1,0.1) 

    # abundance of species j in group i after plane passing (predicted from N-mix), modeled with a mean lambda
    #this is where we might need the negbin with overdispersion parameter for species(group) size distribution
    N[i,j] ~ dpois(lambda[i,j]) 
    
    #the relationship between abundance before and after plane, by species
    mu[i,j] <- beta[j] * lambda[i,j]
    

 #use dsum to sum Ntrue.givenseen over true individuals  
    for(o in 1:n.observers){
      
      #translate observation data to what was seen 
      true[i,k,o] ~ dmultinom(pr.spK.sawJ[tr[i],,j], obs[i,j,o])
      
      #Ntrue.givenseen is a multinomial outcome from the true abundance and the pr that you saw J given K
      #I feel like I might have some logic error in here but my brain is fried
      #weird that we aren't using the outcome of Bayes calculation 
      #Nseen.true[i,j,,o] ~ dmultinom(pr.sawJ.spK[,j],true[i,j,o]) 
      #obs[i,j,o] ~ dsum(Nseen.true[i,j,,o]) 

      
      #true[i,k,o] ~ dbin(p[k,o],N[i,k])
    }  
  }     
}

for(j in 1:n.species){
  #prior on relationship between abundance before/after - allows movement in both directions, but we expect beta > 1
  beta[j] ~ dgamma(0.1,0.1)

  # prior on detection probability of species j by observer o
  for(o in 1:n.observers){
    p[j,o] ~ dbeta(1, 1)
  }
}  
  
}"


data <- list(FF = FF, obs = obs, ID = ID, POV = POV, M = M, M.obs = M.obs, tr = tr, n.transects = n.transects, 
             n.groups = n.groups, n.species = n.species, n.observers = n.observers)

#create smart initial values for beta 
inits.beta <- function(FF,obs){
  max.FF <- FF
  max.FF[which(FF == 0)] <- 1
  max.obs <- apply(obs,1,max)
  max.obs[which(max.obs == 0)] <- 1
  inits.beta <- max.FF/max.obs
  return(inits.beta)
}  

inits.mult<-function(){list(p = runif(2, 0, 1), lambda = apply(obs, 1, max), 
                            beta = runif(1,0,1))}

params.mult <- c("p", "beta")

#Generate samples from the posterior distribution
mult.mod<-run.jags(model = ducks.mult,
                   monitor = params.mult,
                   data = data,
                   n.chains = 3,
                   burnin = 500,
                   sample = 1000,
                   inits = inits.mult)






