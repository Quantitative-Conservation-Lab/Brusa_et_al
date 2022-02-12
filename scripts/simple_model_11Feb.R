library(tidyverse)
library(stringr)
library(jagsUI)
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

M.ID <- apply(ID,c(1,3),sum)

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

#Counts.BM <- Counts %>%
#  pivot_wider(id_cols = c(numeric_tran.grp), names_from = SPECIES, values_from = Count.BM,
#              values_fn = sum, values_fill = 0)

#Counts.BM <- Counts.BM[,c(3, 4, 5, 2)]
#BM <- as.matrix(Counts.BM)

#Counts.TC <- Counts %>%
#  pivot_wider(id_cols = c(numeric_tran.grp), names_from = SPECIES, values_from = Count.TC,
#              values_fn = sum, values_fill = 0)

#Counts.TC <- Counts.TC[,c(3, 4, 5, 2)]
#TC <- as.matrix(Counts.TC)

#obs <- array(c(BM, TC), dim = c(1758, 4, 2))

#M.BM <- Counts.BM

#total.species <- vector()
#for(i in 1:nrow(M.BM)){
#  total.species[i] <- sum(M.BM[i,])
#  M.BM$total.species[i] <- total.species[i]
#}

#M.BM <- M.BM$total.species

#M.TC <- Counts.TC

#total.species <- vector()
#for(i in 1:nrow(M.TC)){
#  total.species[i] <- sum(M.TC[i,])
#  M.TC$total.species[i] <- total.species[i]
#}

#M.TC <- M.TC$total.species

#M.obs <- cbind(M.BM, M.TC)

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

#choose samples to break up data.set 
mis.nos <- sample(c(1:n.groups),round(n.groups/2,dig=0))
mis.nos <- sort(mis.nos)
det.nos <- setdiff(c(1:n.groups),mis.nos)

#Clean up the data processing script
#Want to make sure we're getting the right stuff 

#We need obs, POV, and FF data, plus tr object - seems like there is more than that above? 
#what is obs vs ID calculated above - they seem to be the same 

ID.mis <- ID[mis.nos,,]
ID.det <- ID[det.nos,,]
M.obs.mis <- M.ID[mis.nos,] #only need this for misID 
POV.mis <- POV[mis.nos,] #only need this for misID
M.mis <- M[mis.nos] #only need this for misID 
FF.det <- FF[det.nos,] #only need this for detection
tr.mis <- as.numeric(as.factor(tr[mis.nos])) #only need this for misID 
n.transects.mis <- length(unique(tr.mis)) #only need this for misID
n.groups.mis <- dim(ID.mis)[1]
n.groups.det <- dim(ID.det)[1]
n.species <- dim(ID.mis)[2]
n.observers <- dim(ID.mis)[3]

cat("
model {

#########################SPECIES MIS-ID COMPONENT#########################

#modeled with data cut in half (n.transects.mis and n.groups.mis vs n.transects.det and n.groups.det)

#species composition in observer data with random effects of transect 
for(t in 1:n.transects.mis){
  for(o in 1:n.observers){
    pr.ID[t,o,1] <- exp(m.pr.ID[t,o,1])/(1 + exp(m.pr.ID[t,o,1]) + exp(m.pr.ID[t,o,2]) + exp(m.pr.ID[t,o,3]))
    pr.ID[t,o,2] <- exp(m.pr.ID[t,o,2])/(1 + exp(m.pr.ID[t,o,1]) + exp(m.pr.ID[t,o,2]) + exp(m.pr.ID[t,o,3]))
    pr.ID[t,o,3] <- exp(m.pr.ID[t,o,3])/(1 + exp(m.pr.ID[t,o,1]) + exp(m.pr.ID[t,o,2]) + exp(m.pr.ID[t,o,3]))
    pr.ID[t,o,4] <- 1-pr.ID[t,o,1]-pr.ID[t,o,2]-pr.ID[t,o,3]

    m.pr.ID[t,o,1] ~ dnorm(mn.pr.ID[1],tau.ID)      
    m.pr.ID[t,o,2] ~ dnorm(mn.pr.ID[2],tau.ID)      
    m.pr.ID[t,o,3] ~ dnorm(mn.pr.ID[3],tau.ID)      
  }
}
#overall multinomial logit means 
for(s in 1:(n.species-1)){
  mn.pr.ID[s] ~ dnorm(0,0.001)
}
#overall random effect SDs 
tau.ID <- pow(sigma.ID,-2)
sigma.ID ~ dunif(0,5)

#the species ID data by group and observer are multinomial given the proportion of species J in the observations and the total observed      
for(g in 1:n.groups.mis){
  for(o in 1:n.observers){
    ID.mis[g,1:4,o] ~ dmultinom(pr.ID[tr.mis[g],o,],M.obs.mis[g,o])
  }
}

#back-predict the overall probability of species J in the observer data 
for(s in 1:(n.species-1)){
  pr.ID.pred[s] <- exp(mn.pr.ID[s])/(1+exp(mn.pr.ID[1])+exp(mn.pr.ID[2])+exp(mn.pr.ID[3]))
}
pr.ID.pred[4] <- 1-pr.ID.pred[1]-pr.ID.pred[2]-pr.ID.pred[3]


#########################SPECIES COMPOSITION COMPONENT(POV)#########################

#species composition in POV data with random effects of transect 
for(t in 1:n.transects.mis){
  pr.sp[t,1] <- exp(m.pr.sp[t,1])/(1 + exp(m.pr.sp[t,1]) + exp(m.pr.sp[t,2]) + exp(m.pr.sp[t,3]))
  pr.sp[t,2] <- exp(m.pr.sp[t,2])/(1 + exp(m.pr.sp[t,1]) + exp(m.pr.sp[t,2]) + exp(m.pr.sp[t,3]))
  pr.sp[t,3] <- exp(m.pr.sp[t,3])/(1 + exp(m.pr.sp[t,1]) + exp(m.pr.sp[t,2]) + exp(m.pr.sp[t,3]))
  pr.sp[t,4] <- 1-pr.sp[t,1]-pr.sp[t,2]-pr.sp[t,3]

  m.pr.sp[t,1] ~ dnorm(mn.pr.sp[1],tau.sp)      
  m.pr.sp[t,2] ~ dnorm(mn.pr.sp[2],tau.sp)      
  m.pr.sp[t,3] ~ dnorm(mn.pr.sp[3],tau.sp)      
}

for(s in 1:(n.species-1)){
  mn.pr.sp[s] ~ dnorm(0,0.001)
}
tau.sp <- pow(sigma.sp,-2)
sigma.sp ~ dunif(0,5)

#POV camera data - multinomial sample from the total number of birds in the group, with pr.spK probabilities   
for(i in 1:n.groups.mis){
  POV.mis[i,1:4] ~ dmultinom(pr.sp[tr.mis[i],], M.mis[i])
}

for(i in 1:(n.species-1)){
  pr.sp.pred[i] <- exp(mn.pr.sp[i])/(1+exp(mn.pr.sp[1])+exp(mn.pr.sp[2])+exp(mn.pr.sp[3]))
}
pr.sp.pred[4] <- 1-pr.sp.pred[1]-pr.sp.pred[2]-pr.sp.pred[3]


#########################MIS-ID COMPONENT#########################

#Derived mis-ID ratios 
#ratio of species J in the observer data to species J in the POV data 
for(i in 1:n.species){
  mis.ID[i] <- pr.ID.pred[i]/pr.sp.pred[i]
}

#Correction abundance estimates from the other half of data for mis-ID 
for(i in 1:n.groups.det){
  for(k in 1:n.species){
    N.spK[i,k] <- N.obsJ[i,k]/mis.ID[k]
  }
}

#sum both corrected and uncorrected estimates across groups 
for(k in 1:n.species){
  N.spKall[k] <- sum(N.spK[,k])
  N.obsJall[k] <- sum(N.obsJ[,k])
}

#########################DETECTION COMPONENT#########################
  
#prior for abundance of species k in group i after plane passing (predicted from N-mix)
#probably need a negative binomial here 
for(i in 1:n.groups.det){
  for(k in 1:n.species){
    N.obsJ[i,k] ~ dpois(lambda[i,k]) 
    
  }
}

#estimation of the abundance of putative species J - Nmixture
for(j in 1:n.species){
  for(o in 1:n.observers){
    for(i in 1:n.groups.det){
      ID.det[i,j,o] ~ dbinom(p[j,o],N.obsJ[i,j])
    }
    #prior for detection probability 
    p[j,o] ~ dbeta(1,1)
  }
}  

#########################AVAILABILITY COMPONENT(FF vs OBSERVERS)#########################

for(i in 1:n.groups.det){
  for(k in 1:n.species){
    #FF camera data - abundance of species k in group i prior to plane passing (data), modeled with a mean mu
    #might need a negative binomial here 
    FF.det[i,k] ~ dpois(mu[i,k])
    
    #prior on abundance of species j in group i from FF#
    mu[i,k] ~ dgamma(1,1) 
    
    #the relationship between abundance from observers and abundance before the plane passed, by species
    lambda[i,k] <- beta[k] * mu[i,k] 
  }
}

#prior on beta 
for(k in 1:n.species){
  beta[k] ~ dunif(0,5)
}  

##########################################################################################

}",file = "ducks.txt")

#DATA
data.duck <- list(ID.mis=ID.mis,M.obs.mis=M.obs.mis,ID.det=ID.det,POV.mis=POV.mis,M.mis=M.mis,FF.det=FF.det,tr.mis=tr.mis,n.transects.mis=n.transects.mis,n.groups.mis=n.groups.mis,n.groups.det=n.groups.det,n.observers=n.observers,n.species=n.species)

inits.duck<-function(){list(N.obsJ=(apply(ID.det,c(1,2),sum)+matrix(rep(1,879*4),nrow=879,ncol=4)))}

#PARAMETERS
params.duck <- c("pr.ID.pred","pr.sp.pred","mis.ID")

#Generate samples from the posterior distribution
out.jags = jags(data.duck, inits.duck, params.duck, model.file="ducks.txt",
                n.chains=3, n.iter=20000, n.burnin=5000, n.thin=1)

