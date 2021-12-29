detection.exp <- "model{


###THIS IS THE MIS-ID MODEL 
#priors and constraints 

#priors on the confusion probabiities - mlogit
pr.sawJ.givK[1,1] <- 1-pr.sawJ.givK[2,1]-pr.sawJ.givK[3,1]-pr.sawJ.givK[4,1]
pr.sawJ.givK[2,1] <- exp(m.pr.sawJ.givK[2,1])/(1+exp(m.pr.sawJ.givK[2,1]) + exp(m.pr.sawJ.givK[3,1]) + exp(m.pr.sawJ.givK[4,1]))
pr.sawJ.givK[3,1] <- exp(m.pr.sawJ.givK[3,1])/(1+exp(m.pr.sawJ.givK[2,1]) + exp(m.pr.sawJ.givK[3,1]) + exp(m.pr.sawJ.givK[4,1]))
pr.sawJ.givK[4,1] <- exp(m.pr.sawJ.givK[4,1])/(1+exp(m.pr.sawJ.givK[2,1]) + exp(m.pr.sawJ.givK[3,1]) + exp(m.pr.sawJ.givK[4,1]))

pr.sawJ.givK[1,2] <- exp(m.pr.sawJ.givK[1,2])/(1+exp(m.pr.sawJ.givK[1,2]) + exp(m.pr.sawJ.givK[3,2]) + exp(m.pr.sawJ.givK[4,2]))
pr.sawJ.givK[2,2] <- 1-pr.sawJ.givK[1,2]-pr.sawJ.givK[3,2]-pr.sawJ.givK[4,2]
pr.sawJ.givK[3,2] <- exp(m.pr.sawJ.givK[3,2])/(1+exp(m.pr.sawJ.givK[1,2]) + exp(m.pr.sawJ.givK[3,2]) + exp(m.pr.sawJ.givK[4,2]))
pr.sawJ.givK[4,2] <- exp(m.pr.sawJ.givK[4,2])/(1+exp(m.pr.sawJ.givK[1,2]) + exp(m.pr.sawJ.givK[3,2]) + exp(m.pr.sawJ.givK[4,2]))

pr.sawJ.givK[1,3] <- exp(m.pr.sawJ.givK[1,3])/(1+exp(m.pr.sawJ.givK[1,3]) + exp(m.pr.sawJ.givK[2,3]) + exp(m.pr.sawJ.givK[4,3]))
pr.sawJ.givK[2,3] <- exp(m.pr.sawJ.givK[2,3])/(1+exp(m.pr.sawJ.givK[1,3]) + exp(m.pr.sawJ.givK[2,3]) + exp(m.pr.sawJ.givK[4,3]))
pr.sawJ.givK[3,3] <- 1-pr.sawJ.givK[1,3]-pr.sawJ.givK[2,3]-pr.sawJ.givK[4,3]
pr.sawJ.givK[4,3] <- exp(m.pr.sawJ.givK[4,3])/(1+exp(m.pr.sawJ.givK[1,3]) + exp(m.pr.sawJ.givK[2,3]) + exp(m.pr.sawJ.givK[4,3]))

pr.sawJ.givK[1,4] <- exp(m.pr.sawJ.givK[1,4])/(1+exp(m.pr.sawJ.givK[1,4]) + exp(m.pr.sawJ.givK[2,4]) + exp(m.pr.sawJ.givK[3,4]))
pr.sawJ.givK[2,4] <- exp(m.pr.sawJ.givK[2,4])/(1+exp(m.pr.sawJ.givK[1,4]) + exp(m.pr.sawJ.givK[2,4]) + exp(m.pr.sawJ.givK[3,4]))
pr.sawJ.givK[3,4] <- exp(m.pr.sawJ.givK[3,4])/(1+exp(m.pr.sawJ.givK[1,4]) + exp(m.pr.sawJ.givK[2,4]) + exp(m.pr.sawJ.givK[3,4]))
pr.sawJ.givK[4,4] <- 1-pr.sawJ.givK[1,4]-pr.sawJ.givK[2,4]-pr.sawJ.givK[3,4]

log(m.pr.sawJ.givK[2,1]) ~ dnorm(0,0.01)      
log(m.pr.sawJ.givK[3,1]) ~ dnorm(0,0.01)      
log(m.pr.sawJ.givK[4,1]) ~ dnorm(0,0.01) 
log(m.pr.sawJ.givK[1,2]) ~ dnorm(0,0.01)
log(m.pr.sawJ.givK[3,2]) ~ dnorm(0,0.01)
log(m.pr.sawJ.givK[4,2]) ~ dnorm(0,0.01)
log(m.pr.sawJ.givK[1,3]) ~ dnorm(0,0.01)
log(m.pr.sawJ.givK[2,3]) ~ dnorm(0,0.01)
log(m.pr.sawJ.givK[4,3]) ~ dnorm(0,0.01)
log(m.pr.sawJ.givK[1,4]) ~ dnorm(0,0.01)
log(m.pr.sawJ.givK[2,4]) ~ dnorm(0,0.01)
log(m.pr.sawJ.givK[3,4]) ~ dnorm(0,0.01)


#this is using Bayes theorem to get what you need: the probability the species was K given you saw J  
#you can map transect to group so we model misID at the transect scale and detection at the group scale 

#j is the putative species
for(j in 1:n.species){  
  #first calculate the denominator for each species - Pr(saw J)
  pr.sawJ[j] <- pr.sawJ.givK[j,1]*pr.s[tr[i],1] + pr.sawJ.givK[j,2]*pr.s[tr[i],2] + pr.sawJ.givK[j,3]*pr.s[tr[i],3] + pr.sawJ.givK[j,4]*pr.s[tr[i],4] 

  #then get the probability of the true species given what you saw 
  #k is the true species
  for(k in 1:n.species){  
    pr.K.givsawJ[k,j] <- pr.sawJ.givK[j,k]*pr.s[tr[i],k] / pr.sawJ[j]
  }
}

#prior on species composition in group i 
for(i in 1:n.groups){
  pr.s[i,] ~ ddirch(alpha[])
}
alpha <- c(1,1,1,1)


#this is modeled at the transect level 
for (tr in 1:n.transects){
  for(j in 1:n.species){
    for(o in 1:n.observers){
    
      #matrix of confusion probabilities 
      
      #i haven't figured out the indexing yet for group to transect 
      psi[tr,1,1,o] <- pr.s[tr,1]*pr.sawJ.givK[1,1]  #Pr(saw 1|1)                
      psi[tr,1,2,o] <- pr.s[tr,1]*pr.sawJ.givK[2,1]  #Pr(saw 1|2)   
      psi[tr,1,3,o] <- pr.s[tr,1]*pr.sawJ.givK[3,1]  #Pr(saw 1|3)
      psi[tr,1,4,o] <- pr.s[tr,1]*pr.sawJ.givK[4,1]  #Pr(saw 1|4)  
      psi[tr,2,1,o] <- pr.s[tr,2]*pr.sawJ.givK[1,2]  #Pr(saw 2|1)
      psi[tr,2,2,o] <- pr.s[tr,2]*pr.sawJ.givK[2,2]  #Pr(saw 2|2)    
      psi[tr,2,3,o] <- pr.s[tr,2]*pr.sawJ.givK[3,2]  #Pr(saw 2|3)  
      psi[tr,2,4,o] <- pr.s[tr,2]*pr.sawJ.givK[4,2]  #Pr(saw 2|4)  
      psi[tr,3,1,o] <- pr.s[tr,3]*pr.sawJ.givK[1,3]  #Pr(saw 3|1)
      psi[tr,3,2,o] <- pr.s[tr,3]*pr.sawJ.givK[2,3]  #Pr(saw 3|2)
      psi[tr,3,3,o] <- pr.s[tr,3]*pr.sawJ.givK[3,3]  #Pr(saw 3|3)
      psi[tr,3,4,o] <- pr.s[tr,3]*pr.sawJ.givK[4,3]  #Pr(saw 3|4)  
      psi[tr,4,1,o] <- pr.s[tr,4]*pr.sawJ.givK[1,4]  #Pr(saw 4|1)
      psi[tr,4,2,o] <- pr.s[tr,4]*pr.sawJ.givK[2,4]  #Pr(saw 4|2)
      psi[tr,4,3,o] <- pr.s[tr,4]*pr.sawJ.givK[3,4]  #Pr(saw 4|3)
      psi[tr,4,4,o] <- pr.s[tr,4]*pr.sawJ.givK[4,4]  #Pr(saw 4|4)

      #observer data - multinomial sample from total number of birds in the group and the confusion probabilities 
      ID[t,j,o] ~ dmultinom(psi[i,j,,o], M.obs[i,o])
    }
  }
}
    

###THESE ARE THE DETECTION AND AVAILABILITY MODELS  
for(i in 1:n.groups){
  for(j in 1:n.species){

    # prior on abundance of species j in group i
    lambda[i,j] ~ dgamma(0.1,0.1) 

    ## Process model
    #POV camera data - multinomial sample from the total number of birds in the group, with pr.s probabilities
    
    #I am not totally sure if it is going to work to model pr.s like this here
    #if not, we can just take the mean over transects elsewhere
    #and use the mean over transects in the misID model 
    POV[i,j] ~ dmultinom(pr.s[tr[i],j], M[i])

    #FF camera data - abundance of species j in group i prior to plane passing (data), modeled with a mean mu
    FF[i,j] ~ dpois(mu[i,j]) 

    #abundance of species j in group i after plane passing (predicted from N-mix), modeled with a mean lambda
    #this is where we might need the negbin with overdispersion parameter for species(group) size distribution
    N[i,j] ~ dpois(lambda[i,j]) 

    #the relationship between abundance before and after plane, by species
    mu[i,j] <- beta[j] * lambda[i,j]
  

    #use dsum to sum Ntrue.givenseen over true individuals  
    for(o in 1:n.observer){
      obs[i,j,o] ~ dsum(Nseen.true[i,j,,o]) 
      
      #Ntrue.givenseen is a multinomial outcome from the true abundance and the pr that you saw J given K
      #I feel like I might have some logic error in here but my brain is fried
      #weird that we aren't using the outcome of Bayes calculation 
      Nseen.true[i,j,,o] ~ dmultinom(pr.sawJ.givK[,j],true[i,j,o]) 

      #
      true[i,j,o] ~ dbin(p[j,o],N[i,j])
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

#FF is a matrix of dimension groups, species - this is the count of species j in group i by the FF camera 
#obs is an array of dimension groups, species, observers - this is the count of species j in group i by observer o
#ID is an array of dimension transects, species, observers - this is the sum of obs over groups to the transect level 
#POV is a matrix of dimension groups, species - this is the count of species j in group i by the POV camera
#M is a vector of dimension groups - this is the total count across species in group i by the POV camera
#tr is a vector of dimension groups that maps group to transect number (for group i, what is the transect number t)
#scalars: n.observers, n.species, n.transects, n.groups
data <- list(FF = FF, obs = obs, ID = ID, POV = POV, M = M, tr = tr, n.transects = n.transects, n.groups = n.groups, n.species = n.species, n.observers = n.observers)
  
