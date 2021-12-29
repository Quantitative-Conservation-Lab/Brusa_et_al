
#this is the range of your abundances 
#currently very similar but you could make it really different and you'd get 0s
UL <- 20 #lowest abundance/transect of most rare species
LL <- 30 #highest abundance/transect of most common species 

#number of species in the system 
sp <- 4
#number of transects 
transects <- 100 

#create a confusion matrix (dimensions = sp * sp)
conf.levels <- c(1,0.95,0.9,0.85,0.8,0.75,0.7)
confusion <- matrix(NA,nrow=sp,ncol=sp)
for(i in 1:sp){
  conf <- sample(conf.levels,1)
  for(j in 1:sp){
    tot.prob <- 1-conf 
    if(i == j){
    confusion[i,j] <- conf
    }else confusion[i,j] <- tot.prob/(sp-1)
  }
}

#we simulate true abundances 
#then apply the confusion matrix 
#to get the sightings by transect, true species, and identified species
abund <- matrix(NA,nrow=transects,ncol=sp)
sightings <- array(NA,dim = c(transects,sp,sp))
for(i in 1:transects){
  lambda <- sample(c(LL:UL),sp)
  for(j in 1:sp){
    abund[i,j] <- rpois(1,lambda[j])
    sightings[i,j,] <- rmultinom(1,abund[i,j],confusion[j,])
  }
}
#we then sum over true species to get identified species for each transect 
data <- apply(sightings,c(1,3),sum)

#calculate the table from Aebischer
out <- array(NA,dim = c(transects,sp,sp))
for(i in 1:transects){
  for(j in 1:sp){
    for(k in 1:sp){
      out[i,j,k] <- log(data[i,k]/data[i,j]) #- log(abund[i,k]/abund[i,j])
    }
  }
}
#this is the mean over transects of the table from Aebisher 
#what is the relationship between this and confusion? 
apply(out,c(2,3),mean)

