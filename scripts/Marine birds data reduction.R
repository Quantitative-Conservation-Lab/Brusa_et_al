library(tidyverse)
library(nimble)
library(ggmcmc)
library(coda)
library(MCMCvis)
library(viridis)
library(here)

s.allo_df <- read.csv(here('Data', 's_allo.csv'))


##Prep data for model##

#The specific model data set up and model are currently set up to run one model with only observer records from the
#observer sitting in the middle (also called front) seat and one model with only observer records from the observer
#sitting in the rear seat. The general model is set up for seat assignment as a covariate, but it is not functional as is

#Set up array for observer records
s.allo_df$Transect <- str_trunc(s.allo_df$TransectGroup, 12, "right", ellipsis = "")
ID.obs <- s.allo_df

tran <- unique(ID.obs$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))

ID.obs$numeric_tran <- 0
for(i in dumb$tran){
  ID.obs$numeric_tran[ID.obs$Transect == i] <- match(i, dumb$tran)
}

#Preserve group sizes since they are changed with the reallocation of unknown species
grps <- unique(ID.obs$TransectGroup)

grp.sizes <- ID.obs %>% group_by(TransectGroup) %>%
  summarise(ff.grp.sz = sum(Count.FF))

num.gr <- data.frame(TransectGroup = grps,
                     grp.sizes = grp.sizes$ff.grp.sz)

ID.obs <- ID.obs %>% left_join(x = ID.obs, y = num.gr, by = "TransectGroup")

ID.obs$SPECIES[str_sub(ID.obs$SPECIES,start = 1, end = 1) == "U"] <- "Z.UNK"
ID.obs$sp.group[str_sub(ID.obs$sp.group,start = 1, end = 1) == "U"] <- "Z.UNK"

#Comment/uncomment below to set up model for middle seat observations only or rear seat observations only
#ID.obs <- ID.obs %>% filter(Fr.Det == 1)
ID.obs <- ID.obs %>% filter(Rr.Det == 1)


ID.BM <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = sp.group, values_from = Count.BM,
              values_fn = sum, values_fill = 0)


ID.BM <- as.matrix(ID.BM[,-1])

ID.BM <- ID.BM[,order(colnames(ID.BM))]


ID.TC <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = sp.group, values_from = Count.TC,
              values_fn = sum, values_fill = 0)


ID.TC <- as.matrix(ID.TC[,-1])

ID.TC <- ID.TC[,order(colnames(ID.TC))]

OBS <- array(c(ID.BM, ID.TC), dim = c(nrow(ID.BM), ncol(ID.BM), 2))

BM.total <- apply(ID.BM, 1, sum)

TC.total <- apply(ID.TC, 1, sum)

OBS.total <- cbind(BM.total, TC.total)


ID.POV <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = sp.group, values_from = Count.POV,
              values_fn = sum, values_fill = 0)


ID.POV <- as.matrix(ID.POV[,-1])

ID.POV <- ID.POV[,order(colnames(ID.POV))]

POV <- ID.POV

POV.total <- apply(POV, 1, sum)

ID.FF <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = sp.group, values_from = Count.FF,
              values_fn = sum, values_fill = 0)


ID.FF <- as.matrix(ID.FF[,-1])

ID.FF <- ID.FF[,order(colnames(ID.FF))]

FF <- ID.FF

FF.total <- apply(FF, 1, sum)

# #Set up seat position covariate
# ID.obs$Obs1.Front <- ifelse(ID.obs$Front == "BM", 1, 0)
# Obs1.Front <- ID.obs %>% group_by(Transect) %>% summarise(Obs1.Front = round(mean(Obs1.Front)))
# Obs1.Front <- Obs1.Front$Obs1.Front
# 
# ID.obs$Obs2.Front <- ifelse(ID.obs$Front == "TC", 1, 0)
# Obs2.Front <- ID.obs %>% group_by(Transect) %>% summarise(Obs2.Front = round(mean(Obs2.Front)))
# Obs2.Front <- Obs2.Front$Obs2.Front


ID.obs <- ID.obs[order(ID.obs$sp.group),]
nspecies <- length(unique(ID.obs$sp.group)) -1
nsites <- length(unique(ID.obs$numeric_tran))
nobs <- 2


##Less detailed model##
#-Nimble Code-#

code <- nimbleCode({
  
  #-Priors-#
  
  #Composition of latent abundance (corrected for imperfect detection)

  for(o in 1:nobs){
    phi[1:(nspecies+1), o] ~ ddirch(phi.ones[1:(nspecies+1),o])
    
    for(i in 1:(nspecies+1)){
      phi.ones[i,o] <- 1
    }

  }
  
  #Derived product of movement and detection
  for(o in 1:nobs){
    E.epsilon[o] ~ dnorm(0, 0.01)
  }
  
  #-Likelihood-#
  
  for(j in 1:nsites){
    
    #Front facing camera composition
    FF[j,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[j])
    
    #Front facing camera total abundance
    FF.total[j] ~ dpois(lambda.total)
    
    for(o in 1:nobs){
      
      OBS[j,1:(nspecies+1),o] ~ dmulti(phi[1:(nspecies+1),o], OBS.total[j,o])
      
      OBS.total[j,o] ~ dpois(lambda.total * E.epsilon[o])
      
    }#end o
    
  }#end j
  
  for(i in 1:nspecies){
    
    pi[i] <- lambda[i]/lambda.total
    
   }

     
  for(i in 1:nspecies){
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    

    for(o in 1:nobs){
     
      
    correction[i,o] <- E.epsilon[o] * phi[i,o]/pi[i] 
    
    }#end o
  }#end i
  
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Compile data-#

data <- list(FF = FF,
             FF.total = FF.total,
             # Obs1.Front = Obs1.Front,
             # Obs2.Front = Obs2.Front,
             OBS = OBS,
             OBS.total = OBS.total
)


constants <- list(nspecies = nspecies, nsites = nsites, nobs = 2)


#-Initial values-#

inits <- function(){list(pi = apply(FF/apply(FF, 1, sum), 2, mean, na.rm = TRUE),
                         E.epsilon = rep(sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))
                                         *runif(1, 0.25, 1)), 2),
                         # int.epsilon = rep(sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))), 2),
                         # beta.Obs1.Front <- rnorm(1, 0, 0.1),
                         # beta.Obs2.Front <- rnorm(1, 0, 0.1),
                         phi = apply(apply(OBS, c(1,2,3), max)/apply(apply(OBS, c(1,2,3), max), 1, sum), c(2,3), mean,
                                     na.rm = TRUE)
)}

#-Parameters to save-#

params <- c(
  "pi",
  "phi",
  # "beta.Obs1.Front",
  # "beta.Obs2.Front",
  "lambda",
  "lambda.total",
  "E.epsilon",
  "correction"
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

s.out <- runMCMC(compiled.model$MCMC,
               niter = ni, nburnin = nb,
               nchains = nc, thin = nt,
               samplesAsCodaMCMC = TRUE)


#Diagnostics
s.out.mcmc <- as.mcmc.list(s.out)
s.out.ggs <- ggs(s.out.mcmc)

s.out.diag <- ggs_diagnostics(s.out.ggs)
s.out.gelman <- gelman.diag(s.out.mcmc)
MCMCtrace(s.out.mcmc, params = c("pi",
                               "phi",
                               # "beta.Obs1.Front",
                               # "beta.Obs2.Front",
                               "lambda",
                               "lambda.total",
                               "E.epsilon",
                               "correction"
                               ))

#Work with output
#Check dimensions - 20K iterations minus 10K for burn-in, 194 parameters
setwd("~/Documents/Windsor/UW Postdoc/Sea duck detection")
s.out <- readRDS("specificfront_model.RDS") 
#OR
s.out <- readRDS("specificrear_model.RDS")
dim(s.out[[1]])
s.mcmc.params <- as.mcmc.list(s.out)
s.params.groups <- data.frame(as.matrix(s.mcmc.params)) 

s.params.groups = as.matrix(s.params.groups)

s.out.groups <- data.frame(
  Mean = apply(s.params.groups, 2, mean),
  lcl = apply(s.params.groups, 2, quantile, probs = c(.05), na.rm = TRUE),
  ucl = apply(s.params.groups, 2, quantile, probs = c(.95), na.rm = TRUE),
  SD = apply(s.params.groups, 2, sd))

s.out.g.summ <- MCMCsummary(s.mcmc.params)

FF.s.sums <- apply(FF, 2, sum)

Species <- colnames(FF)

lambdas.s <- s.out.groups[c(11:18),]
lambda.mean.s <- nsites*lambdas.s$Mean

epsilons.s <- s.out.groups[c(3:10),1]
comp.ratio.s <- s.out.groups[c(36:51),1] / s.out.groups[c(60:67),1] #phi/psi

birds.adjust.s <- lambda.mean.s*(1/(comp.ratio.s*epsilons.s))

s.bird.groups <- as.data.frame(cbind(Species[1:8], FF.sums[1:8], round(lambda.mean.s, 2), round(birds.adjust.s, 2)))

colnames(s.bird.groups)[1:4] <- c("Species group","Counts from FF","Unadjusted N estimate", "Adjusted N estimate")
row.names(s.bird.groups) <- NULL



######
#For general model

#This code is currently set up to try to incorporate seat assignment as a covariate...but needs work

g.allo_df <- read.csv(here('Data', 'g_allo.csv'))


##Prep data for model##

#Set up array for observer records
g.allo_df$Transect <- str_trunc(g.allo_df$TransectGroup, 12, "right", ellipsis = "")
ID.obs <- g.allo_df

tran <- unique(ID.obs$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))

ID.obs$numeric_tran <- 0
for(i in dumb$tran){
  ID.obs$numeric_tran[ID.obs$Transect == i] <- match(i, dumb$tran)
}

#Preserve group sizes since they are changed with the reallocation of unknown species
grps <- unique(ID.obs$TransectGroup)

grp.sizes <- ID.obs %>% group_by(TransectGroup) %>%
  summarise(ff.grp.sz = sum(Count.FF))

num.gr <- data.frame(TransectGroup = grps,
                     grp.sizes = grp.sizes$ff.grp.sz)

ID.obs <- ID.obs %>% left_join(x = ID.obs, y = num.gr, by = "TransectGroup")

ID.obs$SPECIES[str_sub(ID.obs$SPECIES,start = 1, end = 1) == "U"] <- "Z.UNK"
ID.obs$group[str_sub(ID.obs$group,start = 1, end = 1) == "U"] <- "Z.UNK"


ID.BM <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = group, values_from = Count.BM,
              values_fn = sum, values_fill = 0)


ID.BM <- as.matrix(ID.BM[,-1])

ID.BM <- ID.BM[,order(colnames(ID.BM))]

ID.TC <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = group, values_from = Count.TC,
              values_fn = sum, values_fill = 0)


ID.TC <- as.matrix(ID.TC[,-1])

ID.TC <- ID.TC[,order(colnames(ID.TC))]

OBS <- array(c(ID.BM, ID.TC), dim = c(nrow(ID.BM), ncol(ID.BM), 2))

BM.total <- apply(ID.BM, 1, sum)

TC.total <- apply(ID.TC, 1, sum)

OBS.total <- cbind(BM.total, TC.total)


ID.POV <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = group, values_from = Count.POV,
              values_fn = sum, values_fill = 0)


ID.POV <- as.matrix(ID.POV[,-1])

ID.POV <- ID.POV[,order(colnames(ID.POV))]

POV <- ID.POV

POV.total <- apply(POV, 1, sum)

ID.FF <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = group, values_from = Count.FF,
              values_fn = sum, values_fill = 0)


ID.FF <- as.matrix(ID.FF[,-1])

ID.FF <- ID.FF[,order(colnames(ID.FF))]

FF <- ID.FF

FF.total <- apply(FF, 1, sum)

#Set up seat position covariate
ID.obs$Obs1.Front <- ifelse(ID.obs$Front == "BM", 1, 0)
Obs1.Front <- ID.obs %>% group_by(Transect) %>% summarise(Obs1.Front = round(mean(Obs1.Front)))
Obs1.Front <- Obs1.Front$Obs1.Front

ID.obs$Obs2.Front <- ifelse(ID.obs$Front == "TC", 1, 0)
Obs2.Front <- ID.obs %>% group_by(Transect) %>% summarise(Obs2.Front = round(mean(Obs2.Front)))
Obs2.Front <- Obs2.Front$Obs2.Front


ID.obs <- ID.obs[order(ID.obs$group),]
nspecies <- length(unique(ID.obs$group)) -1
nsites <- length(unique(ID.obs$numeric_tran))
nobs <- 2


##Less detailed model##
#-Nimble Code-#

code <- nimbleCode({
  
  #-Priors-#
  
  #Composition of latent abundance (corrected for imperfect detection)
  
  for(o in 1:nobs){
    phi[1:(nspecies+1), o] ~ ddirch(phi.ones[1:(nspecies+1),o])
    
    for(i in 1:(nspecies+1)){
      phi.ones[i,o] <- 1
    }
    
  }
  
  #Derived product of movement and detection
  for(o in 1:nobs){
    int.epsilon[o] ~ dnorm(0, 0.01)
  }
  
  beta.Obs1.Front ~ dnorm(0, 0.1)
  beta.Obs2.Front ~ dnorm(0, 0.1)
  
  #-Likelihood-#
  
  for(j in 1:nsites){
    
    #Front facing camera composition
    FF[j,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[j])
    
    #Front facing camera total abundance
    FF.total[j] ~ dpois(lambda.total)
    
    for(o in 1:nobs){
      
      OBS[j,1:(nspecies+1),o] ~ dmulti(phi[1:(nspecies+1),o], OBS.total[j,o])
      
      OBS.total[j,o] ~ dpois(lambda.total * E.epsilon[j,o])
      
    }#end o
    
  }#end j
  
  for(i in 1:nspecies){
    
    pi[i] <- lambda[i]/lambda.total
    
  
    
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
    for(o in 1:nobs){
      for(j in 1:nsites){
      
      log(E.epsilon[j,o]) <- int.epsilon[o] + beta.Obs1.Front*Obs1.Front[j] + beta.Obs2.Front*Obs2.Front[j]
      
      correction[j,i,o] <- E.epsilon[j,o] * phi[i,o]/pi[i]
      
    }
      
    }#end o
  }#end i
  
  
  lambda.total <- sum(lambda[1:nspecies])
  
})

#-Compile data-#

data <- list(FF = FF,
             FF.total = FF.total,
             Obs1.Front = Obs1.Front,
             Obs2.Front = Obs2.Front,
             OBS = OBS,
             OBS.total = OBS.total
)


constants <- list(nspecies = nspecies, nsites = nsites, nobs = 2)


#-Initial values-#

inits <- function(){list(pi = apply(FF/apply(FF, 1, sum), 2, mean, na.rm = TRUE),
                         # E.epsilon = rep(sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))
                         #                     *runif(1, 0.25, 1)), 2),
                         int.epsilon = sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))),
                         beta.Obs1.Front <- rnorm(1, 0, 0.1),
                         beta.Obs2.Front <- rnorm(1, 0, 0.1),
                         phi = apply(apply(OBS, c(1,2,3), max)/apply(apply(OBS, c(1,2,3), max), 1, sum), c(2,3), mean,
                                     na.rm = TRUE)
)}

#-Parameters to save-#

params <- c(
  "pi",
  "phi",
  "lambda",
  "lambda.total",
  "E.epsilon",
  "beta.Obs1.Front",
  "beta.Obs2.Front",
  "correction"
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

g.out <- runMCMC(compiled.model$MCMC,
                 niter = ni, nburnin = nb,
                 nchains = nc, thin = nt,
                 samplesAsCodaMCMC = TRUE)


#Diagnostics
g.out.mcmc <- as.mcmc.list(g.out)
g.out.ggs <- ggs(g.out.mcmc)

g.out.diag <- ggs_diagnostics(g.out.ggs)
g.out.gelman <- gelman.diag(g.out.mcmc)
MCMCtrace(g.out.mcmc, params = c("pi",
                                 "phi",
                                 "E.epsilon",
                                 "lambda",
                                 "lambda.total",
                                 "beta.Obs1.Front",
                                 "beta.Obs2.Front",
                                 "correction"
                                 ))

#Work with output
#Check dimensions - 20K iterations minus 10K for burn-in, 67 parameters
#Bring in model output
setwd("~/Documents/Windsor/UW Postdoc/Sea duck detection")
g.out <- readRDS("general_model.RDS")
dim(g.out[[1]])
g.mcmc.params <- as.mcmc.list(g.out)
g.params.groups <- data.frame(as.matrix(g.mcmc.params)) 

g.params.groups = as.matrix(g.params.groups)

g.out.groups <- data.frame(
  Mean = apply(g.params.groups, 2, mean),
  lcl = apply(g.params.groups, 2, quantile, probs = c(.05), na.rm = TRUE),
  ucl = apply(g.params.groups, 2, quantile, probs = c(.95), na.rm = TRUE),
  SD = apply(g.params.groups, 2, sd))

FF.g.sums <- apply(FF, 2, sum)

Species <- colnames(FF)


#Note, the following code is outdated and will be updated once we decide on a final model/model set
lambdas.g <- g.out.groups[c(11:18),]
lambda.mean.g <- nsites*lambdas.g$Mean

epsilons.g <- g.out.groups[c(3:10),1]
comp.ratio.g <- g.out.groups[c(36:51),1] / g.out.groups[c(60:67),1] #phi/psi

birds.adjust.g <- lambda.mean.g*(1/(comp.ratio.g*epsilons.g))

g.bird.groups <- as.data.frame(cbind(Species[1:8], FF.sums[1:8], round(lambda.mean.g, 2), round(birds.adjust.g, 2)))

colnames(g.bird.groups)[1:4] <- c("Species group","Counts from FF","Unadjusted N estimate", "Adjusted N estimate")
row.names(g.bird.groups) <- NULL






#Plot unknown species for each camera/observer
sb_df <- sb_df %>% filter(sb_df$SPECIES != "HAPO" & sb_df$SPECIES != "HASE" & sb_df$SPECIES != "BLLA"
                          & sb_df$SPECIES != "UNMA")

sb_df$Transect <- str_trunc(sb_df$TransectGroup, 12, "right", ellipsis = "")

sb_df <- sb_df %>% filter(platform == 1)

tran <- unique(sb_df$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))

sb_df$numeric_tran <- 0
for(i in dumb$tran){
  sb_df$numeric_tran[sb_df$Transect == i] <- match(i, dumb$tran)
}

SB.FF <- sb_df %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.FF,
              values_fn = sum, values_fill = 0)
unk.FF <- pivot_longer(data = SB.FF, cols = !numeric_tran, names_to = "Species", 
                                   values_to = "Count")

unk.FF <- unk.FF %>% filter(grepl("UN",Species) | grepl("UM",Species) | grepl("UL",Species) | Species == "USAC")

unk.FF.total <- unk.FF %>% group_by(Species) %>%
  summarise(AFF.Total = sum(Count))

SB.POV <- sb_df %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.POV,
              values_fn = sum, values_fill = 0)
unk.POV <- pivot_longer(data = SB.POV, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

unk.POV <- unk.POV %>% filter(grepl("UN",Species) | grepl("UM",Species) | grepl("UL",Species) | Species == "USAC")

unk.POV.total <- unk.POV %>% group_by(Species) %>%
  summarise(APOV.Total = sum(Count))

SB.BM <- sb_df %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.BM,
              values_fn = sum, values_fill = 0)
unk.BM <- pivot_longer(data = SB.BM, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

unk.BM <- unk.BM %>% filter(grepl("UN",Species) | grepl("UM",Species) | grepl("UL",Species) | Species == "USAC")

unk.BM.total <- unk.BM %>% group_by(Species) %>%
  summarise(BM.Total = sum(Count))

SB.TC <- sb_df %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.TC,
              values_fn = sum, values_fill = 0)
unk.TC <- pivot_longer(data = SB.TC, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

unk.TC <- unk.TC %>% filter(grepl("UN",Species) | grepl("UM",Species) | grepl("UL",Species) | Species == "USAC")

unk.TC.total <- unk.TC %>% group_by(Species) %>%
  summarise(TC.Total = sum(Count))

unknowns <- cbind(unk.FF.total, unk.POV.total, unk.BM.total, unk.TC.total)

unknowns <- pivot_longer(data = unknowns, cols = c(AFF.Total, APOV.Total, BM.Total, TC.Total), names_to = "Totals", 
                        values_to = "Counts")

ggplot() +
  geom_col(data = unknowns, aes(y = Counts, x = Species, fill = Totals), color = "black",
           position = position_dodge(width = 1.0)) +
  scale_fill_viridis(option = "magma", discrete = TRUE, direction = -1, name = "", labels = c("FF", "POV", "Obs1",
                                                                                              "Obs2")) +
  theme_bw() +
  scale_y_continuous(limits=c(0, 750), breaks=c(seq(0, 750, by=50))) +
  xlab("") +
  ylab("Total counts") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 6, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))

ggsave("Unknown species counts.tiff", plot = last_plot(), 
       path = "~/Documents/Windsor/UW Postdoc/Sea duck detection", dpi = 300)


AFF.sums <- apply(FF, 2, sum)
APOV.sums <- apply(POV, 2, sum)
Obs1.sums <- apply(OBS[,,1], 2, sum)
Obs2.sums <- apply(OBS[,,2], 2, sum)

duck.sums <- as.data.frame(cbind(AFF.sums, APOV.sums, Obs1.sums, Obs2.sums))
bg <- rownames(duck.sums)
duck.sums$Group <- bg
duck.sums <- pivot_longer(data = duck.sums, cols = c(AFF.sums, APOV.sums, Obs1.sums, Obs2.sums), 
                          names_to = "CamObs", values_to = "Counts")



ggplot() +
  geom_col(data = duck.sums, aes(y = Counts, x = Group, fill = CamObs), color = "black",
           position = position_dodge(width = 1.0)) +
  scale_fill_viridis(option = "plasma", discrete = TRUE, direction = -1, name = "", labels = c("FF", "POV", "Obs1",
                                                                                              "Obs2")) +
  theme_bw() +
  scale_y_continuous(limits=c(0, 4050), breaks=c(seq(0, 4050, by=150))) +
  xlab("") +
  ylab("Total counts") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 6, color = "black"),
        axis.title.x = element_text(size = 10, color = "black"),
        axis.title.y = element_text(size = 10, color = "black"),
        axis.line.x = element_line(color = "black", size = 0.5),
        axis.line.y = element_line(color = "black", size = 0.5))

ggsave("Species group counts.tiff", plot = last_plot(), 
       path = "~/Documents/Windsor/UW Postdoc/Sea duck detection", dpi = 300)



