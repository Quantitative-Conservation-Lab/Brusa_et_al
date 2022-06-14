library(tidyverse)
library(nimble)
library(ggmcmc)
library(coda)
library(MCMCvis)
library(here)

sb_df <- read.csv(here('Data', 'seaducks.csv'))

#Mis-Identification by group, then by individual, then integrate mis-ID info into detection model
#Remove uninformative records, all shorebirds have been removed
sb_df <- sb_df %>% filter(sb_df$SPECIES != "PEFA" & sb_df$SPECIES != "UNSD" & sb_df$SPECIES != "BAEA" 
                          & sb_df$SPECIES != "BLOY" & sb_df$SPECIES != "GBHE" & sb_df$SPECIES != "NWCR"
                          & sb_df$SPECIES != "HAPO" & sb_df$SPECIES != "HASE" & sb_df$SPECIES != "BLLA"
                          & sb_df$SPECIES != "UNMA" & sb_df$SPECIES !=  "UMSD"
                          & sb_df$SPECIES != "UNDU" & sb_df$SPECIES !=  "USSD")

setwd("~/Documents/Windsor/UW Postdoc/Sea duck detection")
bird.groups <- read.csv("Marine bird groups.csv")

#Work with just species in dataset
bird.groups <- bird.groups %>% filter(Bird.species..1.0. == 1 & In.the.dataset..1.0. == 1)

#Investigate group sizes
bird.group.totals <- bird.groups %>% group_by(Label) %>%
  summarise(Total.FF = sum(Total.FF), Total.POV = sum(Total.POV), Total.Obs1 = sum(Total.Obs.1),
            Total.Obs2 = sum(Total.Obs.2))

Total.FF <- bird.group.totals$Total.FF

bird.group.totals <- bird.group.totals[order(Total.FF),]




#Larger groups
#Probabily will want to remove dove - not a species of interest and not enough observations
groups_df <- sb_df %>%
  mutate(case_when(str_detect(sb_df$SPECIES, "SC$") | str_detect(sb_df$SPECIES, "BM$") | 
                     str_detect(sb_df$SPECIES, "UNGO") | str_detect(sb_df$SPECIES, "COGO") | 
                     str_detect(sb_df$SPECIES, "BAGO") | str_detect(sb_df$SPECIES, "ME$") |
                   str_detect(sb_df$SPECIES, "HADU") | str_detect(sb_df$SPECIES, "DS$") |
                     str_detect(sb_df$SPECIES, "RUDU") | str_detect(sb_df$SPECIES, "SCAU") |
                     str_detect(sb_df$SPECIES, "UNDD") ~ "diving duck", 
                   str_detect(sb_df$SPECIES, "LO$") ~ "loon", str_detect(sb_df$SPECIES, "GR$") ~ "grebe", 
                   str_detect(sb_df$SPECIES, "UNGU") | str_detect(sb_df$SPECIES, "GWGU") | 
                     str_detect(sb_df$SPECIES, "WG$") ~ "gull", str_detect(sb_df$SPECIES, "UNSB") |
                   str_detect(sb_df$SPECIES, "CO$") ~ "sea bird", str_detect(sb_df$SPECIES, "LL$") |
                   str_detect(sb_df$SPECIES, "WI$") |str_detect(sb_df$SPECIES, "PI$") | str_detect(sb_df$SPECIES, "UNPD")
                   ~ "dabbling duck", 
                   str_detect(sb_df$SPECIES, "BR$") | str_detect(sb_df$SPECIES, "CAGO") ~ "waterfowl",
                    str_detect(sb_df$SPECIES, "CAAU") | 
                     str_detect(sb_df$SPECIES, "PAAU") | str_detect(sb_df$SPECIES, "RHAU") | 
                     str_detect(sb_df$SPECIES, "AC$") | str_detect(sb_df$SPECIES, "MU$") | 
                     str_detect(sb_df$SPECIES, "PIGU") | str_detect(sb_df$SPECIES, "ML$") ~ "alcid",  
                    str_detect(sb_df$SPECIES, "DO$") ~ "dove"))


colnames(groups_df)[30] <- "sp.group"

groups_df <- groups_df %>% filter(platform == 1)
groups_df <- groups_df %>% filter(sp.group != "dove")


##Prep data for model##

#Set up array for observer records
groups_df$Transect <- str_trunc(groups_df$TransectGroup, 12, "right", ellipsis = "")
ID.obs <- groups_df

tran <- unique(ID.obs$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))

ID.obs$numeric_tran <- 0
for(i in dumb$tran){
  ID.obs$numeric_tran[ID.obs$Transect == i] <- match(i, dumb$tran)
}


ID.BM <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = sp.group, values_from = Count.BM,
              values_fn = sum, values_fill = 0)


ID.BM <- as.matrix(ID.BM[,-1])

ID.TC <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = sp.group, values_from = Count.TC,
              values_fn = sum, values_fill = 0)


ID.TC <- as.matrix(ID.TC[,-1])

OBS <- array(c(ID.BM, ID.TC), dim = c(nrow(ID.BM), ncol(ID.BM), 2))

BM.total <- apply(ID.BM, 1, sum)

TC.total <- apply(ID.TC, 1, sum)

OBS.total <- cbind(BM.total, TC.total)


ID.POV <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = sp.group, values_from = Count.POV,
              values_fn = sum, values_fill = 0)


ID.POV <- as.matrix(ID.POV[,-1])

POV <- ID.POV

POV.total <- apply(POV, 1, sum)

ID.FF <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = sp.group, values_from = Count.FF,
              values_fn = sum, values_fill = 0)


ID.FF <- as.matrix(ID.FF[,-1])

FF <- ID.FF

FF.total <- apply(FF, 1, sum)


nspecies <- length(unique(ID.obs$sp.group))
nsites <- length(unique(ID.obs$numeric_tran))
nobs <- 2


##Less detailed model##
#-Nimble Code-#

code <- nimbleCode({
  
  #-Priors-#
  
  #Composition of point of view camera
  psi[1:nspecies] ~ ddirch(psi.ones[1:nspecies])
  
  #Composition of latent abundance (corrected for imperfect detection)
  #phi[1:nspecies] ~ ddirch(phi.ones[1:nspecies])
  for(o in 1:nobs){
    phi[1:nspecies, o] ~ ddirch(phi.ones[1:nspecies,o])
  }
  
  #Derived product of movement and detection
  for(o in 1:nobs){
    E.epsilon[o] ~ dnorm(0, 0.01)
  }
  
  #-Likelihood-#
  
  for(j in 1:nsites){
    
    #Front facing camera composition
    FF[j,1:nspecies] ~ dmulti(pi[1:nspecies], FF.total[j])
    
    #Point of view camera composition
    POV[j,1:nspecies] ~ dmulti(psi[1:nspecies], POV.total[j])
    
    #Front facing camera total abundance
    FF.total[j] ~ dpois(lambda.total)
    
    for(o in 1:nobs){
      
      OBS[j,1:nspecies,o] ~ dmulti(phi[1:nspecies,o], OBS.total[j,o])
      
      OBS.total[j,o] ~ dpois(lambda.total * E.epsilon[o])
      
    }#end o
    
  }#end j
  
  for(i in 1:nspecies){
    
    pi[i] <- lambda[i]/lambda.total
    
    psi.ones[i] <- 1
    
    for(o in 1:nobs){
      phi.ones[i,o] <- 1
    }
    
    #phi.ones[i] <- 1 
   }

     
  for(i in 1:nspecies){
    log(lambda[i]) <- lambda0[i]
    lambda0[i] ~ dnorm(0, 0.01)
    
    for(o in 1:nobs){
    epsilon[i] <- psi[i] * E.epsilon[o] / (pi[i])
    misID[i,o] <- phi[i,o]/psi[i]
    }#end o
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
                         # E.epsilon = matrix(ncol = nobs, nrow = 1, 
                         #                    sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))
                         #                 *runif(1, 0.25, 1))),
                         E.epsilon = rep(sum(pi * (rnorm(nspecies, runif(1, 0.5, 1), 0.1)*runif(1,1,1.5))
                                         *runif(1, 0.25, 1)), 2),
                         epsilon = rnorm(nspecies, runif(1, 0.5, 1), 0.1) * runif(1,1,1.5) * runif(1, 0.25, 1),
                         epsilon = rnorm(nspecies, runif(1, 0.5, 1), 0.1) * runif(1,1,1.5) * runif(1, 0.25, 1),
                         psi = apply(POV/apply(POV, 1, sum), 2, mean, na.rm = TRUE),
                         phi = apply(apply(OBS, c(1,2,3), max)/apply(apply(OBS, c(1,2,3), max), 1, sum), c(2,3), mean,
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
  "E.epsilon",
  "misID"
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


#Diagnostics
out.mcmc <- as.mcmc.list(out)
out.ggs <- ggs(out.mcmc)

out.diag <- ggs_diagnostics(out.ggs)
out.gelman <- gelman.diag(out.mcmc)
MCMCtrace(out.mcmc, params = c("pi",
                               "psi",
                               "phi",
                               "lambda",
                               "lambda.total",
                               "epsilon",
                               "E.epsilon",
                               "misID" ))

#Work with output
#Check dimensions - 20K iterations minus 10K for burn-in, 67 parameters
dim(out[[1]])
mcmc.params <- as.mcmc.list(out)
params.groups <- data.frame(as.matrix(mcmc.params)) 

params.groups = as.matrix(params.groups)

out.groups <- data.frame(
  Mean = apply(params.groups, 2, mean),
  lcl = apply(params.groups, 2, quantile, probs = c(.05)),
  ucl = apply(params.groups, 2, quantile, probs = c(.95)),
  SD = apply(params.groups, 2, sd))

out.g.summ <- MCMCsummary(mcmc.params)

FF.g.sums <- apply(FF, 2, sum)

Species <- colnames(FF)

lambdas.g <- out.groups[c(11:18),]
lambda.mean.g <- nsites*lambdas.g$Mean

epsilons.g <- out.groups[c(2:9),1]
comp.ratio.g <- out.groups[c(30:37),1] / out.groups[c(48:55),1]

birds.adjust.g <- lambda.mean.g*(1/(comp.ratio.g*epsilons.g))

bird.groups <- as.data.frame(cbind(Species[1:8], FF.sums[1:8], round(lambda.mean.g, 2), round(birds.adjust.g, 2)))

colnames(bird.groups)[1:4] <- c("Species group","Counts from FF","Unadjusted N estimate", "Adjusted N estimate")
row.names(bird.groups) <- NULL









####
##Data for genus-level modeling##
sb_df_g <- read.csv(here('Data', 'seaducks.csv'))

#Remove uninformative records:
#UNDD - unknown diving duck, HAPO - harbor porpoise, GBHE - great blue heron, UNSD - sandpiper, HASE - harbor seal,
#UNDU - unknown duck, UNPD - unknown duck, UNSB - unknown seabird, BAEA - bald eagle, NWCR - northwestern crow,
#BLOY - black oystercatcher, UMSD - unknown shorebird, PEFA - peregrine falcon, BLLA - black lab, UNMA - ??

sb_df_g <- sb_df_g %>% filter(sb_df_g$SPECIES != "PEFA" & sb_df_g$SPECIES != "UNSD" & sb_df_g$SPECIES != "BAEA" 
                              & sb_df_g$SPECIES != "BLOY" & sb_df_g$SPECIES != "GBHE" & sb_df_g$SPECIES != "NWCR"
                              & sb_df_g$SPECIES != "HAPO" & sb_df_g$SPECIES != "HASE" & sb_df_g$SPECIES != "BLLA"
                              & sb_df_g$SPECIES != "UNMA" & sb_df_g$SPECIES !=  "UNPD" & sb_df_g$SPECIES !=  "UMSD" 
                              & sb_df_g$SPECIES != "UNDD" & sb_df_g$SPECIES != "UNDU" & sb_df_g$SPECIES !=  "USSD"
                              & sb_df_g$SPECIES != "UNSB")


#Smaller groups
small.groups <- sb_df_g %>%
  mutate(case_when(str_detect(sb_df$SPECIES, "SC$") ~ "Scoter", str_detect(sb_df$SPECIES, "LO$") ~ "Loon", 
                   str_detect(sb_df$SPECIES, "GR$") ~ "Grebe", str_detect(sb_df$SPECIES, "ME$") ~ "Merganser", 
                   str_detect(sb_df$SPECIES, "UNGO") | str_detect(sb_df$SPECIES, "COGO") | 
                     str_detect(sb_df$SPECIES, "BAGO") ~ "Goldeneye", str_detect(sb_df$SPECIES, "BM$") ~ "BuBMlehead",
                   str_detect(sb_df$SPECIES, "UNGU") | str_detect(sb_df$SPECIES, "GWGU") | 
                     str_detect(sb_df$SPECIES, "WG$") ~ "Gull",
                   str_detect(sb_df$SPECIES, "MU$") | str_detect(sb_df$SPECIES, "ML$") ~ "Murrelet",
                   str_detect(sb_df$SPECIES, "CO$") ~ "Cormorant", str_detect(sb_df$SPECIES, "LL$") ~ "Mallard",
                   str_detect(sb_df$SPECIES, "WI$") ~ "Wigeon", 
                   str_detect(sb_df$SPECIES, "PI$") ~ "Pintail", str_detect(sb_df$SPECIES, "BR$") ~ "Brant",
                   str_detect(sb_df$SPECIES, "DS$") ~ "Long-tailed Duck",
                   str_detect(sb_df$SPECIES, "RUDU") ~ "Ruddy Duck", str_detect(sb_df$SPECIES, "CAAU") | 
                     str_detect(sb_df$SPECIES, "PAAU") | str_detect(sb_df$SPECIES, "RHAU") | 
                     str_detect(sb_df$SPECIES, "AC$") ~ "Alcid", str_detect(sb_df$SPECIES, "PIGU") ~ "Guillemot",  
                   str_detect(sb_df$SPECIES, "HADU") ~ "Harelquin Duck", str_detect(sb_df$SPECIES, "DO$") ~ "Dove",
                   str_detect(sb_df$SPECIES, "CAGO") ~ "Canadian Goose", str_detect(sb_df$SPECIES, "SCAU") ~ "Scaup"))


#fitting one model for each observer...better than running two models to not have
colnames(small.groups)[30] <- "Spp.Grp"

#Set up array for observer records
sb_df_g$Transect <- str_trunc(sb_df_g$TransectGroup, 12, "right", ellipsis = "")
ID.obs <- sb_df_g

tran <- unique(ID.obs$Transect)
dumb <- data.frame(tran = tran, 
                   numbers = 1:length(tran))

ID.obs$numeric_tran <- 0
for(i in dumb$tran){
  ID.obs$numeric_tran[ID.obs$Transect == i] <- match(i, dumb$tran)
}

ID.BM <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Spp.Grp, values_from = Count.BM,
              values_fn = sum, values_fill = 0)

ID.BM <- as.matrix(ID.BM[,-1])

ID.TC <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Spp.Grp, values_from = Count.TC,
              values_fn = sum, values_fill = 0)

ID.TC <- as.matrix(ID.TC[,-1])

OBS <- array(c(ID.BM, ID.TC), dim = c(nrow(ID.BM), ncol(ID.BM), 2))

BM.total <- apply(ID.BM, 1, sum)

TC.total <- apply(ID.TC, 1, sum)

OBS.total <- cbind(BM.total, TC.total)

ID.POV <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Spp.Grp, values_from = Count.POV,
              values_fn = sum, values_fill = 0)

ID.POV <- as.matrix(ID.POV[,-1])

POV <- ID.POV

POV.total <- apply(POV, 1, sum)

ID.FF <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = Spp.Grp, values_from = Count.FF,
              values_fn = sum, values_fill = 0)

ID.FF <- as.matrix(ID.FF[,-1])

FF <- ID.FF

FF.total <- apply(FF, 1, sum)

nspecies <- length(unique(ID.obs$Spp.Grp))
nsites <- length(unique(ID.obs$numeric_tran))
nobs <- 2



#Plot unknown species for each camera/observer
SB.FF <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.FF,
              values_fn = sum, values_fill = 0)
unk.FF <- pivot_longer(data = SB.FF, cols = !numeric_tran, names_to = "Species", 
                                   values_to = "Count")

unk.FF <- unk.FF %>% filter(grepl("UN",Species))

unk.FF.total <- unk.FF %>% group_by(Species) %>%
  summarise(AFF.Total = sum(Count))

SB.POV <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.POV,
              values_fn = sum, values_fill = 0)
unk.POV <- pivot_longer(data = SB.POV, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

unk.POV <- unk.POV %>% filter(grepl("UN",Species))

unk.POV.total <- unk.POV %>% group_by(Species) %>%
  summarise(APOV.Total = sum(Count))

SB.BM <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.BM,
              values_fn = sum, values_fill = 0)
unk.BM <- pivot_longer(data = SB.BM, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

unk.BM <- unk.BM %>% filter(grepl("UN",Species))

unk.BM.total <- unk.BM %>% group_by(Species) %>%
  summarise(BM.Total = sum(Count))

SB.TC <- ID.obs %>%
  pivot_wider(id_cols = c(numeric_tran), names_from = SPECIES, values_from = Count.TC,
              values_fn = sum, values_fill = 0)
unk.TC <- pivot_longer(data = SB.TC, cols = !numeric_tran, names_to = "Species", 
                       values_to = "Count")

unk.TC <- unk.TC %>% filter(grepl("UN",Species))

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
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        panel.border = element_blank(),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 7, color = "black"),
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
  scale_y_continuous(limits=c(0, 2500), breaks=c(seq(0, 2500, by=150))) +
  xlab("") +
  ylab("Total counts") +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
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



