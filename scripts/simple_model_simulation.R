set.seed(122)

library(tidyverse)
library(stringr)
library(runjags)
library(ggmcmc)
library(ggplot2)
library(VGAM)
library(pander)
library(nimble)

#define parameters and constants 
transects <- 20
species <- 4 
observers <- 2 



#generate a vector that looks like "tr" from model
#groups per transect <- poisson 
no.tr <- rpois(transects,5)
groups <- sum(no.tr)

no.gr <- rep(NA,groups)
for(i in 1:groups){
  no.gr[i] <- max(rpois(1,8),1)
}

#Parameters to provide with values - tau.ID, sigma.ID, tau.sp, sigma.sp, mn.pr.ID, lambda, beta, p


beta = c(1.19, 0.78, 1.07, 0.94)
lambda = 0.49
pr.ID = array(c(0.35, 0.01, 0.18, 0.00, 0.23, 0.97, 0.34, 0.35, 0.21, 0.01, 0.48,
                      0.65, 0.21, 0.01, 0.00, 0.00, 0.80, 0.00, 0.56, 0.37, 0.00, 0.71, 0.16,
                      0.40, 0.01, 0.16, 0.28, 0.00, 0.19, 0.13, 0.00, 0.23), dim = c(species, species, observers))
p = matrix(c(0.72, 0.67, 0.67, 0.70, 0.70, 0.72, 0.75, 0.59), nrow = species, ncol = observers)
alpha.misID = c(0.26, 0.20, 0.33, 0.21)






data_gen_full <- function(pr.spK, beta, species, transects, observers, groups, lambda, 
                          pr.obsJ.spK, p,  alpha.misID, M.obs, M, tr){
  
  
  FF <- POV <- matrix(NA, nrow = groups, ncol = species)
  
  for(i in 1:groups){
    FF[i,] <- rmultinom(1,no.gr[i],pr.spK)
    for(k in 1:species){
      POV[i,k] <- round(FF[i,k]*beta[k]*0.9)
    }
  }
  
  M <- apply(POV,1,sum)
  
  
  #for each species and observer, first redistribute the individuals in group according to a misID process for each 
  #species then do a binomial sample on the misIDed individuals for each observer 
  
  tr <- vector()
  tr <- rep(seq(20), no.tr)
  
  
  
  misID.obs <- array(NA, dim = c(groups, species, observers))
  for(i in 1:groups){
    for(k in 1:species){
      for(o in 1:observers){
        misID.obs[i,k,o] <- rbinom(1, tr[i], p[k,o])
      }
    }
  }
  
  
  pr.obsJ <- matrix(NA, nrow = species, ncol = observers)
  for(o in 1:observers){
    pr.obsJ[1,o] <- pr.spK[1]*pr.obsJ.spK[1,1,o] + pr.spK[2]*pr.obsJ.spK[1,2,o] +
      pr.spK[3]*pr.obsJ.spK[1,3,o] + pr.spK[4]*pr.obsJ.spK[1,4,o]
    
    pr.obsJ[2,o] <- pr.spK[1]*pr.obsJ.spK[2,1,o] + pr.spK[2]*pr.obsJ.spK[2,2,o] +
      pr.spK[3]*pr.obsJ.spK[2,3,o] + pr.spK[4]*pr.obsJ.spK[2,4,o]
    
    pr.obsJ[3,o] <- pr.spK[1]*pr.obsJ.spK[3,1,o] + pr.spK[2]*pr.obsJ.spK[3,2,o] +
      pr.spK[3]*pr.obsJ.spK[3,3,o] + pr.spK[4]*pr.obsJ.spK[3,4,o]
    
    pr.obsJ[4,o] <- pr.spK[1]*pr.obsJ.spK[4,1,o] + pr.spK[2]*pr.obsJ.spK[4,2,o] +
      pr.spK[3]*pr.obsJ.spK[4,3,o] + pr.spK[4]*pr.obsJ.spK[4,4,o]
  }
  
  
  misID1 <- matrix(NA, nrow = groups, ncol = species)
  
  for(i in 1:groups){
    for(k in 1:species){
      misID1[i,k] <- rbinom(1,tr[i],pr.obsJ[k,1]) 
    }
  }
  
  M.obs1 <- apply(misID1,1,sum)
  
  misID2 <- matrix(NA, nrow = groups, ncol = species)
  
  for(i in 1:groups){
    for(k in 1:species){
      misID2[i,k] <- rbinom(1,tr[i],pr.obsJ[k,1]) 
    }
  }
  
  M.obs2 <- apply(misID2,1,sum)
  
  
  misID <- array(c(misID1, misID2), dim = c(groups, species, observers))
  data.sim <- data.frame(misID1, misID2, misID.obs[,,1], misID.obs[,,2], POV, FF, M.obs1, M.obs2, M, tr)
  colnames(data.sim) <- c("Obs1Spp1", "Obs1Spp2", "Obs1Spp3", "Obs1Spp4", "Obs2Spp1", "Obs2Spp2",
                          "Obs2Spp3", "Obs2Spp4", "IDObs1.spp1", "IDObs1.spp2", "IDObs1.spp3", "IDObs1.spp4",
                          "IDObs2.spp1", "IDObs2.spp2", "IDObs2.spp3", "IDObs2.spp4", "POV.spp1", "POV.spp2",
                          "POV.spp3", "POV.spp4", "FF.spp1", "FF.spp2", "FF.spp3", "FF.spp4", "M.obs1", "M.obs2", "M",
                          "tr")
  
  return(data.sim)
}


ID.det <- function(n_iter, pr.spK, beta, species, transects, observers, groups, lambda, 
                   pr.obsJ.spK, p,  alpha.misID, M.obs, M, tr,
                   n.chains = 3, adapt = 5000, burnin = 5000,  sample = 5000
){
  s.mod_summ <- as.list(1:n_iter)
  computation_time <- as.list(1:n_iter)
  mcmc_results <- as.list(1:n_iter)
  sim.dat_list <- as.list(1:n_iter)
  
  # set progress bar
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  
  for(n in 1:n_iter) {
    data.sim <- data_gen_full(beta = beta, 
                              p = p, 
                              species = species,
                              transects = transects,
                              groups = groups,
                              observers = observers,
                              pr.spK = pr.spK,
                              M.obs = M.obs,
                              M = M,
                              tr = tr,
                              pr.obsJ.spK = pr.obsJ.spK)
    FF <- cbind(data.sim$FF.spp1, data.sim$FF.spp2, data.sim$FF.spp3, data.sim$FF.spp4)
    POV <- cbind(data.sim$POV.spp1, data.sim$POV.spp2, data.sim$POV.spp3, data.sim$POV.spp4)
    ID.obs1 <- cbind(data.sim$Obs1Spp1, data.sim$Obs1Spp2, data.sim$Obs1Spp3, data.sim$Obs1Spp4)
    ID.obs2 <- cbind(data.sim$Obs2Spp1, data.sim$Obs2Spp2, data.sim$Obs2Spp3, data.sim$Obs2Spp4)
    misID <- array(c(ID.obs1, ID.obs2), dim = c(groups, species, observers))
    misID.obs1 <- cbind(data.sim$IDObs1.spp1, data.sim$IDObs1.spp2, data.sim$IDObs1.spp3, data.sim$IDObs1.spp4)
    misID.obs2 <- cbind(data.sim$IDObs2.spp1, data.sim$IDObs2.spp2, data.sim$IDObs2.spp3, data.sim$IDObs2.spp4)
    misID.obs <- array(c(misID.obs1, misID.obs2), dim = c(groups, species, observers))
    M.obs <- cbind(data.sim$M.obs1, data.sim$M.obs2)
    M <- data.sim$M
    tr <- data.sim$tr
    
    sim.params <- c("beta","p")
    
    data <- list(ID = misID, M.obs = M.obs, obs = misID.obs, POV = POV, M = M, FF = FF, tr = tr, 
                 n.transects = transects, n.groups = groups, n.observers = observers, n.species = species, 
                 alpha.sp = c(1,1,1,1))
    init.pr.spK <- matrix(NA,nrow=transects,ncol = species)
    for(t in 1:transects){
      init.pr.spK[t,] <- c(0.25,0.25,0.25,0.25)
    }
    init.lambda <- matrix(NA,nrow = groups,ncol = species)
    for(i in 1:groups){
      for(j in 1:species){
        init.lambda[i,j] <- max(misID.obs[i,j,])+1 
      }
    }
    inits.mult<-function(){list(pr.spK = init.pr.spK, beta = rep(1,4), lambda = init.lambda)}
    
    start <- Sys.time()
    setwd("~/Documents/Windsor/UW Postdoc/Sea duck detection")
    mult.mod.sim <- run.jags(model = "MultinomialSimModel.txt",
                             monitor = sim.params,
                             data = data,
                             n.chains = n.chains,
                             adapt = adapt,
                             burnin = burnin,
                             sample = sample, 
                             inits = inits.mult)
    computation_time[n] <- Sys.time() - start 
    
    # store summary output for iteration i in list  
    post_summ <- as.data.frame(summary(mult.mod.sim))
    
    # compare posterior summaries to data generating values to track coverage 
    post_summ$iter <- n
    post_summ$par <- c(rownames(post_summ)) 
    # need to test if the data generating values are between the 
    # bounds for the 95% CrI, not the posterior means...
    dg_values <- c(beta, p)
    
    post_summ$capture <- as.numeric(dg_values >= post_summ$Lower95 & dg_values <= post_summ$Upper95)
    post_summ$cri95_width <- post_summ$Upper95 - post_summ$Lower95
    
    #save the summary output from this realization in the list
    s.mod_summ[[n]] <- post_summ
    mcmc_results[[n]] <- mult.mod.sim$mcmc
    sim.dat_list[[n]] <- data.sim
    setTxtProgressBar(pb, n)
    print(paste("Completed iteration", n, "of", n_iter))
  }
  close(pb)
  
  # return the list of summary output from each run of the JAGS model
  jags_summ <- do.call(rbind, s.mod_summ)
  out <- list(jags_summ = jags_summ, mcmc_results = mcmc_results, 
              sim.dat_list = sim.dat_list, 
              comp_time = computation_time)
}

start <- Sys.time()
sim.det.ID <- ID.det(n_iter = 10, burnin = 5000,  sample = 5000, pr.spK = pr.spK, beta = beta, p = p, lambda = lambda,
                     M.obs = M.obs, M = M, tr = tr, groups = groups, species = species, transects = transects, 
                     observers = observers, pr.obsJ.spK = pr.obsJ.spK)

Sys.time() - start


sim.info <- sim.det.ID$jags_summ %>% 
  group_by(par) %>% 
  summarise(cap_rate = mean(capture), 
            avg_est = mean(Mean), 
            ci_width = mean(Upper95 - Lower95), 
            avg_L95 = mean(Lower95), 
            avg_U95 = mean(Upper95),
            avg_rhat = mean(psrf), 
            avg_neff = mean(SSeff))

sim.info$dg_values <- c(1.19, 0.78, 1.07, 0.94, 0.72, 0.67, 0.67, 0.70, 0.70, 0.72, 0.75, 0.59)
sim.info$avg_bias <- sim.info$avg_est - sim.info$dg_values
pander::pander(sim.info, caption = "simualtion results")









