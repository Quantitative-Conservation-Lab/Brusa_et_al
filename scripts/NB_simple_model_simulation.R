set.seed(122)

library(tidyverse)
library(stringr)
library(jagsUI)
library(ggmcmc)
library(ggplot2)
library(VGAM)
library(pander)
library(here)

#define parameters and constants 
transects <- 20
species <- 4 
observers <- 2 



#generate a vector that looks like "tr" from model
no.tr <- rpois(transects,5)
groups <- sum(no.tr)

no.gr <- rep(NA,groups)
for(i in 1:groups){
  no.gr[i] <- max(rpois(1,8),1)
}


beta = c(1.19, 0.78, 1.07, 0.94)
p = matrix(c(0.72, 0.67, 0.67, 0.70, 0.70, 0.72, 0.75, 0.59), nrow = species, ncol = observers)
pr.sp = c(0.30, 0.15, 0.20, 0.35)
pr.ID1 = c(0.25, 0.10, 0.15, 0.40)
pr.ID2 = c(0.30, 0.20, 0.15, 0.30)



data_gen_full <- function(beta, p, pr.sp, pr.ID1, pr.ID2, species, transects, observers, groups, 
                          M.mis, M.obs.mis1, M.obs.mis2, tr){
  
  
  FF.det <- POV.mis <- matrix(NA, nrow = groups, ncol = species)
  
  for(i in 1:groups){
    FF.det[i,] <- rmultinom(1,no.gr[i],pr.sp)
    for(k in 1:species){
      POV.mis[i,k] <- round(FF.det[i,k]*beta[k]*0.9)
    }
  }
  
  
  M.mis <- apply(POV.mis,1,sum)
  
  
  #for each species and observer, first redistribute the individuals in group according to a misID process for each 
  #species then do a binomial sample on the misIDed individuals for each observer 
  
  tr <- vector()
  tr <- rep(seq(20), no.tr)
  
  
  ID.mis1 <- matrix(NA, nrow = groups, ncol = species)
  
  for(i in 1:groups){
    for(k in 1:species){
      ID.mis1[i,k] <- rbinom(1,tr[i],pr.ID1[k]) 
    }
  }
  
  M.obs.mis1 <- apply(ID.mis1,1,sum)
  
  ID.mis2 <- matrix(NA, nrow = groups, ncol = species)
  
  for(i in 1:groups){
    for(k in 1:species){
      ID.mis2[i,k] <- rbinom(1,tr[i],pr.ID2[k]) 
    }
  }
  
  M.obs.mis2 <- apply(ID.mis2,1,sum)
  
  
  ID.mis <- array(c(ID.mis1, ID.mis2), dim = c(groups, species, observers))
  
  beta.2 <- matrix(beta, nrow = species, ncol = observers)
  
  ID.det <- array(NA, dim = c(groups, species, observers))
  for(i in 1:groups){
    for(j in 1:species){
      for(o in 1:observers){
        ID.det[i,j,o] <- rbinom(1, tr[i], p[j,o])
      }
    }
  }
  
  data.sim <- data.frame(ID.mis1, ID.mis2, ID.det[,,1], ID.det[,,2], POV.mis, FF.det, M.obs.mis1, M.obs.mis2, M.mis, tr)
  colnames(data.sim) <- c("Obs1Spp1.mis", "Obs1Spp2.mis", "Obs1Spp3.mis", "Obs1Spp4.mis", "Obs2Spp1.mis", "Obs2Spp2.mis",
                          "Obs2Spp3.mis", "Obs2Spp4.mis", "Obs1.spp1.det", "Obs1.spp2.det", "Obs1.spp3.det", 
                          "Obs1.spp4.det", "Obs2.spp1.det", "Obs2.spp2.det", "Obs2.spp3.det", "Obs2.spp4.det", 
                          "POV.spp1", "POV.spp2", "POV.spp3", "POV.spp4", "FF.spp1", "FF.spp2", "FF.spp3", "FF.spp4", 
                          "M.obs.mis1", "M.obs.mis2", "M.mis", "tr")
  
  return(data.sim)
}


ID.det <- function(n_iter, beta, p, pr.sp, pr.ID1, pr.ID2, species, transects, observers, groups, 
                   M.mis, M.obs.mis1, M.obs.mis2, tr,
                   n.chains = 3, adapt = 5000, burnin = 5000,  sample = 20000
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
                              pr.sp = pr.sp,
                              pr.ID1 = pr.ID1,
                              pr.ID2 = pr.ID2,
                              species = species,
                              transects = transects,
                              observers = observers,
                              groups = groups,
                              M.mis = M.mis,
                              M.obs.mis1 = M.obs.mis1,
                              M.obs.mis2 = M.obs.mis2,
                              tr = tr)
    FF.det <- cbind(data.sim$FF.spp1, data.sim$FF.spp2, data.sim$FF.spp3, data.sim$FF.spp4)
    POV.mis <- cbind(data.sim$POV.spp1, data.sim$POV.spp2, data.sim$POV.spp3, data.sim$POV.spp4)
    ID.mis1 <- cbind(data.sim$Obs1Spp1.mis, data.sim$Obs1Spp2.mis, data.sim$Obs1Spp3.mis, data.sim$Obs1Spp4.mis)
    ID.mis2 <- cbind(data.sim$Obs2Spp1.mis, data.sim$Obs2Spp2.mis, data.sim$Obs2Spp3.mis, data.sim$Obs2Spp4.mis)
    ID.mis <- array(c(ID.mis1, ID.mis2), dim = c(groups, species, observers))
    ID.det1 <- cbind(data.sim$Obs1.spp1.det, data.sim$Obs1.spp2.det, data.sim$Obs1.spp3.det, data.sim$Obs1.spp4.det)
    ID.det2 <- cbind(data.sim$Obs2.spp1.det, data.sim$Obs2.spp2.det, data.sim$Obs2.spp3.det, data.sim$Obs2.spp4.det)
    ID.det <- array(c(ID.det1, ID.det2), dim = c(groups, species, observers))
    M.obs.mis <- cbind(data.sim$M.obs.mis1, data.sim$M.obs.mis2)
    M.mis <- data.sim$M.mis
    tr <- data.sim$tr
    
    pr.ID.pred <- apply(cbind(pr.ID1, pr.ID2), 1, mean)
    
    pr.sp.pred <- pr.sp
    
    sim.params <- c("beta", "p", "pr.sp.pred", "pr.ID.pred")
    
    data <- list(ID.mis = ID.mis, M.obs.mis = M.obs.mis, ID.det = ID.det, POV.mis = POV.mis, M.mis = M.mis, 
                 FF.det = FF.det, tr.mis = tr, n.transects.mis = transects, n.groups.mis = groups, n.groups.det = groups,
                 n.observers = observers, n.species = species)
    inits.duck<-function(){list(N.obsJ=(apply(ID.det,c(1,2),sum)+matrix(rep(1,groups*species),
                                                                        nrow=groups,ncol=species)))}
    
    start <- Sys.time()
    #Generate samples from the posterior distribution
    simp.mod.sim = jagsUI::jags(data, inits.duck, sim.params, model.file=here("scripts", "NBSimpleModSimModel.txt"),
                                n.chains=n.chains, 
                                n.iter=sample, 
                                n.burnin=burnin, 
                                n.thin=1)
    
    computation_time[n] <- Sys.time() - start 
    
    # store summary output for iteration i in list  
    post_summ <- as.data.frame(simp.mod.sim$summary)
    
    #Remove deviance row
    d<-dim(post_summ)[1]
    post_summ<-post_summ[1:(d-1),]
    
    # compare posterior summaries to data generating values to track coverage 
    post_summ$iter <- n
    post_summ$par <- c(rownames(post_summ)) 
    # need to test if the data generating values are between the 
    # bounds for the 95% CrI, not the posterior means...
    dg_values <- c(beta, p, pr.sp.pred, pr.ID.pred)
    
    post_summ$capture <- as.numeric(dg_values >= post_summ$`2.5%` & dg_values <= post_summ$`97.5%`)
    post_summ$cri95_width <- post_summ$`97.5%` - post_summ$`2.5%`
    
    #save the summary output from this realization in the list
    s.mod_summ[[n]] <- post_summ
    mcmc_results[[n]] <- simp.mod.sim$mcmc
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
sim.det.ID <- ID.det(n_iter = 100, burnin = 5000,  sample = 20000, beta = beta, p = p, pr.sp = pr.sp, pr.ID1 = pr.ID1,
                     pr.ID2 = pr.ID2, M.obs.mis1 = M.obs.mis, M.obs.mis2 = M.obs.mis2, M.mis = M.mis, tr = tr, 
                     groups = groups, species = species, transects = transects, observers = observers)

Sys.time() - start


sim.info <- sim.det.ID$jags_summ %>% 
  group_by(par) %>% 
  summarise(cap_rate = mean(capture), 
            avg_est = mean(mean), 
            ci_width = mean(`97.5%` - `2.5%`), 
            avg_L95 = mean(`2.5%`), 
            avg_U95 = mean(`97.5%`),
            avg_rhat = mean(Rhat), 
            avg_neff = mean(n.eff))

sim.info$dg_values <- c(1.19, 0.78, 1.07, 0.94, 0.72, 0.67, 0.67, 0.70, 0.70, 0.72, 0.75, 0.59,
                        0.275, 0.150, 0.150, 0.350, 0.30, 0.15, 0.20, 0.35)
sim.info$avg_bias <- sim.info$avg_est - sim.info$dg_values
pander::pander(sim.info, caption = "simualtion results")








