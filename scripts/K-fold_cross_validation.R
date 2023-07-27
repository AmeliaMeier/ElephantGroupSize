library(tidyverse)
library(jagsUI)
library(groupdata2)
library(overlapping)
library(lattice)


set.seed(333)

##### Model #####
sink("model.txt")
cat("
    model{
    ### priors and constraints
    ## Detection Model
    mu.focal ~ dunif(0, 10)
    sigma.focal ~ dunif(0, 100)
    tau.focal <- 1/(sigma.focal*sigma.focal)
    
    for (f in 1:n.focals){
       p.alpha[f]~dnorm(mu.focal, tau.focal)
    }
    
    # Site
    p.site ~ dnorm(0, 0.01)
    
    ## GroupSize
    psi <- sum(lam[])/M
    
    # alpha
    alpha~ dnorm(0, 0.01)
    
    # intercept prior for sex(fixed)
    b.sex~dnorm(0, 0.01)
    
    # intercept prior for fruit
    b.fruit ~ dnorm(0, 0.01)
    
    # evi 
    b.evi ~ dnorm(0, 0.01)
    
    # distance to water
    b.wat ~ dnorm(0, 0.01)
    
    #hfi prior
    b.hfi ~ dnorm(0, 0.01)
    
    # ppt prior
    b.ppt ~dnorm(0, 0.01)
    

    ### abundance model
    for(s in 1:S){
    log(lam[s])<- alpha + b.sex*sex[s] + b.fruit*fruit[s] + b.hfi*hfi[s] + b.evi*evi[s] + b.wat*distWat[s] + b.ppt*ppt[s] 
    gprobs[s] <- lam[s]/sum(lam[1:S])
    }
    
    ### likelihood
    for(i in 1:M){
    g[i] ~ dcat(gprobs[])
    z[i] ~ dbern(psi)
    logit(p[i]) <- p.alpha[focal[g[i]]] + p.site*site[g[i]] 
    y[i] ~ dbin(mu[i], K[i]) # K[i] allows to account for different number of segments
    mu[i] <- z[i]*p[i]
    e[i]<- (y[i]^0.5-p[i]*K[i])^2
    resid[i]<- (y[i]- p[i]*K[i])*z[i]
 
    y.new[i] ~ dbin(mu[i], K[i])
    g.new[i] ~ dcat(gprobs[])
    z.new[i] ~ dbern(psi)
    e.new[i] <- ((y.new[i])^0.5- p[i]*K[i])^2  # Freeman Tukey 
    resid.new[i]<- (y.new[i]- p[i]*K[i])*z[i]
    }

    #### derived quantities
    e.sum<- sum(e[])
    e.sum.new<- sum(e.new[])

    N <- sum(z[1:M])
    N.new<- sum(z.new[1:M])

    for(i in 1:M){
      group.out[i] <- g[i]*z[i]
      group.out.new[i]<- g.new[i]*z.new[i]
  
      #This will allow us to count the number of individuals in each replicate
      for(j in 1:S){
        g.N[j,i] <- step(0.01*(j-group.out[i])-0.02*(j-group.out[i])*(j-group.out[i])+0.001)
        g.N.new[j,i] <- step(0.01*(j-group.out.new[i])-0.02*(j-group.out.new[i])*(j-group.out.new[i])+0.001)
        }
      }

    for(j in 1:S){
      G.N[j]<- sum(g.N[j,])
      G.N.new[j]<- sum(g.N.new[j,])
  
      ge[j]<- ((G.N[j] - exp(lam[j]))/exp(lam[j]))^2 # chi square test
      ge.new[j]<- ((G.N.new[j] - exp(lam[j]))/exp(lam[j]))^2 # chi square test
    }
    
    ge.sum<- sum(ge[])
    ge.sum.new<- sum(ge.new[])
    
    }
    ",fill = TRUE)
sink()

##### Set Parameters #####

inits <- function(){list(z = rep(1,M))}
parameters <- c("alpha", "b.sex", "b.fruit","b.evi", "b.hfi", "b.wat", "b.ppt", "psi",
                "p.alpha", "p.site", "mu.focal", "sigma.focal",
                "N", "N.new", "e.sum", "e.sum.new", "ge.sum", "ge.sum.new", "resid", "resid.new","y.new", "z.new", 'g.new', 'G.N', 'G.N.new')

ni <- 40000
nt <- 3 
nb <- 1000
nc <- 3
na <- 20000

##### Data #####
dat<- read.csv("./Data/IndivCaptEvents_wCovs.csv", header = T)%>%
  mutate_if(is.factor, as.character)%>%
  mutate(y = rowSums(.[4:20], na.rm = T))%>%
  mutate(K = round(nSeg))%>%
  separate(Follow, into = c("Focal", "Follow_No"), sep = "_", remove =F)%>%
  dplyr::select(Follow, Focal, Site, Sex, Fruit, evi, distWat, hfi, ppt, y, K)%>%
  mutate(Follow = as.factor(Follow))%>%
  mutate(Sex = gsub("Male", 0, Sex),
         Sex = gsub("Female", 1, Sex),
         Sex = as.numeric(Sex),
         Site = gsub("Ivindo", 0, Site),
         Site = gsub("Wonga Wongue", 1, Site),
         Site = as.numeric(Site))


dat$groupID<- as.numeric(dat$Follow)
dat$Fruit[dat$Fruit == 0]<- 0.001
dat$Fruit<- log(dat$Fruit)


##### Creating folds #####
k = 5
dat <- fold(dat, k = k,  id_col = 'Follow')%>%
  ungroup()
A = 1000
##### Making an empty list for results ######
cov.results<- c()
gs.results <- c()


##### Cross Validation #####
for (fold in 1:k){
  # Separating the dataset
  training_set <- dat[dat$.folds != fold,]
  training_setIDs<- unique(training_set$groupID)
  training_set$groupID<- as.numeric(as.factor(training_set$groupID))
  
  
  testing_set <- dat[dat$.folds == fold,]%>%
    ungroup()%>%
    dplyr::select(groupID, Focal, Site, Sex, Fruit, evi, distWat, hfi, ppt, y, K)
  
  ### Model 
  training_set<- training_set%>% 
    ungroup()%>%
    dplyr::select(groupID, Focal, Site, Sex, Fruit, evi, distWat, hfi, ppt, y, K)
  
  covs<- training_set%>%
    dplyr::select(groupID, Focal, Site, Sex, Fruit, evi, distWat, hfi, ppt)%>%
    mutate(Focal = as.numeric(as.factor(Focal)))%>%
    unique()
  
  n.focals<- max(covs$Focal)
  
  
  M <- 1000
  A <- M- nrow(training_set)
  
  training_set<-  rbind(training_set, data.frame(groupID = rep(NA, A), Focal = rep(NA, A), Site = rep(NA,A),
                                                 Sex = rep(NA, A),
                                                 Fruit = rep(NA,A), evi = rep(NA, A), distWat = rep(NA,A),  hfi = rep(NA,A), 
                                                 ppt = rep(NA, A), y = rep(0, A), K = rep(8, A)))
  
  jags.data <- list(y = training_set$y, g = training_set$groupID, S = nrow(covs), M = M, K = training_set$K, 
                    site = as.vector(covs$Site),
                    fruit = as.vector(covs$Fruit), 
                    sex = as.vector(covs$Sex),
                    hfi = as.vector(covs$hfi), 
                    evi = as.vector(covs$evi),
                    distWat = as.vector(covs$distWat),
                    ppt = as.vector(covs$ppt),
                    focal = as.vector(covs$Focal),
                    n.focals = n.focals)
  
  mod <- jags(jags.data, inits, parameters, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt = na) 
  
  train.res.covs<- data.frame(fold = fold, dataset = "training", as.data.frame(mod$mean[1:7]), deviance = mod$mean$deviance)
  Xtimes<- length(mod$q50$G.N)
  train_res.gs<- data.frame(fold = rep(fold, Xtimes), dataset = rep("training", Xtimes), 
                            groupID = training_setIDs, groupSize = mod$q50$G.N)
  
  ### Prediction
  testing_set.gs<- testing_set%>%
    ungroup()%>%
    mutate(groupSize = exp(train.res.covs$alpha + train.res.covs$b.sex*testing_set$Sex + 
                             train.res.covs$b.fruit*testing_set$Fruit + train.res.covs$b.hfi*testing_set$hfi +
                             train.res.covs$b.ppt*testing_set$ppt),
           fold = fold,
           dataset = "test")%>%
    dplyr::select(fold, dataset, groupID, groupSize)
  
  
  ### Dataframes to save
  cov.results<- rbind(cov.results, train.res.covs)
  gs<- rbind(train_res.gs, testing_set.gs)
  gs.results <- rbind(gs.results, gs)
  
  
  mod.it <- paste0("mod.",fold)
  assign(mod.it, mod)
  print(fold)
}


cov.results

gs2<- gs.results%>%
  unique()%>%
  group_by(dataset, groupID)%>%
  summarise(gs = mean(groupSize), .groups = "keep")%>%
  ungroup()%>%
  pivot_wider(id_cols = "groupID", names_from = "dataset", values_from= "gs")

gs2<- gs2%>%
  pivot_longer(cols = c("test", "training"), names_to = "dataset")

ggplot(gs2, aes(x = value, fill = dataset))+
  geom_density(alpha = 0.4)

ggplot(old, aes(x = groupSize, fill = dataset))+
  geom_density(alpha = 0.4)


ggplot(gs2, aes(x = dataset, y = value))+
  geom_boxplot()


x<- list(Test = gs2$test, Training = gs2$training)
out<- overlap(x, plot = T)
out