
rm(list=ls())
setwd("~/Users/new_mfpca/mfpca.face")
library(refund)
library(MASS)
library(Matrix)
library(mgcv)
library(simex)
library(splines)
library(rARPACK)
library(dplyr)

source("GeneData.R")
source("backup.R")
source("mfpca.face.R")

########################################################################
## generate multilevel functional data
########################################################################
I <- 200
J <- 2
L <- 1000
K1 <- 4
K2 <- 4
design <- "regular"
sigma <- 1

## create a data frame storing simulation results
nsim <- 20
#sim_res <- matrix(NA, nrow = nsim*2, ncol = 9)
sim_res <- data.frame(matrix(NA, nrow = nsim*2, ncol = 8), rep(NA,nsim*2))
colnames(sim_res) <- c("I", "J", "L", "iteration", 
                       "comptime", "MISE(Y)", "MISE(Phi)", "MISE(Psi)", "method")
## start simulation
ind <- 1
for(iter in 1:nsim){
  set.seed(iter)
  data <- GeneData(I = I, J = J, L = L, design = design, sigma = sigma, balanced = F, level = 0.4)
  Y <- data$Y
  ## true eigenvalues and eigenfunctions
  evalues_true <- data$evalues 
  eigenf_true <- data$eigenfunctions
  ## other parameters for estimation
  id <- data$id
  
  ## fit MFPCA using mfpca.fast()
  ptm <- proc.time()
  fit_fast <- mfpca.face(Y = data$Y, id = id, weight = "obs")
  #fit_fast <- mfpca.fast2(Y = data$Y, id = id)
  time_fast <- proc.time() - ptm
  # MISE of observations
  diff2 <- 0
  num <- 0
  for(i in 1:nrow(Y)){
    idx = which(!is.na(Y[i, ]))
    num = num + length(idx)
    diff2 = diff2 + sum(abs(fit_fast$Xhat[i,idx]-Y[i,idx])^2)
  }
  MISE2_Y <- diff2/num
  # MISE of eigenfucntions
  MISE2_eigen1 <- sum(unlist(lapply(1:K1, function(x){
    min(sum((eigenf_true[[1]][,x]-fit_fast$efunctions[[1]][,x])^2),
        sum((eigenf_true[[1]][,x]+fit_fast$efunctions[[1]][,x])^2))})))/(K1*L)
  MISE2_eigen2 <- sum(unlist(lapply(1:K2, function(x){
    min(sum((eigenf_true[[2]][,x]-fit_fast$efunctions[[2]][,x])^2),
        sum((eigenf_true[[2]][,x]+fit_fast$efunctions[[2]][,x])^2))})))/(K2*L)
  sim_res[ind,1:8] <- round(c(I, J, L, iter, time_fast[3], MISE2_Y, MISE2_eigen1, MISE2_eigen2),4)
  sim_res[ind,9] <- "mfpca.fast"
  ind <- ind + 1
  
  ## fit MFPCA using mfpca.sc()
  ptm <- proc.time()
  fit_fast2 <- mfpca.face(Y = data$Y, id = id, weight = "subj")
  time_sc <- proc.time() - ptm
  # MISE of observations
  diff1 <- 0
  num <- 0
  for(i in 1:nrow(Y)){
    idx <- which(!is.na(Y[i, ]))
    num <- num + length(idx)
    diff1 <- diff1 + sum(abs(fit_fast2$Xhat[i,idx]-Y[i,idx])^2)
  }
  MISE1_Y <- diff1/num
  # MISE of eigenfucntions
  MISE1_eigen1 <- sum(unlist(lapply(1:K1, function(x){
    min(sum((eigenf_true[[1]][,x]-fit_fast2$efunctions[[1]][,x])^2),
        sum((eigenf_true[[1]][,x]+fit_fast2$efunctions[[1]][,x])^2))})))/(K1*L)
  MISE1_eigen2 <- sum(unlist(lapply(1:K2, function(x){
    min(sum((eigenf_true[[2]][,x]-fit_fast2$efunctions[[2]][,x])^2),
        sum((eigenf_true[[2]][,x]+fit_fast2$efunctions[[2]][,x])^2))})))/(K2*L)
  sim_res[ind,1:8] <- round(c(I, J, L, iter, time_sc[3], MISE1_Y, MISE1_eigen1, MISE1_eigen2),4)
  sim_res[ind,9] <- "mfpca.fast2"
  ind <- ind + 1
  
  print(iter)
  
}

## summarise the results of nsim iterations
sim_res_lite <- sim_res %>% 
  group_by(I, J, L, method) %>% 
  summarise(comptime = median(comptime), `MISE(Y)` = median(`MISE(Y)`),
            `MISE(Phi)` = median(`MISE(Phi)`), `MISE(Psi)` = median(`MISE(Psi)`))
sim_res_lite
