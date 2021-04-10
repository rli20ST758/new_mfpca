## Compare the computing time between mfpca.fast() and mfpca.sc() under different simulation settings

rm(list=ls())
## load packages
library(refund)
library(MASS)
library(splines)
library(simex)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(tidyverse)
source("./code/GeneData.2.R")
source("./code/fbps.cov.R")
source("./code/mfpca.fast.R")
source("./code/mfpca.sc.R")


########################################################################
## generate multilevel functional data
########################################################################
I <- 100
J <- 2
L <- 100
K1 <- 4
K2 <- 4
design <- "regular"
sigma <- 1

## create a data frame storing simulation results
nsim <- 1
sim_res <- matrix(NA, nrow = nsim*2, ncol = 9)
colnames(sim_res) <- c("I", "J", "L", "iteration", 
                       "comptime", "MISE(Y)", "MISE(Phi)", "MISE(Psi)", "method")
## start simulation
ind <- 1
for(iter in 1:nsim){
  set.seed(iter)
  data <- GeneData(I = I, J = J, L = L, design = design, sigma = sigma, balanced = FALSE)
  Y <- data$Y
  ## true eigenvalues and eigenfunctions
  evalues_true <- data$evalues 
  eigenf_true <- data$eigenfunctions
  ## other parameters for estimation
  id <- data$id
  
  ## fit MFPCA using mfpca.fast()
  ptm <- proc.time()
  fit_fast <- mfpca.fast(Y = data$Y, id = id)
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
  fit_sc <- mfpca.sc0(Y = data$Y, id = id, twoway = TRUE)
  time_sc <- proc.time() - ptm
  # MISE of observations
  diff1 <- 0
  num <- 0
  for(i in 1:nrow(Y)){
    idx <- which(!is.na(Y[i, ]))
    num <- num + length(idx)
    diff1 <- diff1 + sum(abs(fit_sc$Yhat[i,idx]-Y[i,idx])^2)
  }
  MISE1_Y <- diff1/num
  # MISE of eigenfucntions
  MISE1_eigen1 <- sum(unlist(lapply(1:K1, function(x){
    min(sum((eigenf_true[[1]][,x]-fit_sc$efunctions[[1]][,x])^2),
        sum((eigenf_true[[1]][,x]+fit_sc$efunctions[[1]][,x])^2))})))/(K1*L)
  MISE1_eigen2 <- sum(unlist(lapply(1:K2, function(x){
    min(sum((eigenf_true[[2]][,x]-fit_sc$efunctions[[2]][,x])^2),
        sum((eigenf_true[[2]][,x]+fit_sc$efunctions[[2]][,x])^2))})))/(K2*L)
  sim_res[ind,1:8] <- round(c(I, J, L, iter, time_sc[3], MISE1_Y, MISE1_eigen1, MISE1_eigen2),4)
  sim_res[ind,9] <- "mfpca.sc"
  ind <- ind + 1
  
  print(iter)
  
}
sim_res <- as.data.frame(sim_res)
sim_res[,1:8] <- lapply(sim_res[,1:8],as.numeric)

## summarise the results of nsim iterations
sim_res_lite <- sim_res %>%
  group_by(I, J, L, method) %>%
  summarise(comptime = median(comptime), `MISE(Y)` = median(`MISE(Y)`),
            `MISE(Phi)` = median(`MISE(Phi)`), `MISE(Psi)` = median(`MISE(Psi)`))

# ## "Smoothing failed! The code is: 52" a lot when I = 100
# 
# ## organize the results
# sim_new <- as.data.frame(sim_new)
# sim_new$method <- "fmfpca"
# 
# 
# ## create a data frame storing simulation results in mfpca.sc()
# # values <- c(100, 200, 500) ## I
# # values <- c(100, 200) ## L
# values <- c(3, 20, 50) ## J
# nsim <- 50
# sim_old <- matrix(NA, nrow = length(values)*nsim, ncol = 8)
# colnames(sim_old) <- c("I", "J", "L", "iteration", 
#                        "comptime", "MISE(Y)", "MISE(Phi)", "MISE(Psi)")
# ## start simulation
# ind <- 1
# for(J in values){
#   for(iter in 1:nsim){
#     data <- GeneData(M = I, J = J, N = L, design = "regular", sigma = sigma)
#     Y <- data$Y
#     ## true eigenvalues and eigenfunctions
#     evalues_true <- data$evalues 
#     eigenf_true <- data$eigenfunctions
#     ## other parameters for estimation
#     id <- rep(1:I, each=J)
#     ## fit MFPCA using mfpca.sc()
#     ptm <- proc.time()
#     fit1 <- mfpca.sc(Y=data$Y, id=id, twoway=TRUE)
#     time1 <- proc.time() - ptm
#     # MISE of observations
#     diff1 <- 0
#     num <- 0
#     for(i in 1:nrow(Y)){
#       idx <- which(!is.na(Y[i, ]))
#       num <- num + length(idx)
#       diff1 <- diff1 + sum(abs(fit1$Yhat[i,idx]-Y[i,idx])^2)
#     }
#     MISE1_Y <- diff1/num
#     # MISE of eigenfucntions
#     MISE1_eigen1 <- sum(unlist(lapply(1:K1, function(x){
#       min(sum((eigenf_true[[1]][,x]-fit1$efunctions[[1]][,x])^2),
#           sum((eigenf_true[[1]][,x]+fit1$efunctions[[1]][,x])^2))})))/(K1*L)
#     MISE1_eigen2 <- sum(unlist(lapply(1:K2, function(x){
#       min(sum((eigenf_true[[2]][,x]-fit1$efunctions[[2]][,x])^2),
#           sum((eigenf_true[[2]][,x]+fit1$efunctions[[2]][,x])^2))})))/(K2*L)
#     sim_old[ind,] <- round(c(I, J, L, iter, time1[3], MISE1_Y, MISE1_eigen1, MISE1_eigen2),4)
#     ind <- ind + 1
#   }
#   print(paste0("value = ", J))
# }
# sim_old <- as.data.frame(sim_old)
# sim_old$method <- "mfpca.sc"
# 
# sim_res <- rbind(sim_new, sim_old)
# saveRDS(sim_res, file = "./results/regular_I100_J_L100.rds")
# 
# 
# 
# ## plot the results under different I
# sim_res_I <- readRDS("./results/regular_I_J3_L100.rds")
# sim_res_simp_I <- sim_res_I %>%
#   group_by(I, J, L, method) %>%
#   summarise(comptime = median(comptime), `MISE(Y)` = median(`MISE(Y)`),
#             `MISE(Phi)` = median(`MISE(Phi)`), `MISE(Psi)` = median(`MISE(Psi)`))
# 
# plist <- list()
# plist[["I_Time"]] <- ggplot(sim_res_simp_I) +
#   theme_classic() +
#   geom_line(aes(x = I, y = comptime, color = method, linetype = method)) +
#   geom_point(aes(x = I, y = comptime, color = method)) +
#   scale_x_continuous(breaks = c(100,500,1000,1500,2000,2500,3000)) +
#   labs(x = "Number of subjects (I)", y = "Computing time (s)") +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# plist[["I_MISE(Phi)"]] <- ggplot(sim_res_simp_I) +
#   theme_classic() +
#   geom_line(aes(x = I, y = `MISE(Phi)`, color = method, linetype = method)) +
#   geom_point(aes(x = I, y = `MISE(Phi)`, color = method)) +
#   ylim(0, 0.08) +
#   scale_x_continuous(breaks = c(100,500,1000,1500,2000,2500,3000)) +
#   labs(x = "Number of subjects (I)", y = expression(paste("MISE(", Phi, ")"))) +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# plist[["I_MISE(Psi)"]] <- ggplot(sim_res_simp_I) +
#   theme_classic() +
#   geom_line(aes(x = I, y = `MISE(Psi)`, color = method, linetype = method)) +
#   geom_point(aes(x = I, y = `MISE(Psi)`, color = method)) +
#   ylim(0, 0.03) +
#   scale_x_continuous(breaks = c(100,500,1000,1500,2000,2500,3000)) +
#   labs(x = "Number of subjects (I)", y = expression(paste("MISE(", Psi, ")"))) +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# # pdf("./results/regular_I_J3_L100.pdf", width = 9, height = 3)
# # do.call(grid.arrange, c(plist, nrow = 1))
# # dev.off()
# 
# 
# ## plot the results under different J
# sim_res_J <- readRDS("./results/regular_I100_J_L100.rds")
# sim_res_simp_J <- sim_res_J %>%
#   group_by(I, J, L, method) %>%
#   summarise(comptime = median(comptime), `MISE(Y)` = median(`MISE(Y)`),
#             `MISE(Phi)` = median(`MISE(Phi)`), `MISE(Psi)` = median(`MISE(Psi)`))
# 
# plist[["J_Time"]] <- ggplot(sim_res_simp_J) +
#   theme_classic() +
#   geom_line(aes(x = J, y = comptime, color = method, linetype = method)) +
#   geom_point(aes(x = J, y = comptime, color = method)) +
#   scale_x_continuous(breaks = c(3,20,50,100,150)) +
#   labs(x = "Number of visits per subject (J)", y = "Computing time (s)") +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# plist[["J_MISE(Phi)"]] <- ggplot(sim_res_simp_J) +
#   theme_classic() +
#   geom_line(aes(x = J, y = `MISE(Phi)`, color = method, linetype = method)) +
#   geom_point(aes(x = J, y = `MISE(Phi)`, color = method)) +
#   scale_x_continuous(breaks = c(3,20,50,100,150)) +
#   ylim(0, 0.08) +
#   labs(x = "Number of visits per subject (J)", y = expression(paste("MISE(", Phi, ")"))) +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# plist[["J_MISE(Psi)"]] <- ggplot(sim_res_simp_J) +
#   theme_classic() +
#   geom_line(aes(x = J, y = `MISE(Psi)`, color = method, linetype = method)) +
#   geom_point(aes(x = J, y = `MISE(Psi)`, color = method)) +
#   scale_x_continuous(breaks = c(3,20,50,100,150)) +
#   ylim(0, 0.03) +
#   labs(x = "Number of visits per subject (J)", y = expression(paste("MISE(", Psi, ")"))) +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# 
# 
# ## plot the results under different L
# sim_res_L <- readRDS("./results/regular_I100_J3_L.rds")
# sim_res_simp_L <- sim_res_L %>%
#   group_by(I, J, L, method) %>%
#   summarise(comptime = median(comptime), `MISE(Y)` = median(`MISE(Y)`),
#             `MISE(Phi)` = median(`MISE(Phi)`), `MISE(Psi)` = median(`MISE(Psi)`))
# 
# plist[["L_Time"]] <- ggplot(sim_res_simp_L) +
#   theme_classic() +
#   geom_line(aes(x = L, y = comptime, color = method, linetype = method)) +
#   geom_point(aes(x = L, y = comptime, color = method)) +
#   scale_x_continuous(breaks = c(100,500,1000,1500)) +
#   labs(x = "Dimension of the functional domain (L)", y = "Computing time (s)") +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# plist[["L_MISE(Phi)"]] <- ggplot(sim_res_simp_L) +
#   theme_classic() +
#   geom_line(aes(x = L, y = `MISE(Phi)`, color = method, linetype = method)) +
#   geom_point(aes(x = L, y = `MISE(Phi)`, color = method)) +
#   scale_x_continuous(breaks = c(100,500,1000,1500)) +
#   ylim(0, 0.08) +
#   labs(x = "Dimension of the functional domain (L)", y = expression(paste("MISE(", Phi, ")"))) +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# plist[["L_MISE(Psi)"]] <- ggplot(sim_res_simp_L) +
#   theme_classic() +
#   geom_line(aes(x = L, y = `MISE(Psi)`, color = method, linetype = method)) +
#   geom_point(aes(x = L, y = `MISE(Psi)`, color = method)) +
#   scale_x_continuous(breaks = c(100,500,1000,1500)) +
#   ylim(0, 0.03) +
#   labs(x = "Dimension of the functional domain (L)", y = expression(paste("MISE(", Psi, ")"))) +
#   theme(legend.position="top", axis.title = element_text(size = 10))
# 
# 
# pdf("./figures/sim_regular.pdf", width = 9, height = 9)
# do.call(grid.arrange, c(plist, nrow = 3))
# dev.off()
# 
# 
# ## single run code
# data <- GeneData(M = I, J = J, N = L, design = "regular", sigma = sigma)
# Y <- data$Y
# 
# ## true eigenvalues and eigenfunctions
# evalues_true <- data$evalues 
# eigenf_true <- data$eigenfunctions
# 
# ## other parameters for estimation
# id <- rep(1:I, each=J)
# 
# 
# ########################################################################
# ## fit MFPCA using different methods
# ########################################################################
# 
# # ## 1. mfpca.sc
# # ptm <- proc.time()
# # fit1 <- mfpca.sc(Y=data$Y, id=id, twoway=TRUE)
# # time1 <- proc.time() - ptm
# # # MISE of observations
# # diff1 <- 0
# # num <- 0
# # for(i in 1:nrow(Y)){
# #   idx <- which(!is.na(Y[i, ]))
# #   num <- num + length(idx)
# #   diff1 <- diff1 + sum(abs(fit1$Yhat[i,idx]-Y[i,idx])^2)
# # }
# # MISE1_Y <- diff1/num
# # # MISE of eigenfucntions
# # MISE1_eigen1 <- sum(unlist(lapply(1:K1, function(x){
# #   min(sum((eigenf_true[[1]][,x]-fit1$efunctions[[1]][,x])^2),
# #       sum((eigenf_true[[1]][,x]+fit1$efunctions[[1]][,x])^2))})))/(K1*L)
# # MISE1_eigen2 <- sum(unlist(lapply(1:K2, function(x){
# #   min(sum((eigenf_true[[2]][,x]-fit1$efunctions[[2]][,x])^2),
# #       sum((eigenf_true[[2]][,x]+fit1$efunctions[[2]][,x])^2))})))/(K2*L)
# # result1 <- round(c(time1[3], MISE1_Y, MISE1_eigen1, MISE1_eigen2),4)
# 
# 
# ## 2. fmfpca
# ptm <- proc.time()
# fit2 <- fmfpca(Y = data$Y, id = id)
# time2 <- proc.time() - ptm
# # print(time2)
# 
# ### I = 100, J = 3: L = 100 (0.3s), 500 (2.5s), 1000 (14s), 2000 (90s) slow in step 6
# ### J = 3, L = 100: I = 100 (0.3s), 500 (1.2s), 1000 (3.7s), 2000 (15s), 5000 (105s) slow in step 4
# ### I = 100, L = 100: J = 3 (0.3s), 20 (0.42s), 100 (1.5s), 200 (3s) slightly slow in data generation
# 
# # MISE of observations
# diff2 <- 0
# num <- 0
# for(i in 1:nrow(Y)){
#   idx = which(!is.na(Y[i, ]))
#   num = num + length(idx)
#   diff2 = diff2 + sum(abs(fit2$Xhat[i,idx]-Y[i,idx])^2)
# }
# MISE2_Y <- diff2/num
# # MISE of eigenfucntions
# MISE2_eigen1 <- sum(unlist(lapply(1:K1, function(x){
#   min(sum((eigenf_true[[1]][,x]-fit2$efunctions[[1]][,x])^2),
#       sum((eigenf_true[[1]][,x]+fit2$efunctions[[1]][,x])^2))})))/(K1*L)
# MISE2_eigen2 <- sum(unlist(lapply(1:K2, function(x){
#   min(sum((eigenf_true[[2]][,x]-fit2$efunctions[[2]][,x])^2),
#       sum((eigenf_true[[2]][,x]+fit2$efunctions[[2]][,x])^2))})))/(K2*L)
# result2 <- round(c(time2[3], MISE2_Y, MISE2_eigen1, MISE2_eigen2),4)
# 
# 
# # result <- rbind(result1,result2)
# # colnames(result) <- c("computing time", "MISE(Y)", "MISE(Phi)", "MISE(Psi)")
# # print("The result from mfpca.sc is")
# # print(result[1,])
# # print("The result from fmfpca is")
# # print(result[2,])
# 
