
setwd("~/Users/new_mfpca")
library(refund)
library(MASS)
library(Matrix)
library(mgcv)
library(simex)
library(splines)
library(rARPACK)

source("GeneData.2.R")
source("codes/fbps.cov.R")
source("codes/backup.R")
source("mfpca.fast.R")
source("mfpca.fast2.R")
source("mfpca.fast3.R")



set.seed(3)
Nsub=100; J=5; L=5000
design="regular"
data <- GeneData(I=Nsub, J=J, L=L,  design=design, level=0.95, sigma=1, balanced=T)
Y <- data$Y
visit <- rep(1:J,times=Nsub)
id <- rep(1:Nsub,each=J)


# True values
evalues_true <- data$evalues
eigenf_true <- data$eigenfunctions
# Parameters
K1 <- 4
K2 <- 4


#########################################################
####step 1: mfpca.fast
#########################################################
s1 <- Sys.time()
fit1 <- mfpca.fast(Y=Y, id=id)
s2 <- Sys.time()
time1 <- difftime(s2, s1, units="mins")


# MISE of observations
diff1 = 0
num = 0
for(i in 1:nrow(Y)){
  idx = which(!is.na(Y[i, ]))
  num = num + length(idx)
  diff1 = diff1 + sum(abs(fit1$Yhat[i,idx]-Y[i,idx])^2)
}
MISE1_Y <- diff1/num

# MISE of eigenfucntions
MISE1_eigen1 <- sum(unlist(lapply(1:K1, function(x){
  min(sum((eigenf_true[[1]][,x]-fit1$efunctions[[1]][,x])^2),sum((eigenf_true[[1]][,x]+fit1$efunctions[[1]][,x])^2))
})))/(K1*L)
MISE1_eigen2 <- sum(unlist(lapply(1:K2, function(x){
  min(sum((eigenf_true[[2]][,x]-fit1$efunctions[[2]][,x])^2),sum((eigenf_true[[2]][,x]+fit1$efunctions[[2]][,x])^2))
})))/(K2*L)
mpca1 <- c(time1, MISE1_Y, MISE1_eigen1, MISE1_eigen2)




#########################################################
####step 2:  mfpca.fast2
#########################################################
s1 <- Sys.time()
fit2 <-   mfpca.fast2(Y=Y, id=id)
s2 <- Sys.time()
time2 <- difftime(s2, s1, units="mins")


# MISE of observations
diff2 = 0
num = 0
for(i in 1:nrow(Y)){
  idx = which(!is.na(Y[i, ]))
  num = num + length(idx)
  diff2 = diff2 + sum(abs(fit2$Yhat[i,idx]-Y[i,idx])^2)
}
MISE2_Y <- diff2/num

# MISE of eigenfucntions
MISE2_eigen1 <- sum(unlist(lapply(1:K1, function(x){
  min(sum((eigenf_true[[1]][,x]-fit2$efunctions[[1]][,x])^2),sum((eigenf_true[[1]][,x]+fit2$efunctions[[1]][,x])^2))
})))/(K1*L)
MISE2_eigen2 <- sum(unlist(lapply(1:K2, function(x){
  min(sum((eigenf_true[[2]][,x]-fit2$efunctions[[2]][,x])^2),sum((eigenf_true[[2]][,x]+fit2$efunctions[[2]][,x])^2))
})))/(K2*L)
mpca2 <- c(time2, MISE2_Y, MISE2_eigen1, MISE2_eigen2)




#########################################################
####step 3:  mfpca.fast3
#########################################################
s1 <- Sys.time()
fit3 <-   mfpca.fast3(Y=Y, id=id)
s2 <- Sys.time()
time3 <- difftime(s2, s1, units="mins")


# MISE of observations
diff3 = 0
num = 0
for(i in 1:nrow(Y)){
  idx = which(!is.na(Y[i, ]))
  num = num + length(idx)
  diff3 = diff3 + sum(abs(fit2$Yhat[i,idx]-Y[i,idx])^2)
}
MISE3_Y <- diff3/num

# MISE of eigenfucntions
MISE3_eigen1 <- sum(unlist(lapply(1:K1, function(x){
  min(sum((eigenf_true[[1]][,x]-fit3$efunctions[[1]][,x])^2),sum((eigenf_true[[1]][,x]+fit3$efunctions[[1]][,x])^2))
})))/(K1*L)
MISE3_eigen2 <- sum(unlist(lapply(1:K2, function(x){
  min(sum((eigenf_true[[2]][,x]-fit3$efunctions[[2]][,x])^2),sum((eigenf_true[[2]][,x]+fit3$efunctions[[2]][,x])^2))
})))/(K2*L)
mpca3 <- c(time3, MISE3_Y, MISE3_eigen1, MISE3_eigen2)




#########################################################
####step 3: mfpca.cs in refund
#########################################################
# s1 <- Sys.time()
# fit3 <-  mfpca.sc(Y=Y, id=id, twoway=TRUE)
# s2 <- Sys.time()
# time3 <- difftime(s2, s1, units="mins")
# 
# 
# # MISE of observations
# diff3 = 0
# num = 0
# for(i in 1:nrow(Y)){
#   idx = which(!is.na(Y[i, ]))
#   num = num + length(idx)
#   diff3 = diff3 + sum(abs(fit3$Yhat[i,idx]-Y[i,idx])^2)
# }
# MISE3_Y <- diff3/num
# 
# # MISE of eigenfucntions
# MISE3_eigen1 <- sum(unlist(lapply(1:K1, function(x){
#   min(sum((eigenf_true[[1]][,x]-fit3$efunctions[[1]][,x])^2),sum((eigenf_true[[1]][,x]+fit3$efunctions[[1]][,x])^2))
# })))/(K1*L)
# MISE3_eigen2 <- sum(unlist(lapply(1:K2, function(x){
#   min(sum((eigenf_true[[2]][,x]-fit3$efunctions[[2]][,x])^2),sum((eigenf_true[[2]][,x]+fit3$efunctions[[2]][,x])^2))
# })))/(K2*L)
# mpca3 <- c(time3, MISE3_Y, MISE3_eigen1, MISE3_eigen2)




################################################################
# MISE of eigenfunctions bwtween mfpca.fast and mfpca.fast2
################################################################
MISE_eigen1 <- sum(unlist(lapply(1:K1, function(x){
  min(sum((fit1$efunctions[[1]][,x]-fit2$efunctions[[1]][,x])^2),sum((fit1$efunctions[[1]][,x]+fit2$efunctions[[1]][,x])^2))
})))/(K1*L)
MISE_eigen2 <- sum(unlist(lapply(1:K2, function(x){
  min(sum((fit1$efunctions[[2]][,x]-fit2$efunctions[[2]][,x])^2),sum((fit1$efunctions[[2]][,x]+fit2$efunctions[[2]][,x])^2))
})))/(K2*L)

mpca <- c(MISE_eigen1, MISE_eigen2)



round(mpca1,4)
round(mpca2,4)
round(mpca3,4)
round(mpca,4)













