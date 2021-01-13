
setwd("/Applications/Books/research/paper/Di et al.(2009)/Report/codes")
library(refund)
library(MASS)
library(mgcv)
library(simex)
library(splines)

source("GeneData.R")
source("fbps.cov.R")
source("mfpca.sc2.R")

################################
# DTI data
################################
data(DTI)
DTI = subset(DTI, Nscans < 6)  
id  = DTI$ID
Y = DTI$cca
D = nrow(Y)

mfpca.DTI =  mfpca.sc(Y=Y, id = id, twoway=FALSE)
mfpca.DTI$evalues

set.seed(1)
new_idx = sample(1:D)
Y_new = Y[new_idx,]
id_new = id[new_idx]
mfpca.DTI2 =  mfpca.sc(Y=Y_new, id = id_new, twoway=FALSE)
mfpca.DTI2$evalues


#### 1. mfpca.cs in refund
s1 <- Sys.time()
fit1 <-  mfpca.sc(Y=Y, id=id, twoway=TRUE)
s2 <- Sys.time()
time1 <- difftime(s2, s1, units = "mins")
# MISE of observations
diff1=0;num=0
for(i in 1:nrow(Y)){
  idx = which(!is.na(Y[i, ]))
  num = num + length(idx)
  diff1 = diff1 + sum(abs(fit1$Yhat[i,idx]-Y[i,idx])^2)
}
MISE1_Y <- diff1/num
mpca1 <- c(time1, MISE1_Y)

#### 2. the modified mfpca function: mfpca.cs2 
s1 <- Sys.time()
fit2 <-  mfpca.sc2(Y=Y, id=id, twoway=TRUE, design="irregular")
s2 <- Sys.time()
time2 <- difftime(s2, s1, units = "mins")
# MISE of observations
diff2=0; num=0
for(i in 1:nrow(Y)){
  idx = which(!is.na(Y[i, ]))
  num = num + length(idx)
  diff2 = diff2 + sum(abs(fit2$Yhat[i,idx]-Y[i,idx])^2)
}
MISE2_Y <- diff2/num
mpca2 <- c(time2, MISE2_Y)

round(mpca1,4)
round(mpca2,4)



################################
# Simulation data
################################
Nsub=100; J=3; D=500
design="irregular"
data <- GeneData(M=Nsub, J=J, N=D,  design=design, level=0.1, sigma=0)
Y <- data$Y

# True values
evalues_true <- data$evalues
eigenf_true <- data$eigenfunctions

# Parameters
K1 <- 4
K2 <- 4
id <- rep(1:Nsub, each=J)
twoway <- TRUE


#########################################################
####step 1: mfpca.cs in refund
#########################################################
s1 <- Sys.time()
fit1 <-  mfpca.sc(Y=Y, id=id, twoway=twoway)
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
})))/(K1*D)
MISE1_eigen2 <- sum(unlist(lapply(1:K2, function(x){
  min(sum((eigenf_true[[2]][,x]-fit1$efunctions[[2]][,x])^2),sum((eigenf_true[[2]][,x]+fit1$efunctions[[2]][,x])^2))
})))/(K2*D)
mpca1 <- c(time1, MISE1_Y, MISE1_eigen1, MISE1_eigen2)


#########################################################
####step 2: the revised mfpca function: mfpca.cs2 
#########################################################
s1 <- Sys.time()
fit2 <-  mfpca.sc2(Y=Y, id=id, twoway=twoway, design=design)
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
})))/(K1*D)
MISE2_eigen2 <- sum(unlist(lapply(1:K2, function(x){
  min(sum((eigenf_true[[2]][,x]-fit2$efunctions[[2]][,x])^2),sum((eigenf_true[[2]][,x]+fit2$efunctions[[2]][,x])^2))
})))/(K2*D)
mpca2 <- c(time2, MISE2_Y, MISE2_eigen1, MISE2_eigen2)


round(mpca1,4)
round(mpca2,4)


