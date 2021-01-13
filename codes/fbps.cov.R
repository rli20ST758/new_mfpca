recov <-function(Sig,ratio){
  
  ## bin a covariance matrix while ignoring its diagonal elements
  
  ## Sig: an estimated (or true) covariance matrix 
  ## ratio: the bin size
  
  m=dim(Sig)[1]
  diag(Sig)=rep(NA,m)
  m1=sqrt(length(Sig))/ratio
  bin_cov=matrix(0,m1,m1)
  
  if(ratio==1){
    bin_cov = Sig
    Diag = rep(0,m)
    for(i in 2:(m-1))
      Diag[i] = (bin_cov[i-1,i]+bin_cov[i+1,i]+bin_cov[i,i-1]+bin_cov[i,i+1])/4
    Diag[1] = (bin_cov[1,2]+bin_cov[2,1])/2
    Diag[m] = (bin_cov[m-1,m]+bin_cov[m,m-1])/2
    diag(bin_cov) = Diag
  }
  
  if(ratio>=2){
    
    for (i in 1:(m1-1))
      for (j in (i+1):m1)
        bin_cov[i,j]=mean(Sig[(i-1)*ratio+(1:ratio),(j-1)*ratio+(1:ratio)])
      
      bin_cov=bin_cov+t(bin_cov)
      
      for (i in 1:m1)
        bin_cov[i,i]=mean(Sig[(i-1)*ratio+(1:ratio),(i-1)*ratio+(1:ratio)], na.rm=1)
  }
  
  return(bin_cov)
}




diff <-function(m,K){
  
  # parameter  m: difference order
  # parameter  K: size  
  
  M = matrix(0,nrow=K-m,ncol=K)
  c = rep(0,m+1)
  
  for(i in 0:m)
    c[i+1] = (-1)^(i+1)*factorial(m)/(factorial(i)*factorial(m-i))
  
  for(i in 1:(K-m)) M[i,i:(i+m)] = c
  
  return(M)
}


pspline.setting <-function(x,knots=35,p=3,m=2){
  
  # x: the marginal data points
  # K: the list of knots or the numbers of knots
  # p: degrees for B-splines, with defaults values 3
  # m: orders of difference penalty, with default values 2
  #library(splines)
  if(length(knots)==1)
  {
    K = knots
    knots=seq(-p,K+p,length=K+1+2*p)/K
    knots = knots*(max(x)-min(x)) + min(x)
  }
  
  if(length(knots)>1) 
  { knots = knots
  K = length(knots)-2*p-1
  }
  
  
  P = diff(m,K+p)
  P = t(P)%*%P
  
  ### design matrix and some pre-calculation 
  ### for the penalty without interaction
  ### The idea about pre-calculation, when lambda
  ## is changed, only eigenvalues change.
  
  B = spline.des(knots=knots, x=x, ord = p+1, outer.ok = TRUE)$design
  
  Sig = t(B)%*%B
  eSig = eigen(Sig)
  V = eSig$vectors
  E = eSig$values
  
  Sigi_sqrt =V%*%diag(1/sqrt(E))%*%t(V)
  Sigi = V%*%diag(1/E)%*%t(V)
  
  tUPU = t(Sigi_sqrt)%*%P%*%Sigi_sqrt
  
  Esig = eigen(tUPU)
  U = Esig$vectors
  s = Esig$values
  s[(K+p-m+1):(K+p)]=0
  A = B%*%Sigi_sqrt%*%U
  
  List = list(
    "A" = A,
    "s" = s,
    "Sigi.sqrt"=Sigi_sqrt,
    "U" = U)
  
  return(List)
}




fbps.cov <-function(data,x = NULL, diag.remove = TRUE, knots=35, p=3, m=2, lambda=NULL, method="L-BFGS-B",lower= -20, upper=20, control=NULL){

# return a smooth and symmetric covariance matrix using fbps

# data: a symmetric matrix 
# x: the collection of data points for each dimension for the covariance matrix
# diag.remove: if diagonal elements should be ignored
# knots: to specify either the number of  knots  or the vector of knots for each dimension; defaults to 35
# p: the degrees of B-splines
# m: the order of difference penalty
# lambda: the user-selected smoothing parameter  
# method: see "optim"
# lower, upper, control: see "optim"
  
#library(splines)
#source("pspline.setting.R")

## make sure the data is symmetric
data = (data+t(data))/2
n = dim(data)[1] 

if(is.null(x)==TRUE) x=(1:n)/n-1/2/n ## if NULL, assume equally distributed 

#######################################################################################

Y = data 

if(diag.remove==TRUE){
 ## use nearest neighboor approach to fill in the diagonals
  Y = recov(data,1)   # ignore the diagnoal elements and bin the raw covariance 
  }

###################   precalculation for bps smoothing  ##########################################66
if(length(knots)==1)
{
  K = knots
  knots=seq(-p,K+p,length=K+1+2*p)/K
  knots = knots*(max(x)-min(x)) + min(x)
  }

if(length(knots)>1) 
  { knots = knots
    K = length(knots)-2*p-1
}


List = pspline.setting(x,knots,p,m)
A = List$A
s = List$s
Sigi_sqrt = List$Sigi.sqrt
U = List$U


##B is the B-spline design matrix for the raw covariance matrix
B = spline.des(knots=knots, x=x, ord = p+1, outer.ok=TRUE)$design
#################select optimal penalty ####################################
tr <-function(A){ return(sum(diag(A)))} ## the trace of a square matrix

Ytilde = t(A)%*%Y%*%A
Y_sum = tr(t(Y)%*%Y)
## redefine fbps_cov_gcv and fbps_cov_est as functions of 
## the smoothing parameter

fbps_cov_gcv =function(x){
lambda=exp(x)
Sigma = diag(1/(1+lambda*s))
Sigma1 = diag(1/(1+lambda*s)^(1/2))
Ysigma = Sigma%*%Ytilde%*%Sigma
Ysigma1 = Sigma1%*%Ytilde%*%Sigma1
gcv = tr(Ysigma%*%Ysigma)- 2*tr(Ysigma1%*%Ysigma1)+Y_sum
trace = sum(1/(1+lambda*s))^2
gcv = gcv/(1-trace/(n^2))^2
return(gcv)
}

fbps_cov_est =function(x){
  
lambda=exp(x)
Sigma = diag(1/(1+lambda*s))
Sigma1 = diag(1/(1+lambda*s)^(1/2))
Ysigma = Sigma%*%Ytilde%*%Sigma
Ysigma1 = Sigma1%*%Ytilde%*%Sigma1
gcv = tr(Ysigma%*%Ysigma)- 2*tr(Ysigma1%*%Ysigma1)+Y_sum
trace = sum(1/(1+lambda*s))^2
gcv = gcv/(1-trace/(n^2))^2

Temp = Sigi_sqrt%*%U
Theta = Temp%*%Ysigma%*%t(Temp)
hatY = B%*%Theta%*%t(B)

result=list(lambda=lambda,hatY=hatY,trace=trace,gcv=gcv,Theta=Theta)
return(result)
}

if(is.null(lambda)==TRUE){
fit = optim(0,fbps_cov_gcv,method=method,control=control,
  lower=lower,upper=upper)
if(fit$convergence>0) {
  expression = paste("Smoothing failed! The code is:",fit$convergence)
  print(expression)
}
lambda = exp(fit$par)
}
est = fbps_cov_est(log(lambda))
fbps_cov=as.matrix(B%*%(est$Theta)%*%t(B))

ft=diag(data-fbps_cov)  
var.noise=sum(ft[ft>0])/length(ft)   # calculate the variance of noise term 
if(diag.remove==F) var.noise = NULL

result=list(cov=fbps_cov, var=var.noise)
return(result)
}
