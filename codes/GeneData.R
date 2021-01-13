
################################################################################
# M: the number of subjects.
# J: the number of visit for each subject.
# N: the number of observations in each curve.
# design: generate regular or irregular spaced data, default as "regular". 
# level: the level of sparse
# sigma: the standard deviation of random errors
################################################################################
GeneData <- function(M=100, J=3, N=100, design="regular", level=0.1, sigma=1){
  
  K1 <- 4
  K2 <- 4
  K <- K1 + K2
  lambda1 <- 0.5^(0:(K1-1))
  lambda2 <- 0.5^(0:(K2-1))
  tlength <- 1
  t <- seq(0,tlength,length=N)
  
  # Eigenfunctions
  f1 <- matrix(0,nrow=K1,ncol=N)
  for ( i in 1:(K1/2) ) {
    f1[2*i-1,] <- sqrt(2/tlength)*sin(i*t*2*pi/tlength)
    f1[2*i,] <- sqrt(2/tlength)*cos(i*t*2*pi*tlength)
  }
  tt <- t/tlength
  f2 <- matrix(0, nrow=K2, ncol=N)
  f2[1,] <- rep(1, N)*sqrt(1/tlength)
  f2[2,] <- sqrt(3/tlength) * (2*tt - 1)
  f2[3,] <- sqrt(5/tlength) * (6*tt^2 - 6 * tt + 1)
  f2[4,] <- sqrt(7/tlength) * (20*tt^3 - 30*tt^2 + 12 * tt -1)
  
  
  # Generate scores
  si1 <- matrix(0, nrow=M, ncol=K1)
  si2 <- matrix(0, nrow=M*J, ncol=K2)
  for(k in 1:K1) {
    si1[,k] <- rnorm(M, sd=sqrt(lambda1[k]))
  }
  for(k in 1:K2) {
    si2[,k] <- rnorm(M*J, sd=sqrt(lambda2[k]))
  }
  # Generate errors
  epsilon <- matrix(rnorm(M*J*N,sd=sigma),nc=N)
  
  
  # Generate dense data
  Y0 <- matrix(0,nrow=M*J,ncol=N)
  for(m in 1:M) {
    temp <- apply( ( si1[m,] %*% t(rep(1,N)) ) * f1, 2,sum)
    for(j in 1:J) {
      Y0[(m-1)*J+j ,] <- temp + apply( ( si2[(m-1)*J + j ,] %*% t(rep(1,N)) ) * f2, 2,sum) + epsilon[(m-1)*J+j,]
    }
  }
  
  # Generate sparse data
  if (design == "regular") {
    Y <- Y0
  } else {
    nobs <- floor(level*N)
    Y <- matrix(NA,nrow=M*J,ncol=N)
    for (i in 1:(M*J)) {
      idx <- sample(1:N, nobs)
      Y[i,idx] <- Y0[i,idx]
    }
  }
  
  # return values
  evalues <- list(level1=lambda1, level2=lambda2)
  eigenfunctions <- list(level1=t(f1), level2=t(f2))
  
  return(list(Y=Y, evalues=evalues, eigenfunctions=eigenfunctions))
}



