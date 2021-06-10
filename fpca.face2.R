fpca.face2 <- function(Y=NULL, ydata=NULL, Y.pred=NULL, argvals=NULL, pve=0.99, npc=NULL,
           var=FALSE, simul=FALSE, sim.alpha=0.95,
           center=TRUE, knots=35, p=3, m=2, lambda=NULL, alpha=1,
           search.grid=TRUE, search.length=100,
           method="L-BFGS-B", lower=-20, upper=20, control=NULL){
    
    ## data: Y, I by J data matrix, functions on rows
    ## argvals:  vector of J
    ## knots: to specify either the number of knots or the vectors of knots;
    ##        defaults to 35;
    ## p: the degree of B-splines;
    ## m: the order of difference penalty
    ## lambda: user-selected smoothing parameter
    ## method: see R function "optim"
    ## lower, upper, control: see R function "optim"
    #require(Matrix)
    #source("pspline.setting.R")
    stopifnot(!is.null(Y))
    stopifnot(is.matrix(Y))
    data_dim <- dim(Y)
    I <- data_dim[1] ## number of subjects
    J <- data_dim[2] ## number of obs per function
    
    if(is.null(argvals))  argvals <- (1:J)/J-1/2/J ## if NULL, assume equally spaced
    
    meanX <- rep(0,J)
    if(center) {##center the functions
      meanX <- colMeans(Y, na.rm=TRUE)
      meanX <- smooth.spline(argvals,meanX,all.knots =TRUE)$y
      Y <- t(t(Y)- meanX)
    }
    
    ## specify the B-spline basis: knots
    p.p <- p
    m.p <- m
    if(length(knots)==1){
      if(knots+p.p>=J) cat("Too many knots!\n")
      stopifnot(knots+p.p<J)
      
      K.p <- knots
      knots <- seq(-p.p,K.p+p.p,length=K.p+1+2*p.p)/K.p
      knots <- knots*(max(argvals)-min(argvals)) + min(argvals)
    }
    if(length(knots)>1) K.p <- length(knots)-2*p.p-1
    if(K.p>=J) cat("Too many knots!\n")
    stopifnot(K.p <J)
    c.p <- K.p + p.p
    
    ######### precalculation for smoothing #############
    List <- pspline.setting(argvals,knots,p.p,m.p)
    B <- List$B
    Bt <- Matrix(t(as.matrix(B)))
    s <- List$s
    Sigi.sqrt <- List$Sigi.sqrt
    U <- List$U
    A0 <- Sigi.sqrt%*%U
    MM <- function(A,s,option=1){
      if(option==2)
        return(A*(s%*%t(rep(1,dim(A)[2]))))
      if(option==1)
        return(A*(rep(1,dim(A)[1])%*%t(s)))
    }
    G <- crossprod(B) / nrow(B)
    eig_G <- eigen(G, symmetric = T)
    G_half <- eig_G$vectors %*% diag(sqrt(eig_G$values)) %*% t(eig_G$vectors)
    G_invhalf <- eig_G$vectors %*% diag(1/sqrt(eig_G$values)) %*% t(eig_G$vectors)
    Bnew <- as.matrix(B %*% G_invhalf)
    Anew <- G_half %*% A0
    
    ######## precalculation for missing data ########
    imputation <- FALSE
    Niter.miss <- 1
    
    Index.miss <- is.na(Y)
    num.miss <- rowSums(is.na(Y))
    if(sum(Index.miss)>0){
      for(i in 1:I){
        if(num.miss[i]>0){
          y <- Y[i,]
          seq <- (1:J)[!is.na(y)]
          seq2 <-(1:J)[is.na(y)]
          t1 <- argvals[seq]
          t2 <- argvals[seq2]
          fit <- smooth.spline(t1,y[seq])
          temp <- predict(fit,t2,all.knots=TRUE)$y
          if(max(t2)>max(t1)) temp[t2>max(t1)] <- mean(y[seq])#y[seq[length(seq)]]
          if(min(t2)<min(t1)) temp[t2<min(t1)] <- mean(y[seq])#y[seq[1]]
          Y[i,seq2] <- temp
        }
      }
      imputation <- TRUE
      Niter.miss <- 100
    }
    convergence.vector <- rep(0,Niter.miss);
    iter.miss <- 1
    lambda.input <- lambda
    totalmiss <- mean(Index.miss)
    
    while(iter.miss <= Niter.miss&&convergence.vector[iter.miss]==0) {
      ###################################################
      ######## Transform the Data           #############
      ###################################################
      Ytilde <- as.matrix(t(A0)%*%as.matrix(Bt%*%t(Y)))
      C_diag <- rowSums(Ytilde^2)

      ###################################################
      ########  Select Smoothing Parameters #############
      ###################################################
      Y_square <- sum(Y^2)
      Ytilde_square <- sum(Ytilde^2)
      face_gcv <- function(x) {
        lambda <- exp(x)
        lambda_s <- (lambda*s)^2/(1 + lambda*s)^2
        gcv <- sum(C_diag*lambda_s) - Ytilde_square + Y_square
        trace <- sum(1/(1+lambda*s))
        gcv <- gcv/(1-alpha*trace/J/(1-totalmiss))^2
        return(gcv)
      }
      
      
      
      if(is.null(lambda.input) && iter.miss<=2) {
        if(!search.grid){
          fit <- optim(0,face_gcv,method=method,lower=lower,upper=upper,control=control)
          if(fit$convergence>0) {
            expression <- paste("Smoothing failed! The code is:",fit$convergence)
            print(expression)
          }
          
          lambda <- exp(fit$par)
        } else {
          Lambda <- seq(lower,upper,length=search.length)
          Length <- length(Lambda)
          Gcv <- rep(0,Length)
          for(i in 1:Length)
            Gcv[i] <- face_gcv(Lambda[i])
          i0 <- which.min(Gcv)
          lambda <- exp(Lambda[i0])
        }
      }
      
      YS <- MM(Ytilde,1/(1+lambda*s),2)
      
      ###################################################
      ####  Eigendecomposition of Smoothed Data #########
      ###################################################

      temp0 <- YS%*%t(YS)/I
      temp <- as.matrix(Anew%*%as.matrix(temp0%*%t(Anew)))
      Eigen <- eigen(temp,symmetric=TRUE)
      A = Eigen$vectors
      Phi = Bnew %*% A
      Sigma = Eigen$values

      if(iter.miss>1&&iter.miss< Niter.miss) {
        diff <- norm(YS-YS.temp,"F")/norm(YS,"F")
        if(diff <= 0.02)
          convergence.vector[iter.miss+1] <- 1
      }
      
      YS.temp <- YS
      iter.miss <- iter.miss + 1
      N <- min(I,c.p)
      d <- Sigma[1:N]
      d <- d[d>0]
      per <- cumsum(d)/sum(d)
      N <- ifelse (is.null(npc), min(which(per>pve)), min(npc, length(d)))
      
      #########################################
      #######     Principal  Scores   #########
      ########   data imputation      #########
      #########################################
      
      if(imputation) {
        Phi.N <- Phi[,1:N]
        A.N <- G_invhalf %*% A[,1:N]
        d <- Sigma[1:N]
        sigmahat2  <-  max(mean(Y[!Index.miss]^2)-sum(Sigma),0)
        #Xi1 = t(Y%*%Phi.N)
        if(N>1){
          Xi <- solve(t(Phi.N)%*%Phi.N + diag(sigmahat2/d)) %*% t(as.matrix(Y%*%B) %*% A.N)
        } else{
          Xi <- solve(t(Phi.N)%*%Phi.N + sigmahat2/d) %*% t(as.matrix(Y%*%B) %*% A.N)
        }
        Yhat <- t(Phi.N %*% Xi)
        Y <- Y*(1-Index.miss) + Yhat*Index.miss
        if(sum(is.na(Y))>0)
          print("error")
      }
      #if(iter.miss%%10==0) print(iter.miss)
    } ## end of while loop
    
    ### now calculate scores
    if(is.null(Y.pred)) {
      Y.pred = Y
    } else {
      Y.pred = t(t(as.matrix(Y.pred))-meanX)
      Index.miss = is.na(Y.pred)
      num.miss <- rowSums(is.na(Y.pred))
      I = nrow(Y.pred)
    }
    
    N <- ifelse (is.null(npc), min(which(per>pve)), npc)
    if (N>ncol(Phi)) {
      warning(paste0("The requested npc of ", npc,
                     " is greater than the maximum allowable value of ",
                     ncol(Phi), ". Using ", ncol(Phi), "."))
      N <- ncol(Phi)
    }
    npc <- N
    Phi.N <- Phi[,1:N]
    A.N <- G_invhalf %*% A[,1:N]
    d <- Sigma[1:N]
    sigmahat2 <- max(mean(Y[!Index.miss]^2) -sum(Sigma),0)
  
    Y.pred[Index.miss] = 0
    temp = t(as.matrix(Y.pred%*%B) %*% A.N)
    if(N>1) {
      Xi <- solve(t(Phi.N)%*%Phi.N + diag(sigmahat2/d)) %*% temp
    } else {
      Xi <- solve(t(Phi.N)%*%Phi.N + sigmahat2/d) %*% temp
    }
    
    
    # If the data are incomplete
    if(sum(num.miss) > 0){
      if(N>1){
        for (i in 1:I) {
          if(num.miss[i] > 0){
            Phi.obs <- Phi.N[!Index.miss[i,],]
            Xi[,i] <- solve(t(Phi.obs)%*%Phi.obs + diag(sigmahat2/d)) %*% temp[,i]
          }
        }
      } else{
        for (i in 1:I) {
          if(num.miss[i] > 0){
            Phi.obs <- Phi.N[!Index.miss[i,],]
            Xi[,i] <- solve(t(Phi.obs)%*%Phi.obs + sigmahat2/d) %*% temp[,i]
          }
        }
      }
    }
    

    efunctions = Phi.N
    evalues = Sigma[1:N]
    Yhat <- Phi.N %*% Xi
    Yhat <- t(Yhat + meanX)
    
    scores <- t(Xi)
    mu <- meanX
    
    ret.objects <- c("Yhat", "Y", "scores", "mu", "efunctions", "evalues", "npc")
    if(var) {
      sigma2 = sigmahat2
      VarMats = vector("list",I)
      diag.var = matrix(NA, nrow=I,ncol=J)
      crit.val = rep(0,I)
      for(i.subj in 1:I){
        if(N >1){
          temp = sigma2*efunctions%*%solve(t(efunctions)%*%efunctions + sigma2*diag(evalues))%*%t(efunctions)
        } else{
          temp = sigma2*efunctions%*%solve(t(efunctions)%*%efunctions + sigma2*evalues)%*%t(efunctions)
        }
        VarMats[[i.subj]] = temp
        diag.var[i.subj,] = diag(temp)
        if (simul & sigma2 != 0) {
          norm.samp = mvrnorm(2500, mu = rep(0, J), Sigma = VarMats[[i.subj]])/matrix(sqrt(diag(VarMats[[i.subj]])),
                                                                                      nrow = 2500, ncol = J, byrow = TRUE)
          crit.val[i.subj] = quantile(apply(abs(norm.samp), 1, max), sim.alpha)
        }
      }
      ret.objects = c(ret.objects,"sigma2","diag.var","VarMats")
      if (simul) {
        #require(MASS)
        ret.objects = c(ret.objects, "crit.val")
      }
    }
    ret = lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
    names(ret) = ret.objects
    class(ret) = "fpca"
    return(ret)
    
}




pspline.setting <- function(x,knots=select.knots(x,35),
                            p=3,m=2,weight=NULL,type="full",
                            knots.option="equally-spaced"){
  
  # x: the marginal data points
  # knots: the list of interior knots or the numbers of interior knots
  # p: degrees for B-splines, with defaults values 3
  # m: orders of difference penalty, with default values 2
  # knots.option: type of knots placement, with default values "equally-spaced"
  
  #require(splines)
  #require(Matrix)
  
  ### design matrix 
  K = length(knots)-2*p-1
  B = spline.des(knots=knots, x=x, ord = p+1,outer.ok = TRUE,sparse=TRUE)$design
  
  bs = "ps"
  
  if(knots.option == "quantile"){
    bs = "bs"
  }
  
  s.object = s(x=x, bs=bs, k=K+p,m=c(p-1,2), sp=NULL)
  object  = smooth.construct(s.object,data = data.frame(x=x),knots=list(x=knots))
  P =  object$S[[1]]
  if(knots.option == "quantile") P = P / max(abs(P))*10 # rescaling
  
  if(is.null(weight)) weight <- rep(1,length(x))
  
  if(type=="full"){
    
    Sig = crossprod(matrix.multiply(B,weight,option=2),B)
    eSig = eigen(Sig)
    V = eSig$vectors
    E = eSig$values
    if(min(E)<=0.0000001) {#cat("Warning! t(B)%*%B is singular!\n");
      #cat("A small identity matrix is added!\n");
      E <- E + 0.000001;
      
    }
    Sigi_sqrt = matrix.multiply(V,1/sqrt(E))%*%t(V)
    
    tUPU = Sigi_sqrt%*%(P%*%Sigi_sqrt)
    Esig = eigen(tUPU,symmetric=TRUE)
    U = Esig$vectors
    s = Esig$values
    s[(K+p-m+1):(K+p)]=0
    A = B%*%(Sigi_sqrt%*%U)
  }
  
  if(type=="simple"){
    
    A = NULL
    s = NULL
    Sigi_sqrt = NULL
    U = NULL
  }
  List = list(
    "A" = A,
    "B" = B,
    "s" = s,
    "Sigi.sqrt" = Sigi_sqrt,
    "U" = U,
    "P" = P)
  
  return(List)
}


matrix.multiply <- function(A,s,option=1){
  if(option==2)
    return(A*(s%*%t(rep(1,dim(A)[2]))))
  if(option==1)
    return(A*(rep(1,dim(A)[1])%*%t(s)))
}