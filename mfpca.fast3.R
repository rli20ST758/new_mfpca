#' This function performs multilevel FPCA (Di et al.2009) on multilevel functional data
#' It is based on the "mfpca.sc2()" function from Ruonan Li and Luo Xiao, 
#' and the "mfpca.sc()" function from Julia Wrobel in the R refund package.
#' 
#' @param Y a multilevel/longitudinal functional dataset on a regular grid stored in a matrix
#' @param id a vector containing the id information
#' @param group a vector containing information used to identify groups/visits
#' @param twoway logical, indicating whether to carry out twoway ANOVA and calculate visit-specific means. Defaults to \code{TRUE}.
#' @param weight: the way of calculating total covariance, "obs" indicating that the sample covariance is weighted by observations and "subj" indicating that the sample covariance is weighted equally by subjects. Defaults to \code{"obs"}.
#' @param smooth: the way of calculating within covariance, "separate" indicating that smooth two parts of within covariance separately and "joint" indicating that smooth within covariance jointly. Defaults to \code{"separate"}.
#' @param argvals a vector containing observed locations on the functional domain
#' @param pve proportion of variance explained: used to choose the number of principal components
#' @param npc prespecified value for the number of principal components (if given, this overrides \code{pve})
#' @param p integer; the degree of B-splines functions to use
#' @param m integer; the order of difference penalty to use
#' @param knots number of knots to use or the vectors of knots; defaults to 35
#' 
#' @export
#' @import refund
#' @import splines
#' @import mgcv
#' @import MASS
#' @import simex
#' @import Matrix

mfpca.fast3 <- function(Y, id, group = NULL, twoway = TRUE, weight = "obs", smooth = "separate",
                        argvals = NULL, pve = 0.99, npc = NULL, p = 3, m = 2, knots = 35, silent = TRUE){
  ## required packages
  library(refund)
  library(splines)
  library(mgcv)
  library(MASS)
  library(simex)
  library(Matrix)
  #source("./code/backup.R")
  
  ##################################################################################
  ## Organize the input
  ##################################################################################
  if(silent == FALSE) print("Step 0: Organize the input -- may take some time")
  
  stopifnot((!is.null(Y) & !is.null(id)))
  stopifnot(is.matrix(Y))
  
  ## specify group variable if not provided
  if (!is.null(group)){ 
    group <- as.factor(group)
  }else{ ## if group is not provided, assume the group id are 1,2,... for each subject
    group <- as.factor(ave(id, id, FUN=seq_along))
  }
  id <- as.factor(id) ## convert id into a factor
  
  
  ## organize data into one data frame
  df <- data.frame(id = id, group = group, Y = I(Y))
  rm(id, group, Y)
  
  ## derive several variables that will be used later
  J <- length(levels(df$group)) ## number of groups
  L <- ncol(df$Y) ## number of observations along the domain
  nGroups <- data.frame(table(df$id))  ## calculate number of groups for each subject
  colnames(nGroups) = c("id", "numGroups")
  ID = sort(unique(df$id)) ## id of each subject
  I <- length(ID) ## number of subjects
  ## assume observations are equally-spaced on [0,1] if not specified
  if (is.null(argvals))  argvals <- seq(0, 1, length.out=L) 
  
  ##################################################################################
  ## Estimate population mean function (mu)
  ##################################################################################
  if(silent == FALSE) print("Step 1: Estimate population mean function (mu)")
  
  meanY <- colMeans(df$Y, na.rm = TRUE)
  fit_mu <- gam(meanY ~ s(argvals))
  mu <- as.vector(predict(fit_mu, newdata = data.frame(argvals = argvals)))
  # mu <- smooth.spline(argvals, meanY, all.knots = FALSE)$y
  rm(meanY)
  
  
  ##################################################################################
  ## Estimate group-specific mean function (eta)
  ##################################################################################
  if(silent == FALSE) print("Step 2: Estimate group-specific mean function (eta)")
  
  mueta = matrix(0, L, J) 
  eta = matrix(0, L, J) ## matrix to store visit-specific means
  colnames(mueta) <- colnames(eta) <- levels(df$group)
  Ytilde <- matrix(NA, nrow = nrow(df$Y), ncol = ncol(df$Y))
  if(twoway==TRUE) {
    for(j in 1:J) {
      ind_j <- which(df$group == levels(df$group)[j])
      if(length(ind_j) > 1){
        meanYj <- colMeans(df$Y[ind_j,], na.rm=TRUE)
      }else{
        meanYj <- df$Y[ind_j,]
      }
      fit_mueta <- gam(meanYj ~ s(argvals))
      mueta[,j] <- predict(fit_mueta, newdata = data.frame(argvals = argvals))
      # mueta[,j] <- smooth.spline(argvals, meanYj, all.knots = FALSE)$y
      eta[,j] <- mueta[,j] - mu
      Ytilde[ind_j,] <- df$Y[ind_j,] - matrix(mueta[,j], nrow = length(ind_j), ncol = L, byrow = TRUE)
    }
    rm(meanYj, ind_j, j)
  } else{
    Ytilde <- df$Y - matrix(mu, nrow = nrow(df$Y), ncol = L, byrow = TRUE)
  }

  df$Ytilde <- I(Ytilde) ## Ytilde is the centered multilevel functional data
  rm(Ytilde)
  
  ##################################################################################
  ## specify the B-spline basis: knots
  ##################################################################################
  if(length(knots)==1){
    if(knots+p>=L) cat("Too many knots!\n")
    stopifnot(knots+p<L)
    
    K.p <- knots
    knots <- seq(-p, K.p+p, length=K.p+1+2*p)/K.p
    knots <- knots*(max(argvals)-min(argvals)) + min(argvals)
  }
  if(length(knots)>1) K.p <- length(knots)-2*p-1
  if(K.p>=L) cat("Too many knots!\n")
  stopifnot(K.p < L)
  c.p <- K.p + p #the number of B-spline basis functions
  
  ## Precalculation for smoothing
  List <- pspline.setting(argvals, knots, p, m)
  B <- List$B #B is the J ×c design matrix
  Bt <- t(as.matrix(B))
  Sigi.sqrt <- List$Sigi.sqrt #(t(B)B)^(-1/2)
  s <- List$s #eigenvalues of Sigi_sqrt%*%(P%*%Sigi_sqrt)
  U <- List$U #eigenvectors of Sigi_sqrt%*%(P%*%Sigi_sqrt)
  A0 <- Sigi.sqrt %*% U
  
  
  ##################################################################################
  ## impute missing data of Y using FACE approach
  ##################################################################################
  Ji <- as.numeric(table(df$id))
  diagD <- rep(Ji, Ji)
  smooth.Gt = face.Cov(Y=unclass(df$Ytilde), argvals, A0, Bt, s, c.p)
  ## impute missing data of Y using FACE approach
  if(sum(is.na(df$Ytilde))>0){
    df$Ytilde[which(is.na(df$Ytilde))] <- smooth.Gt$Yhat[which(is.na(df$Ytilde))]
  }
  if(weight=="subj"){
    YH <- unclass(df$Ytilde)*sqrt(nrow(df$Ytilde)/(I*diagD))
    smooth.Gt <- face.Cov(Y=YH, argvals, A0, Bt, s, c.p)
    rm(YH)
  }
  diag_Gt <- colMeans(df$Ytilde^2)
  
  
  ##################################################################################
  ## Estimate principal components of Kw (the within covariance)
  ##################################################################################
  if(silent == FALSE) print("Step 3: Estimate the within (Kw) subject covariance matrix")
  
  
  if(smooth=="joint"){
    ## method 1: smooth Gw jointly
    inx_row_ls <- split(1:nrow(df$Ytilde), f=factor(df$id, levels=unique(df$id)))
    weight <- sqrt(nrow(df$Ytilde)/(sum(diagD) - nrow(df$Ytilde)))
    Ysubm <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x,,drop=FALSE],na.rm=TRUE), numeric(L)))
    YR <-  do.call("rbind",lapply(1:I, function(x) {
      weight * t(t(sqrt(Ji[x])*df$Ytilde[inx_row_ls[[x]],,drop=FALSE]) - Ysubm[x,]/sqrt(Ji[x]))
    }))
    smooth.Gw <- face.Cov(Y=YR, argvals, A0, Bt, s, c.p)
    rm(Ji, diagD, inx_row_ls, weight, Ysubm, YR, Bt, s, Sigi.sqrt, U)
  } else{
    ## method 2: smooth two parts of Gw separately
    ### first part of formula (Shou et al.2015): t(Y) %*% D %*% Y
    #weight1 <- sqrt(nrow(df$Ytilde)/(sum(diagD) - nrow(df$Ytilde)))
    YD <- unclass(df$Ytilde)*sqrt(diagD)
    smooth.Gw1 <- face.Cov(Y = YD, argvals, A0, Bt, s, c.p)
    ### second part of formula: t(E %*% Y) %*% E %*% Y
    inx_row_ls <- split(1:nrow(df$Ytilde), f=factor(df$id, levels=unique(df$id)))
    #weight2 <- sqrt(I/(sum(diagD) - nrow(df$Ytilde)))
    YE <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x,,drop=FALSE],na.rm=TRUE), numeric(L)))
    smooth.Gw2 <- face.Cov(Y = YE, argvals, A0, Bt, s, c.p)
    # Obtain Kw
    temp = (nrow(YD)*smooth.Gw1$decom - nrow(YE)*smooth.Gw2$decom)/(sum(diagD) - nrow(df$Ytilde))
    Eigen <- eigen(temp,symmetric=TRUE)
    A <- Eigen$vectors
    Sigma <- Eigen$values
    d <- Sigma[1:c.p]
    d <- d[d>0]
    per <- cumsum(d)/sum(d)
    N <- ifelse (is.null(npc), min(which(per>pve)), min(npc, length(d)))
    A.N <- A[,1:N]
    smooth.Gw <- list(decom=temp, evalues=Sigma[1:N], evectors=as.matrix(B%*%A0%*%A.N))
    rm(Ji, diagD, inx_row_ls, YD, YE, smooth.Gw1, smooth.Gw2, Bt, s, temp, Eigen, A, Sigma, d, per, N, A.N)
  }

  
  
  ##################################################################################
  ## Estimate principal components of Kb (the between covariance)
  ##################################################################################
  if(silent == FALSE) print("Step 4: Estimate the between (Kb) subject covariance matrix")
  
  temp = smooth.Gt$decom - smooth.Gw$decom
  Eigen <- eigen(temp,symmetric=TRUE)
  A <- Eigen$vectors
  Sigma <- Eigen$values
  d <- Sigma[1:c.p]
  d <- d[d>0]
  per <- cumsum(d)/sum(d)
  N <- ifelse (is.null(npc), min(which(per>pve)), min(npc, length(d)))
  A.N <- A[,1:N]
  smooth.Gb <- list(evalues=Sigma[1:N], evectors=as.matrix(B%*%A0%*%A.N))
  rm(smooth.Gt, temp, Eigen, A, Sigma, d, per, N, B, A0, A.N)
  
  # Ji <- as.numeric(table(df$id))
  # diagD <- rep(Ji, Ji)
  # #weight <- sqrt(I/(sum(diagD) - nrow(df$Ytilde)))
  # weight <- sqrt(I/sum(diagD))
  # inx_row_ls <- split(1:nrow(df$Ytilde), f=factor(df$id, levels=unique(df$id)))
  # Ysubm <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x,,drop=FALSE],na.rm=TRUE), numeric(L))) 
  # #YH1 <- sqrt((Ji-1)/Ji) * weight * Ysubm
  # YH1 <- weight * Ysubm
  # smooth.Gb <- face.Cov(Y=YH1, argvals, A0, Bt, s, c.p)
  

  ###########################################################################################
  ## Estimate eigen values and eigen functions at two levels by calling the 
  ## eigen function (in R "base" package) on discretized covariance matrices.
  ###########################################################################################
  if(silent == FALSE) print("Step 6: Estimate eigen values and eigen functions at two levels")
  
  w <- quadWeights(argvals, method = "trapezoidal")
  Wsqrt <- sqrt(w)
  #efunctions <- list(level1=smooth.Gb$evectors/Wsqrt, level2=smooth.Gw$evectors/Wsqrt)
  efunctions <- list(level1=sqrt(L)*smooth.Gb$evectors, level2=sqrt(L)*smooth.Gw$evectors)
  evalues <- list(level1=smooth.Gb$evalues/L,level2=smooth.Gw$evalues/L)
  npc <- list(level1=length(evalues[[1]]), level2=length(evalues[[2]]))
  names(efunctions) <- names(evalues) <- names(npc) <- c("level1", "level2")
  rm(w, Wsqrt, smooth.Gb, smooth.Gw)
  
  

  ###################################################################
  # Estimate the measurement error variance (sigma^2)
  ###################################################################
  if(silent == FALSE) print("Step 7: Estimate the measurement error variance (sigma^2)")
  
  cov.hat <- lapply(c("level1", "level2"), function(x) colSums(t(efunctions[[x]]^2)*evalues[[x]]))
  T.len <- argvals[L] - argvals[1]
  T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))
  T1.max <- max(which(argvals <= argvals[L] - 0.25 * T.len))
  DIAG <- (diag_Gt - cov.hat[[1]] - cov.hat[[2]])[T1.min:T1.max]
  w2 <- quadWeights(argvals[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0) ## estimated measurement error variance
  rm(cov.hat, T.len, T1.min, T1.max, DIAG, w2)
  
  
  ###################################################################
  # Estimate the principal component scores
  ###################################################################
  if(silent == FALSE) print("Step 8: Estimate principal component scores -- may take some time")
  
  ## estimated subject-visit and subject-specific underlying smooth curves
  Xhat <- Xhat.subject <- matrix(0, nrow(df$Y), L) 
  phi1 <- efunctions[[1]] ## extract eigenfunctions for simpler notations
  phi2 <- efunctions[[2]]
  score1 <- matrix(0, I, npc[[1]]) ## matrices storing scores of two levels
  score2 <- matrix(0, nrow(df$Y), npc[[2]])
  
  unGroups <- unique(nGroups$numGroups) ## unique number of groups
  if(length(unGroups) < I){
    # ptm <- proc.time()
    for(j in 1:length(unGroups)){
      Jm <- unGroups[j]
      ## calculate block matrices
      if(sigma2 < 1e-4){
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      }else{
        if(length(evalues[[1]]==1)){
          A <- Jm * (t(phi1) %*% phi1) / sigma2 + 1 / evalues[[1]]
        }else{
          A <- Jm * (t(phi1) %*% phi1) / sigma2 + diag(1 / evalues[[1]])
        }
        B = matrix(rep(t(phi1) %*% phi2 / sigma2, Jm), nrow = npc[[1]])
        if(length(evalues[[2]]==1)){
          temp = ginv(t(phi2) %*% phi2 / sigma2 + 1 / evalues[[2]])
        }else{
          temp = ginv(t(phi2) %*% phi2 / sigma2 + diag(1 / evalues[[2]]))
        }
      }
      C <- t(B)
      invD = diag.block(temp, Jm)
      
      ## calculate inverse of each block components separately
      MatE <- ginv(A-B%*%invD%*%C)
      MatF <- -invD%*%C%*%MatE
      MatG <- -MatE%*%B%*%invD
      MatH <- invD - invD%*%C%*%MatG
      Mat1 <- cbind(MatE,MatG)
      Mat2 <- cbind(MatF,MatH)
      
      ## estimate the principal component scores using MME
      ind.Jm <- nGroups$id[which(nGroups$numGroups == Jm)]
      YJm <- matrix(df$Ytilde[which(df$id %in% ind.Jm),], ncol = L)
      int1 <- rowsum(df$Ytilde[which(df$id %in% ind.Jm),] %*% phi1, rep(1:length(ind.Jm), each = Jm))
      int2 <- t(matrix(t(df$Ytilde[which(df$id %in% ind.Jm),] %*% phi2), nrow = npc[[2]]*Jm))
      int <- cbind(int1, int2)
      if(sigma2 >= 1e-4){
        int <- int / sigma2
      }
      score1[which(nGroups$id %in% ind.Jm),] <- int %*% t(Mat1)
      score2[which(df$id %in% ind.Jm),] <- t(matrix(Mat2 %*% t(int), nrow = npc[[2]]))
      
      temp <- score1[which(nGroups$id %in% ind.Jm),] %*% t(phi1)
      Xhat.subject[which(df$id %in% ind.Jm),] <- temp[rep(1:length(ind.Jm), each = Jm),]
      Xhat[which(df$id %in% ind.Jm),] <- Xhat.subject[which(df$id %in% ind.Jm),] + 
        score2[which(df$id %in% ind.Jm),] %*% t(phi2)
    }
    for(g in 1:length(levels(df$group))){
      ind.group <- which(df$group == levels(df$group)[g])
      Xhat.subject[ind.group,] <- t(t(Xhat.subject[ind.group,]) + mu + eta[,levels(df$group)[g]])
      Xhat[ind.group,] <- t(t(Xhat[ind.group,]) + mu + eta[,levels(df$group)[g]])
    }
    # proc.time() - ptm ## 91 s
  }else{
    # ptm <- proc.time()
    for(m in 1:I){
      Jm <- nGroups[m, 2]  ## number of visits for mth subject
      ## calculate block matrices
      if(sigma2 < 1e-4){
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      }else{
        if(length(evalues[[1]]==1)){
          A <- Jm * (t(phi1) %*% phi1) / sigma2 + 1 / evalues[[1]]
        }else{
          A <- Jm * (t(phi1) %*% phi1) / sigma2 + diag(1 / evalues[[1]])
        }
        B = matrix(rep(t(phi1) %*% phi2 / sigma2, Jm), nrow = npc[[1]])
        if(length(evalues[[2]]==1)){
          temp = ginv(t(phi2) %*% phi2 / sigma2 + 1 / evalues[[2]])
        }else{
          temp = ginv(t(phi2) %*% phi2 / sigma2 + diag(1 / evalues[[2]]))
        }
      }
      C <- t(B)
      invD = diag.block(temp, Jm)
      
      ## calculate inverse of each block components separately
      MatE <- ginv(A-B%*%invD%*%C)
      MatF <- -invD%*%C%*%MatE
      MatG <- -MatE%*%B%*%invD
      MatH <- invD - invD%*%C%*%MatG
      Mat1 <- cbind(MatE,MatG)
      Mat2 <- cbind(MatF,MatH)
      
      ## estimate the principal component scores
      int1 <- colSums(matrix(df$Ytilde[df$id==ID[m],], ncol = L) %*% phi1)
      int2 <- matrix(df$Ytilde[df$id==ID[m],], ncol = L) %*% phi2
      if(sigma2 < 1e-4){
        int <- c(int1, as.vector(t(int2)))
      }else{
        int <- c(int1, as.vector(t(int2))) / sigma2
      }
      score1[m,] <- Mat1 %*% int
      score2[which(df$id==ID[m]),] <- matrix(Mat2 %*% int, ncol = npc[[2]], byrow=TRUE)
      for (j in which(df$id==ID[m])) {
        Xhat.subject[j,] <- as.matrix(mu) + eta[,df$group[j]] + as.vector(phi1 %*% score1[m,])
        Xhat[j,] <- Xhat.subject[j,] + as.vector(phi2 %*% score2[j,])
      }
    }
    # proc.time() - ptm ## 1752 s
  }
  scores <- list(level1 = score1, level2 = score2)
  
  rm(A, B, C, int, int1, int2, invD, Mat1, Mat2, MatE, MatF, MatG, MatH, temp, YJm,
     g, ind.group, ind.Jm, j, Jm, unGroups, phi1, phi2, score1, score2)
  
  
  ###################################################################
  # Organize the results
  ###################################################################
  if(silent == FALSE) print("Step 9: Organize the results")
  
  res <- list(Xhat = Xhat, Xhat.subject = Xhat.subject, mu = mu, eta = eta, scores = scores, 
              efunctions = efunctions, evalues = evalues, npc = npc, sigma2 = sigma2, Y = df$Y)
  
  rm(df, efunctions, eta, evalues, mueta, nGroups, npc, scores, Xhat, Xhat.subject, 
     argvals, diag_Gt, I, ID, J, mu, pve, L, sigma2)
  
  return(res)
}

