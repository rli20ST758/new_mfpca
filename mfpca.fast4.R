#' This function implements fast MFPCA on multilevel functional data
#' 
#' @param Y a multilevel functional dataset on a regular/irregular grid stored in a matrix
#' @param id a vector containing the id information
#' @param group a vector containing information used to identify groups/visits
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
#' @import MASS
#' @import simex
#' @import Matrix
#' @import rARPACK

mfpca.fast4 <- function(Y, id, group = NULL, argvals = NULL, pve = 0.99, npc = NULL, 
                        p = 3, m = 2, knots = 35, silent = TRUE){
  ## required packages
  library(refund)
  library(splines)
  library(MASS)
  library(simex)
  library(Matrix)
  library(rARPACK)
  # source("./code/backup.R")
  
  ##################################################################################
  ## Organize the input
  ##################################################################################
  if(silent == FALSE) print("Organize the input")
  
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
  ## Estimate population mean function (mu) and group-specific mean function (eta)
  ##################################################################################
  if(silent == FALSE) print("Estimate population and group-specific mean functions")
  
  meanY <- colMeans(df$Y, na.rm = TRUE)
  fit_mu <- gam(meanY ~ s(argvals))
  mu <- as.vector(predict(fit_mu, newdata = data.frame(argvals = argvals)))
  rm(meanY)
  
  mueta <- eta <- matrix(NA, L, J) ## eta is a matrix to store visit-specific means
  colnames(mueta) <- colnames(eta) <- levels(df$group)
  Ytilde <- matrix(NA, nrow = nrow(df$Y), ncol = ncol(df$Y))
  for(j in 1:J){
    ind_j <- which(df$group == levels(df$group)[j])
    if(length(ind_j) > 1){
      meanYj <- colMeans(df$Y[ind_j,], na.rm=TRUE)
    }else{
      meanYj <- df$Y[ind_j,]
    }
    fit_mueta <- gam(meanYj ~ s(argvals))
    mueta[,j] <- predict(fit_mueta, newdata = data.frame(argvals = argvals))
    eta[,j] <- mueta[,j] - mu
    Ytilde[ind_j,] <- df$Y[ind_j,] - matrix(mueta[,j], nrow = length(ind_j), ncol = L, byrow = TRUE)
  }
  df$Ytilde <- I(Ytilde) ## Ytilde is the centered multilevel functional data
  rm(Ytilde, meanYj, ind_j, j)
  
  
  ##################################################################################
  ## FACE preparation: see Xiao et al. (2016) for details
  ##################################################################################
  
  ## Specify the knots of B-spline basis
  if(length(knots) == 1){
    if(knots + p >= L) cat("Too many knots!\n")
    stopifnot(knots + p < L)
    K.p <- knots
    knots <- seq(-p, K.p+p, length = K.p+1+2*p) / K.p
    knots <- knots * (max(argvals) - min(argvals)) + min(argvals)
  }
  if(length(knots) > 1) K.p <- length(knots) - 2*p - 1
  if(K.p >= L) cat("Too many knots!\n")
  stopifnot(K.p < L)
  c.p <- K.p + p ## the number of B-spline basis functions
  
  ## Precalculation for smoothing
  List <- pspline.setting(argvals, knots, p, m)
  B <- List$B ## the L × c design matrix
  Bt <- t(as.matrix(B))
  Sigi.sqrt <- List$Sigi.sqrt ## (t(B)B)^(-1/2)
  s <- List$s ## eigenvalues of Sigi_sqrt %*% (P %*% Sigi_sqrt)
  U <- List$U ## eigenvectors of Sigi_sqrt %*% (P %*% Sigi_sqrt)
  A0 <- Sigi.sqrt %*% U
  
  
  ##################################################################################
  ## Estimate the total covariance matrix (Kt)
  ##################################################################################
  if(silent == FALSE) print("Estimate the total covariance matrix")
  
  ## impute missing data of Y using FACE approach
  if(sum(is.na(df$Ytilde))>0){
    impute <- face.Cov(Y = unclass(df$Ytilde), argvals, A0, Bt, s, c.p)
    df$Ytilde[which(is.na(df$Ytilde))] <- impute$Yhat[which(is.na(df$Ytilde))]
    rm(impute)
  }
  
  Ji <- nGroups$numGroups
  diagD <- rep(Ji, Ji)
  YH <- unclass(df$Ytilde)*sqrt(nrow(df$Ytilde)/(I*diagD))
  smooth.Gt <- face.Cov(Y = YH, argvals, A0, Bt, s, c.p, Cov=T)
  Kt <- (t(smooth.Gt$Ktilde) + smooth.Gt$Ktilde) / 2
  diag_Gt <- colMeans(df$Ytilde^2)
  rm(YH, smooth.Gt)
  
  
  ##################################################################################
  ## Estimate the level-specific covariance matrix (Kb and Kw)
  ##################################################################################
  if(silent == FALSE) print("Estimate the level-specific covariance matrices")
  
  ## method 1: smooth Gw jointly
  # inx_row_ls <- split(1:nrow(df$Ytilde), f = factor(df$id, levels = unique(df$id)))
  # weight <- sqrt(nrow(df$Ytilde) / (sum(diagD) - nrow(df$Ytilde)))
  # Ysubm <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x,,drop=FALSE],na.rm=TRUE), numeric(L))) # The visit mean per subject 
  # YR <-  do.call("rbind", lapply(1:I, function(x) {
  #   weight * t(t(sqrt(Ji[x])*df$Ytilde[inx_row_ls[[x]],,drop=FALSE]) - Ysubm[x,]/sqrt(Ji[x]))
  #   }))
  # smooth.Gw <- face.Cov(Y = YR, argvals, A0, Bt, s, c.p, Cov=T)
  # Kw <- (t(smooth.Gw$Ktilde) + smooth.Gw$Ktilde) / 2
  # Kb <- Kt - Kw  ## the smoothed between-subject covariance matrix
  # rm(Ji, diagD, inx_row_ls, weight, Ysubm, YR, smooth.Gw, Kt, B, Bt, s, Sigi.sqrt, U, A0)
  
  
  # ## method 2: smooth two parts of Gw separately
  # ### first part of formula (Shou et al.2015): t(Y) %*% D %*% Y
  YD <- unclass(df$Ytilde)*sqrt(diagD)
  smooth.Gw1 <- face.Cov(Y = YD, argvals, A0, Bt, s, c.p, Cov=T)
  Kw1 <- nrow(YD) * smooth.Gw1$Ktilde
  ### second part of formula: t(E %*% Y) %*% E %*% Y
  inx_row_ls <- split(1:nrow(df$Ytilde), f=factor(df$id, levels=unique(df$id)))
  YE <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x,,drop=FALSE],na.rm=TRUE), numeric(L)))
  smooth.Gw2 <- face.Cov(Y = YE, argvals, A0, Bt, s, c.p, Cov=T)
  Kw2 <- nrow(YE) * smooth.Gw2$Ktilde
  Kw <- (Kw1 - Kw2) / (sum(diagD) - nrow(df$Ytilde))
  Kw <- (t(Kw) + Kw)/2
  Kb <- Kt - Kw  ## the smoothed between-subject covariance matrix
  rm(Ji, diagD, inx_row_ls, YD, YE, Kw1, Kw2, smooth.Gw1, smooth.Gw2, Kt, B, Bt, s, Sigi.sqrt, U, A0)
  
  
  ###########################################################################################
  ## Estimate the level-specific eigenvalues and eigenfunctions
  ###########################################################################################
  if(silent == FALSE) print("Estimate level-specific eigenvalues and eigenfunctions")
  
  w <- quadWeights(argvals, method = "trapezoidal")
  Wsqrt <- sqrt(w)
  npc.0wb <- list(level1 = Kb, level2 = Kw)  
  W0 <- (Wsqrt) %*% t(Wsqrt)
  V <- lapply(npc.0wb, function(x) W0*x)
  
  ecomp <- lapply(V, function(x) eigen(x, only.values = TRUE, symmetric = TRUE)) ## get npc
  evalues <- lapply(ecomp, function(x) replace(x$values, which(x$values <= 0), 0))
  npc <- lapply(evalues, function(x) ifelse(is.null(npc), min(which(cumsum(x)/sum(x) > pve)), npc))
  ecomp <- lapply(names(V), function(x) eigs_sym(V[[x]], npc[[x]], which="LM")) ## get the first npc eigenvectors
  efunctions <- lapply(1:2, function(x) 
    matrix((1/Wsqrt)*(ecomp[[x]])$vectors, nrow=L, ncol=npc[[x]]))
  evalues <- lapply(names(V), function(x) (evalues[[x]])[1:npc[[x]]])
  names(efunctions) <- names(evalues) <- names(npc) <- c("level1", "level2")
  rm(Kb, Kw, w, Wsqrt, npc.0wb, V, ecomp, W0)
  
  
  ###################################################################
  # Estimate the measurement error variance (sigma^2)
  ###################################################################
  if(silent == FALSE) print("Estimate the measurement error variance")
  
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
  if(silent == FALSE) print("Estimate principal component scores")
  
  ## estimated subject-visit and subject-specific underlying smooth curves
  Xhat <- Xhat.subject <- matrix(0, nrow(df$Y), L) 
  phi1 <- efunctions[[1]] ## extract eigenfunctions for simpler notations
  phi2 <- efunctions[[2]]
  score1 <- matrix(0, I, npc[[1]]) ## matrices storing scores of two levels
  score2 <- matrix(0, nrow(df$Y), npc[[2]])
  
  unGroups <- unique(nGroups$numGroups) ## unique number of groups
  if(length(unGroups) < I){
    for(j in 1:length(unGroups)){
      Jm <- unGroups[j]
      ## calculate block matrices
      if(sigma2 < 1e-4){
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      }else{
        A <- Jm * (t(phi1) %*% phi1) / sigma2 + diag(1 / evalues[[1]])
        B = matrix(rep(t(phi1) %*% phi2 / sigma2, Jm), nrow = npc[[1]])
        temp = ginv(t(phi2) %*% phi2 / sigma2 + diag(1 / evalues[[2]]))
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
  }else{
    for(m in 1:I){
      Jm <- nGroups[m, 2]  ## number of visits for mth subject
      ## calculate block matrices
      if(sigma2 < 1e-4){
        A <- Jm * (t(phi1) %*% phi1)
        B <- matrix(rep(t(phi1) %*% phi2, Jm), nrow = npc[[1]])
        temp <- ginv(t(phi2) %*% phi2)
      }else{
        A <- Jm * (t(phi1) %*% phi1) / sigma2 + diag(1 / evalues[[1]])
        B = matrix(rep(t(phi1) %*% phi2 / sigma2, Jm), nrow = npc[[1]])
        temp = ginv(t(phi2) %*% phi2 / sigma2 + diag(1 / evalues[[2]]))
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
  }
  scores <- list(level1 = score1, level2 = score2)
  
  rm(A, B, C, int, int1, int2, invD, Mat1, Mat2, MatE, MatF, MatG, MatH, temp, YJm,
     g, ind.group, ind.Jm, j, Jm, unGroups, phi1, phi2, score1, score2)
  
  
  ###################################################################
  # Organize the results
  ###################################################################
  if(silent == FALSE) print("Organize the results")
  
  res <- list(Xhat = Xhat, Xhat.subject = Xhat.subject, mu = mu, eta = eta, scores = scores, 
              efunctions = efunctions, evalues = evalues, npc = npc, sigma2 = sigma2, Y = df$Y)
  
  rm(df, efunctions, eta, evalues, mueta, nGroups, npc, scores, Xhat, Xhat.subject, 
     argvals, diag_Gt, I, ID, J, mu, pve, L, sigma2)
  
  return(res)
}

