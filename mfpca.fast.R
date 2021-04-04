#' This function performs multilevel FPCA (Di et al.2009) on multilevel functional data
#' It is based on the "mfpca.sc2()" function from Ruonan Li and Luo Xiao, 
#' and the "mfpca.sc()" function from Julia Wrobel in the R refund package.
#' 
#' @param Y a multilevel/longitudinal functional dataset on a regular grid stored in a matrix
#' @param id a vector containing the id information
#' @param group a vector containing information used to identify groups/visits
#' @param argvals a vector containing observed locations on the functional domain
#' @param pve proportion of variance explained: used to choose the number of principal components.
#' @param npc prespecified value for the number of principal components (if given, this overrides \code{pve}).
#' 
#' @export
#' @import refund
#' @import splines
#' @import MASS
#' @import simex

mfpca.fast <- function(Y, id, group = NULL, argvals = NULL, pve = 0.99, npc = NULL, silent = TRUE){
  ## required packages
  library(refund)
  library(splines)
  library(MASS)
  library(simex)
  source("./code/fbps.cov.R")
  
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
  
  ## impute missing data of Y using FACE approach
  if(!is.null(is.na(Y))){
    Y[which(is.na(Y))] <- fpca.face(Y)$Yhat[which(is.na(Y))]
  }
  
  ## organize data into one data frame
  df <- data.frame(id = id, group = group, Y = I(Y))
  rm(id, group, Y)
  
  ## derive several variables that will be used later
  J <- length(levels(df$group)) ## number of groups
  S <- ncol(df$Y) ## number of observations along the domain
  nknots <- min(100, round(0.1*S)) ## number of knots used for general smoothing
  nGroups <- data.frame(table(df$id))  ## calculate number of groups for each subject
  colnames(nGroups) = c("id", "numGroups")
  ID = sort(unique(df$id)) ## id of each subject
  I <- length(ID) ## number of subjects
  ## assume observations are equally-spaced on [0,1] if not specified
  if (is.null(argvals))  argvals <- seq(0, 1, length.out=S) 
  
  
  ##################################################################################
  ## Estimate population mean function (mu)
  ##################################################################################
  if(silent == FALSE) print("Step 1: Estimate population mean function (mu)")
  
  meanY <- colMeans(df$Y, na.rm = TRUE)
  mu <- smooth.spline(argvals, meanY, all.knots = FALSE)$y
  rm(meanY)
  
  
  ##################################################################################
  ## Estimate group-specific mean function (eta)
  ##################################################################################
  if(silent == FALSE) print("Step 2: Estimate group-specific mean function (eta)")
  
  mueta = matrix(NA, S, J) 
  eta = matrix(NA, S, J) ## matrix to store visit-specific means
  colnames(mueta) <- colnames(eta) <- levels(df$group)
  Ytilde <- matrix(NA, nrow = nrow(df$Y), ncol = ncol(df$Y))
  for(j in 1:J){
    ind_j <- which(df$group == levels(df$group)[j])
    if(length(ind_j) > 1){
      meanYj <- colMeans(df$Y[ind_j,], na.rm=TRUE)
    }else{
      meanYj <- df$Y[ind_j,]
    }
    mueta[,j] <- smooth.spline(argvals, meanYj, all.knots = FALSE)$y
    eta[,j] <- mueta[,j] - mu
    Ytilde[ind_j,] <- df$Y[ind_j,] - matrix(mueta[,j], nrow = length(ind_j), ncol = S, byrow = TRUE)
  }
  df$Ytilde <- I(Ytilde) ## Ytilde is the centered multilevel functional data
  rm(Ytilde, meanYj, ind_j, j)
  
  
  ##################################################################################
  ## Estimate Kt, the total covariance
  ##################################################################################
  if(silent == FALSE) print("Step 3: Estimate the total covariance matrix (Kt)")
  
  ptm <- proc.time()
  cov.sum <- crossprod(unclass(df$Ytilde))
  Gt <- cov.sum/nrow(df$Ytilde) ## estimate Gt using MoM estimator
  
  diag_Gt <- diag(Gt) 
  diag(Gt) <- NA
  npc.0 <- fbps.cov(Gt, diag.remove=T, knots=nknots)$cov
  Kt <- (npc.0 + t(npc.0))/2 ## the smoothed (total) covariance matrix
  proc.time() - ptm ## 106 s
  
  rm(cov.sum, npc.0, Gt)
  
  ##################################################################################
  ## Estimate Kb (the between covariance) and Kw (the within covariance)
  ##################################################################################
  if(silent == FALSE) print("Step 4: Estimate the between (Kb) and within (Kw) subject covariance matrix")
  
  ## first part of formula (Shou et al.2015): t(Y) %*% D %*% Y
  Ji <- table(df$id)
  diagD <- rep(Ji, Ji)
  m1 <- crossprod(unclass(df$Ytilde)*sqrt(diagD))
  
  ## second part of formula: t(E %*% Y) %*% E %*% Y
  inx_row_ls <- split(1:nrow(df$Ytilde), f=factor(df$id, levels=unique(df$id)))
  Ytilde_subj <- t(vapply(inx_row_ls, function(x) colSums(df$Ytilde[x,,drop=FALSE],na.rm=TRUE), numeric(S)))
  m2 <- crossprod(Ytilde_subj)
  
  Gw <- (m1 - m2) / (sum(diagD) - nrow(df$Ytilde)) ## estimate Gw using MoM estimator
  
  npc.0w <- fbps.cov(Gw, diag.remove=T, knots=nknots)$cov
  Kw <- (npc.0w + t(npc.0w))/2 ## the smoothed within covariance matrix
  Kb <- (Kt - Kw + t(Kt - Kw))/2 ## the smoothed between covariance matrix
  
  rm(Ji, diagD, m1, m2, inx_row_ls, Ytilde_subj, Gw, npc.0w)
  
  # ##################################################################################
  # ## Estimate Kb, the between covariance
  # ### Use the possible pairs of same-subject visits to calculate the between covariance function
  # ##################################################################################
  # if(silent == FALSE) print("Step 4: Estimate the between subject covariance matrix (Kb) -- may take several minutes")
  # 
  # # ptm <- proc.time()
  # cov.sum <- matrix(0, S, S)
  # ids.Kb <- nGroups[nGroups$numGroups > 1, c("id")] ## subject ids with at least 2 visits
  # for (m in 1:I) {
  #   if (ID[m] %in% ids.Kb) {
  #     mat_m <- matrix(df$Ytilde[df$id==ID[m],], ncol=S)
  #     cov.sum <- cov.sum + tcrossprod(colSums(mat_m)) - crossprod(mat_m)
  #   }
  # }
  # # proc.time() - ptm ## 682 s
  # 
  # Gb <- cov.sum/sum(nGroups[,2]*(nGroups[,2]-1)) ## between covariance
  # npc.0b <- fbps.cov(Gb, diag.remove=T, knots=nknots)$cov
  # Kb <- (npc.0b + t(npc.0b))/2 ## smoothed (between) covariance matrix
  # 
  # rm(cov.sum, npc.0b, ids.Kb, mat_m, m, Gb)
  # 
  # 
  # ###########################################################################################
  # ## Estimate Kw, the within covariance
  # ############################################################################################
  # if(silent == FALSE) print("Step 5: Estimate the within subject covariance matrix (Kw)")
  # 
  # Kw <- (Kt - Kb + t(Kt - Kb))/2 ## to ensure symmetric
  # 
  
  ###########################################################################################
  ## Estimate eigen values and eigen functions at two levels by calling the 
  ## eigen function (in R "base" package) on discretized covariance matrices.
  ###########################################################################################
  if(silent == FALSE) print("Step 6: Estimate eigen values and eigen functions at two levels")
  
  w <- quadWeights(argvals, method = "trapezoidal")
  Wsqrt <- diag(sqrt(w))
  Winvsqrt <- diag(1/(sqrt(w)))
  
  npc.0wb <- list(level1 = Kb, level2 = Kw)   
  V <- lapply(npc.0wb, function(x) Wsqrt %*% x %*% Wsqrt)
  evalues <- lapply(V, function(x) eigen(x, symmetric = TRUE, only.values = TRUE)$values)
  evalues <- lapply(evalues, function(x) replace(x, which(x <= 0), 0))
  npc <- lapply(evalues, function(x) ifelse(is.null(npc), min(which(cumsum(x)/sum(x) > pve)), npc ))
  efunctions <- lapply(names(V), function(x) 
    matrix(Winvsqrt%*%eigen(V[[x]], symmetric=TRUE)$vectors[,seq(len=npc[[x]])], nrow=S, ncol=npc[[x]]))
  evalues <- lapply(names(V), function(x) eigen(V[[x]], symmetric=TRUE, only.values=TRUE)$values[1:npc[[x]]])
  names(efunctions) <- names(evalues) <- names(npc) <- c("level1", "level2")
  
  rm(w, Wsqrt, Winvsqrt, npc.0wb, V)
  
  
  ###################################################################
  # Estimate the measurement error variance (sigma^2)
  ###################################################################
  if(silent == FALSE) print("Step 7: Estimate the measurement error variance (sigma^2)")
  
  cov.hat <- lapply(c("level1", "level2"), 
                    function(x) efunctions[[x]] %*% diag(evalues[[x]], npc[[x]], npc[[x]]) %*% t(efunctions[[x]]))
  T.len <- argvals[S] - argvals[1]
  T1.min <- min(which(argvals >= argvals[1] + 0.25 * T.len))
  T1.max <- max(which(argvals <= argvals[S] - 0.25 * T.len))
  DIAG <- (diag_Gt - diag(cov.hat[[1]])-diag(cov.hat[[2]]))[T1.min:T1.max]
  w2 <- quadWeights(argvals[T1.min:T1.max], method = "trapezoidal")
  sigma2 <- max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0) ## estimated measurement error variance
  
  rm(cov.hat, T.len, T1.min, T1.max, DIAG, w2)
  
  
  ###################################################################
  # Estimate the principal component scores
  ###################################################################
  if(silent == FALSE) print("Step 8: Estimate principal component scores -- may take some time")
  
  ## estimated subject-visit and subject-specific underlying smooth curves
  Xhat <- Xhat.subject <- matrix(0, nrow(df$Y), S) 
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
      YJm <- matrix(df$Ytilde[which(df$id %in% ind.Jm),], ncol = S)
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
      int1 <- colSums(matrix(df$Ytilde[df$id==ID[m],], ncol = S) %*% phi1)
      int2 <- matrix(df$Ytilde[df$id==ID[m],], ncol = S) %*% phi2
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
  
  rm(df, efunctions, eta, evalues, Kb, Kt, Kw, mueta, nGroups, npc, scores, Xhat, Xhat.subject, 
     argvals, diag_Gt, I, ID, J, mu, nknots, pve, S, sigma2)
  
  return(res)
}




## a function supposed to be in the refund package
quadWeights<- function(argvals, method = "trapezoidal")
{
  ret <- switch(method,
                trapezoidal = {D <- length(argvals)
                1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])},
                midpoint = c(0,diff(argvals)),
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
  
  return(ret)  
}

