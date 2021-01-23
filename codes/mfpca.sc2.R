

################################################################################################
## @param design: default as "irregular". 
##                "regular" design refers to dense and regularly sampled functional data;
##                "irregular" design refers to sparse and irregularly sampled functional data.
## @param weight: default as "observations". It only works for regular data.
##                "obs": the sample covariance is weighted equally by observations ;
##                "subj": the sample covariance is weighted equally by subjects;.
## @param nknots: the number of knots with default value of 35. 
################################################################################################

mfpca.sc2 <- function(Y = NULL, id = NULL, visit = NULL, twoway = FALSE, design = "irregular",
                      weight = "obs", argvals = NULL, nbasis = 10, nknots = 35, pve = 0.99, 
                      npc = NULL, center = TRUE, integration = "trapezoidal") {
  
  library(simex)
  if (!is.null(visit)){
    visit = as.integer(factor(visit))
  }else{ visit = ave(id, id, FUN=seq_along)}
  ## for irregular data, don't consider weights by subjects
  if (design=="irregular") weight = "obs" 
  Y.df = data.frame(id=id, visit=visit)
  Y.df$Y = Y
  J = length(unique(visit))  ## gets max number of visits
  ID = sort(unique(id))
  M = length(ID) ## number of subjects
  D = NCOL(Y)  
  I = NROW(Y) 
  nVisits = data.frame(table(id))  ## calculate number of visitis for each subject
  colnames(nVisits) = c("id", "numVisits")
  if (is.null(argvals))  argvals = seq(0, 1, length.out=D)  
  if (center && weight=="obs") {
    meanY = colMeans(Y, na.rm=TRUE)
    mu = smooth.spline(argvals,meanY)$y
    #d.vec = rep(argvals, each = I)
    #gam0 = gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    #mu = predict(gam0, newdata = data.frame(d.vec = argvals))  
  } else if (center && weight=="subj") {
    meanY = Reduce("+", lapply(1:M, function(m){
      temp = subset(Y.df, id==ID[m])$Y
      if (nrow(temp>1)) colMeans(temp)
      else temp}))/M  
    mu = smooth.spline(argvals,meanY)$y
  } else {
    mu = rep(0, D)
  }
  
  mueta = matrix(0, D, J) ## matrix to store visit-specific means
  eta = matrix(0, D, J)
  if (center && twoway) {
    Y.tilde <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
    for(j in 1:J){
      ind.j <- which(Y.df$visit==j)
      if(length(ind.j)>1) {
        fit.eta = colMeans(Y.df$Y[ind.j,], na.rm=TRUE)
      } else {
        fit.eta = Y.df$Y
      }
      mueta[, j] = smooth.spline(argvals,fit.eta)$y
      eta[,j] = mueta[,j] - mu
      Y.tilde[ind.j,] <- Y.df$Y[ind.j,] - matrix(mueta[,j], nrow = length(ind.j), ncol = D, byrow = TRUE)
    }
  } else{
    Y.tilde = Y - matrix(mu, I, D, byrow = TRUE)  ## subtract the mean function from Y
  }
  
  ## split Y.tilde by subjects
  Y.tilde.split = lapply(1:M, function(m){matrix(Y.tilde[id==ID[m],],nc=D)}) 

  
  ##################################################################################
  #      "regular" design
  ###     Smoothing sample covariances by sandwich smoother.
  ##################################################################################
  if (design == "regular") {
    
    ##################################################################################
    #      CALCULATE Kt, THE TOTAL COVARIANCE
    ##################################################################################
    if (weight=="obs") {
      cov.sum = crossprod(Y.tilde)
      G.0 = cov.sum/I
    } else {
      cov.sum = lapply(1:M, function(m){crossprod(Y.tilde.split[[m]])/nVisits[m,2]})
      G.0 = Reduce("+",cov.sum)/M
    } ## end of if(weight=="obs")
    diag.G0 = diag(G.0) 
    diag(G.0) = NA
    npc.0 = fbps.cov(G.0, diag.remove=T, knots=nknots)$cov
    npc.0 = (npc.0 + t(npc.0))/2 ## the smoothed (total) covariance matrix
    
    ##################################################################################
    #      CALCULATE Kb, THE BETWEEN COVARIANCE
    ###     Calculate the pairs to calculate the between covariance function
    ###     first, calculate the possible pairs of same-subject visits
    ##################################################################################
    cov.sum = matrix(0, D, D)
    ids.KB = nVisits[nVisits$numVisits > 1, c("id")]
    if (weight=="obs") {
      for (m in 1:M) {
        if (ID[m] %in% ids.KB) { ## check if mth subject has at least 2 visits
          cov.sum = cov.sum  + tcrossprod(colSums(Y.tilde.split[[m]])) - crossprod(Y.tilde.split[[m]])
        }
      }
      G.0b = cov.sum/sum(nVisits[,2]*(nVisits[,2]-1)) ## between covariance
    } else {
      for (m in 1:M) {
        if (ID[m] %in% ids.KB) { ## check if mth subject has at least 2 visits
          cov.sum = cov.sum + (tcrossprod(colSums(Y.tilde.split[[m]]))-crossprod(Y.tilde.split[[m]]))/(nVisits[m,2]*(nVisits[m,2]-1))
        }
      }
      G.0b = cov.sum/M ## between covariance
    } ## end of if(weight=="obs")
    
    npc.0b = fbps.cov(G.0b, diag.remove=T, knots=nknots)$cov
    npc.0b = (npc.0b + t(npc.0b))/2 ##  smoothed (between) covariance matrix
    
  } ## end of if(design == "regular")
  
  
  ##################################################################################
  #      "irregular" design
  ###     Smoothing sample covariances by gam function.
  ##################################################################################
  if (design == "irregular") {
    
    ##################################################################################
    #      CALCULATE Kt, THE TOTAL COVARIANCE
    ##################################################################################
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    for (m in 1:M) {
      for (j in 1:nVisits[m,2]) {
        obs.points = which(!is.na(Y.tilde.split[[m]][j,]))
        cov.count[obs.points, obs.points] <- cov.count[obs.points, obs.points] + 1
        cov.sum[obs.points, obs.points] <- cov.sum[obs.points, obs.points] + tcrossprod(Y.tilde.split[[m]][j, obs.points])
      }
    } 
    
    G.0 = ifelse(cov.count == 0, NA, cov.sum/cov.count)  
    diag.G0 = diag(G.0) 
    diag(G.0) = NA
    row.vec = rep(argvals, each = D)
    col.vec = rep(argvals, D)
    s.npc.0 =predict(gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis), weights = as.vector(cov.count)), 
                     newdata = data.frame(row.vec = row.vec, col.vec = col.vec))
    npc.0 = matrix(s.npc.0, D, D)
    npc.0 = (npc.0 + t(npc.0))/2  ## the smoothed (total) covariance matrix
    
    
    ##################################################################################
    #      CALCULATE Kb, THE BETWEEN COVARIANCE
    ###     Calculate the pairs to calculate the between covariance function
    ###     first, calculate the possible pairs of same-subject visits
    ##################################################################################
    cov.sum = cov.count = cov.mean = matrix(0, D, D)
    ids.KB = nVisits[nVisits$numVisits > 1, c("id")]
    
    for (m in 1:M) {
      if (ID[m] %in% ids.KB) { ## check if mth subject has at least 2 visits
        for (j in 1:(nVisits[m, 2]-1)) {
          obs.points1 = which(!is.na(Y.tilde.split[[m]][j,]))
          for(k in (j+1):nVisits[m, 2]) {
            obs.points2 = which(!is.na(Y.tilde.split[[m]][k,]))
            cov.count[obs.points1, obs.points2] <- cov.count[obs.points1, obs.points2] + 1
            cov.sum[obs.points1, obs.points2] <- cov.sum[obs.points1, obs.points2] + tcrossprod(Y.tilde.split[[m]][j, obs.points1], Y.tilde.split[[m]][k, obs.points2]) 
            cov.count[obs.points2, obs.points1] <- cov.count[obs.points2, obs.points1] + 1
            cov.sum[obs.points2, obs.points1] <- cov.sum[obs.points2, obs.points1] + tcrossprod(Y.tilde.split[[m]][k, obs.points2], Y.tilde.split[[m]][j, obs.points1])                                                        
          }
        }
      }
    }
    
    G.0b <- ifelse(cov.count==0, NA,  cov.sum/cov.count)  ## between covariance
    row.vec = rep(argvals, each = D)  
    col.vec = rep(argvals, D)    
    s.npc.0b =predict(gam(as.vector(G.0b) ~ te(row.vec, col.vec, k = nbasis), weights = as.vector(cov.count)), 
                      newdata = data.frame(row.vec = row.vec, col.vec = col.vec))
    npc.0b = matrix(s.npc.0b, D, D)
    npc.0b = (npc.0b + t(npc.0b))/2  ##  smoothed (between) covariance matrix
    
  } ## end of if(design == "irregular")
  
  
  ###########################################################################################
  #      CALCULATE Kw, THE WITHIN COVARIANCE
  ### use smoothed covariance total and covariance between to calculated covariance within
  ############################################################################################
  npc.0w = npc.0 - npc.0b
  npc.0w = (npc.0w + t(npc.0w))/2  ## smoothed within covariance matrix
  
  
  ###########################################################################################
  ###     Estimate eigen values and eigen functions at two levels by calling the 
  ###     eigen function (in R "base" package) on discretized covariance matrices.
  ###########################################################################################
  w = quadWeights(argvals, method=integration)
  Wsqrt = diag(sqrt(w))
  Winvsqrt = diag(1/(sqrt(w)))
  
  npc.0wb = list(level1 = npc.0b, level2 = npc.0w)   
  V = lapply(npc.0wb, function(x) Wsqrt %*% x %*% Wsqrt)
  evalues = lapply(V, function(x) eigen(x, symmetric = TRUE, only.values = TRUE)$values)
  evalues = lapply(evalues, function(x) replace(x, which(x <= 0), 0))
  npc = lapply(evalues, function(x) ifelse(is.null(npc), min(which(cumsum(x)/sum(x) >  pve)), npc ))
  efunctions = lapply(names(V), function(x) 
    matrix(Winvsqrt%*%eigen(V[[x]], symmetric=TRUE)$vectors[,seq(len=npc[[x]])],nrow=D,ncol=npc[[x]]))
  evalues = lapply(names(V), function(x) eigen(V[[x]], symmetric=TRUE, only.values=TRUE)$values[1:npc[[x]]])
  names(efunctions) = names(evalues) = names(npc) = c("level1", "level2")
  
  
  ###################################################################
  # calculate the measurement error variance
  ###################################################################
  cov.hat = lapply(c("level1", "level2"), 
                   function(x) efunctions[[x]] %*% diag(evalues[[x]], npc[[x]], npc[[x]]) %*% t(efunctions[[x]]))
  T.len = argvals[D] - argvals[1]
  T1.min = min(which(argvals >= argvals[1] + 0.25 * T.len))
  T1.max = max(which(argvals <= argvals[D] - 0.25 * T.len))
  DIAG = (diag.G0 - diag(cov.hat[[1]])-diag(cov.hat[[2]]))[T1.min:T1.max]
  w2 = quadWeights(argvals[T1.min:T1.max], method = integration)
  sigma2 = max(weighted.mean(DIAG, w = w2, na.rm = TRUE), 0) ### estimated measurement error variance
  
 
  Yhat = Yhat.subject = matrix(0, I, D)
  Z1 = efunctions[[1]]
  Z2 = efunctions[[2]]
  score1 = matrix(0, M, npc[[1]])
  score2 = matrix(0, I, npc[[2]])
  
  
  ###################################################################
  ### Finally, calculate the principal component scores
  ###################################################################
  row.split = lapply(1:M, function(m){which(id==ID[m])}) 
  if (design == "regular") {
    # if sigma2 is very small
    if (sigma2<1e-4) {
      for(m in 1:M) {
        Jm = nVisits[m, 2]  ## number of visits for mth subject
        A = Jm*(t(Z1)%*%Z1)
        B = matrix(rep(t(Z1)%*%Z2,Jm), nr=npc[[1]])
        C = t(B)
        temp = ginv(t(Z2)%*%Z2)
        #invD = as.matrix(bdiag(lapply(1:Jm, function(j){temp})))
        invD = diag.block(temp, Jm)
        MatE = ginv(A-B%*%invD%*%C)
        MatF = -invD%*%C%*%MatE
        MatG = -MatE%*%B%*%invD
        MatH = invD - invD%*%C%*%MatG
        Mat1 = cbind(MatE,MatG)
        Mat2 = cbind(MatF,MatH)
        
      
        int1 = colSums(Y.tilde.split[[m]]%*%Z1)
        int2 = Y.tilde.split[[m]]%*%Z2
        int = c(int1,as.vector(t(int2)))
        score1[m, ] = Mat1 %*% int
        score2[row.split[[m]], ] = matrix(Mat2%*%int,nc=npc[[2]],byrow=T)
        for (j in row.split[[m]]) {
          Yhat.subject[j, ] = as.matrix(mu) + eta[ ,visit[j]] + as.vector(Z1%*%score1[m, ])
          Yhat[j, ] = Yhat.subject[j, ] + as.vector(Z2%*%score2[j, ])
        }
      } #end m
    } else {
      for(m in 1:M) {
        Jm = nVisits[m, 2]  ## number of visits for mth subject
        A =  Jm*(t(Z1)%*%Z1)/sigma2 + diag(1/evalues[[1]])
        B = matrix(rep(t(Z1)%*%Z2/sigma2,Jm),nr=npc[[1]])
        C = t(B)
        temp = ginv(t(Z2)%*%Z2/sigma2 + diag(1/evalues[[2]]))
        #invD = as.matrix(bdiag(lapply(1:Jm, function(j){temp})))
        invD = diag.block(temp, Jm)
        MatE = ginv(A-B%*%invD%*%C)
        MatF = -invD%*%C%*%MatE
        MatG = -MatE%*%B%*%invD
        MatH = invD - invD%*%C%*%MatG
        Mat1 = cbind(MatE,MatG)
        Mat2 = cbind(MatF,MatH)
        
  
        int1 = colSums(Y.tilde.split[[m]]%*%Z1)
        int2 = Y.tilde.split[[m]]%*%Z2
        int = c(int1,as.vector(t(int2)))/sigma2
        score1[m,] = Mat1 %*% int
        score2[row.split[[m]], ] = matrix(Mat2%*%int,nc=npc[[2]],byrow=T)
        for (j in row.split[[m]]) {
          Yhat.subject[j, ] = mu + eta[,visit[j]] + as.vector(Z1%*%score1[m, ])
          Yhat[j, ] =  Yhat.subject[j, ] + as.vector(Z2%*%score2[j, ])
        }
      } #end m
    } #end (if sigma2 > 1e-4)
  } #end if (design == "regular")
  
  
  
  if (design == "irregular") {
    # if sigma2 is very small
    if (sigma2<1e-4) {
      for(m in 1:M) {
        Jm = nVisits[m, 2]  ## number of visits for mth subject
        obs.points = lapply(1:Jm, function(j) which(!is.na(Y.tilde.split[[m]][j, ])))
        Yij.center = lapply(1:Jm, function(j) {Y.tilde.split[[m]][j,obs.points[[j]]]})
        Yi.center = matrix(unlist(Yij.center))
        
        A = t(Z1[unlist(obs.points),])%*%Z1[unlist(obs.points),]
        B = do.call("cbind",lapply(1:Jm, function(j){t(Z1[obs.points[[j]],])%*%Z2[obs.points[[j]],]}))
        C = t(B)
        temp = lapply(1:Jm, function(j){ginv(t(Z2[obs.points[[j]],])%*%Z2[obs.points[[j]],])})
        #invD = as.matrix(bdiag(temp))
        invD = diag.block(temp)
        MatE = ginv(A-B%*%invD%*%C)
        MatF = -invD%*%C%*%MatE
        MatG = -MatE%*%B%*%invD
        MatH = invD - invD%*%C%*%MatG
        Mat1 = cbind(MatE,MatG)
        Mat2 = cbind(MatF,MatH)
        
        
        int1 = t(Z1[unlist(obs.points),])%*%Yi.center
        int2 = do.call('rbind',lapply(1:Jm,function(j){t(t(Z2[obs.points[[j]],])%*%Yij.center[[j]])}))
        int = c(int1,as.vector(t(int2)))
        score1[m,] = Mat1 %*% int
        score2[row.split[[m]],] = matrix(Mat2%*%int,nc=npc[[2]],byrow=T)
        for (j in row.split[[m]]) {
          Yhat.subject[j, ] = mu + eta[ ,visit[j]] + as.vector(Z1%*%score1[m,])
          Yhat[j, ] =  Yhat.subject[j, ] + as.vector(Z2%*%score2[j,])
        }
      } #end m
    } else {
      for(m in 1:M) {
        Jm = nVisits[m, 2]  ## number of visits for mth subject
        obs.points = lapply(1:Jm, function(j) which(!is.na(Y.tilde.split[[m]][j, ])))
        Yij.center = lapply(1:Jm, function(j) {Y.tilde.split[[m]][j,obs.points[[j]]]})
        Yi.center = matrix(unlist(Yij.center))
        
        
        A = t(Z1[unlist(obs.points),])%*%Z1[unlist(obs.points),]/sigma2 + diag(1/evalues[[1]])
        B = do.call("cbind",lapply(1:Jm, function(j){t(Z1[obs.points[[j]],])%*%Z2[obs.points[[j]],]}))/sigma2
        C = t(B)
        temp = lapply(1:Jm, function(j){
          ginv(t(Z2[obs.points[[j]],])%*%Z2[obs.points[[j]],]/sigma2 + diag(1/evalues[[2]]))})
        #invD = as.matrix(bdiag(temp))
        invD = diag.block(temp)
        MatE = ginv(A-B%*%invD%*%C)
        MatF = -invD%*%C%*%MatE
        MatG = -MatE%*%B%*%invD
        MatH = invD - invD%*%C%*%MatG
        Mat1 = cbind(MatE,MatG)
        Mat2 = cbind(MatF,MatH)
        
       
        int1 = t(Z1[unlist(obs.points),])%*%Yi.center
        int2 = do.call('rbind',lapply(1:Jm,function(j){t(t(Z2[obs.points[[j]],])%*%Yij.center[[j]])}))
        int = c(int1,as.vector(t(int2)))/sigma2
        score1[m,] = Mat1 %*% int
        score2[row.split[[m]],] = matrix(Mat2%*%int, nc=npc[[2]], byrow=T)
        for (j in row.split[[m]]) {
          Yhat.subject[j, ] = mu + eta[ ,visit[j]] + as.vector(Z1%*%score1[m,])
          Yhat[j, ] = Yhat.subject[j, ] + as.vector(Z2%*%score2[j, ])
        }
      } #end m
    } # end if (sigma2>1e-4)
  } #end if (design == "irregular")
 
  
  scores <- list(level1 = score1, level2 = score2)
  ###     Return the results from multilevel FPCA as a list
  return(list(Yhat=Yhat, Yhat.subject=Yhat.subject, Y.df=Y.df, scores=scores, mu=mu, 
              efunctions=efunctions, evalues=evalues, sigma2=sigma2, npc=npc, eta=eta))
}



quadWeights<- function(argvals, method="trapezoidal")
{
  ret <- switch(method,
                trapezoidal = {D <- length(argvals)
                1/2*c(argvals[2] - argvals[1], argvals[3:D] -argvals[1:(D-2)], argvals[D] - argvals[D-1])},
                midpoint = c(0,diff(argvals)),
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule"))
  
  return(ret)  
}







