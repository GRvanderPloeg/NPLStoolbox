calculateZ = function(X, u, DimX) {

  missing = any(is.na(X))
  J = ncol(X)

  # Find Z matrix
  # If there is missing data, exclude it from calculating covariance.
  if(missing){
    w = rep(NA, J)
    for(j in 1:J){
      mask = !is.na(X[,j])
      if(t(u[mask]) %*% u[mask] != 0){
        ww = t(X[mask,j]) %*% u[mask] / (t(u[mask]) %*% u[mask])
      } else{
        ww = t(X[mask,j]) %*% u[mask]
      }

      if(length(ww) == 0){
        w[j] = 0
      } else{
        w[j] = ww
      }
    }
  } else{
    w = t(X) %*% u
  }

  if(length(DimX) > 2){
    Z = array(w, dim = c(DimX[2], prod(DimX[3:length(DimX)])))
  } else{
    Z = c(w)
  }

  return(Z)
}

calculateScores = function(X, Wj, Wk){

  # if(length(dim(X)) > 2){
  #   X = rTensor::k_unfold(rTensor::as.tensor(X), 1)@data
  #   I = nrow(X)
  # } else{
  #   X = array(X, dim=c(1,dim(X)))
  #   X = rTensor::k_unfold(rTensor::as.tensor(X), 1)@data
  #   I = 1
  # }
  I = dim(X)[1]
  Wkron = kronecker(Wk, Wj)
  t = rep(NA, I)

  if(any(is.na(X))){
    for(i in 1:I){
      mask = !is.na(X[i,])
      t[i] = X[i,mask] %*% Wkron[mask] / (t(Wkron[mask]) %*% Wkron[mask])
    }
  } else{
    t = X %*% Wkron
  }

  return(t)
}

calculateCore = function(X, Fac) {

  LMatTmp = Fac[[1]]
  RMatTmp = t(pracma::kron(Fac[[3]], Fac[[2]]))

  RedData = missmult(pracma::pinv(LMatTmp), X)
  RedData = missmult(RedData, pracma::pinv(RMatTmp))

  G = RedData
  return(G)
}

missmult = function(A, B) {

  ia = nrow(A)
  ja = ncol(A)
  ib = nrow(B)
  jb = ncol(B)

  X = matrix(0L, ia, jb)
  one_array = matrix(1, ia, 1)

  for (j in 1:jb) {
    p = one_array %*% t(B[,j])
    tmpMat = A * p
    X[,j] = t(misssum(t(tmpMat)))
  }

  return(X)
}

misssum = function(X) {

  mask = is.na(X)
  X[mask] = 0

  n_real = length(c(X)) - sum(mask)
  weight = length(c(X))

  # Handle cases where an entire column is missing
  i = which(n_real == 0)
  if(length(i)==0){
    mm = (weight %*% colSums(X)) / n_real
  } else{
    n_real[i] = 1
    mm = weight * sum(X) / n_real
    mm[i] = NaN
  }

  return(mm)
}

safePseudoInverse = function(M, mu=1e-6){
  attempt = try(pracma::pinv(M), silent=TRUE)
  if(inherits(attempt, "try-error")){
    regM = M + mu * diag(1, nrow(M), ncol(M))  # Ensure correct dimensions
    inv = pracma::pinv(regM)
  } else{
    inv = attempt
  }

  return(inv)
}


safeSolve = function(M, mu=1e-6){

  attempt = try(solve(M), silent=TRUE)
  if(inherits(attempt, "try-error")){
    regM = M + mu * diag(ncol(M))
    inv = solve(regM)
  } else{
    inv = attempt
  }

  return(inv)
}

