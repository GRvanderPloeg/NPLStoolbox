#' Tri-PLS1: three-way PLS regressed onto a y vector
#'
#' @param X Centered tensor of independent data
#' @param y Centered dependent variable
#' @param numComponents Number of components to fit
#' @param tol Relative change in loss for the model to converge (default 1e-10).
#' @param maxIter Maximum number of iterations (default 100).
#'
#' @return Model
#' @export
#'
#' @examples
#' set.seed(123)
#' X <- array(rnorm(100 * 5 * 4), dim = c(100, 5, 4))  # Random tensor (100 samples, 5 vars, 4 vars)
#' y <- rnorm(100)  # Random response variable
#' model <- triPLS1(X, y, numComponents = 2)
triPLS1 = function(X, y, numComponents, tol=1e-10, maxIter=100) {

  I = dim(X)[1]  # Number of samples
  J = dim(X)[2]  # Number of second-mode variables
  K = dim(X)[3]  # Number of third-mode variables

  # Initialize storage
  T_scores = matrix(0, I, numComponents)
  WJ = matrix(0, J, numComponents)
  WK = matrix(0, K, numComponents)
  B = matrix(0L, nrow=numComponents, ncol=numComponents)

  Yres = y
  unfoldedX = rTensor::k_unfold(rTensor::as.tensor(X), 1)@data
  Xres = unfoldedX
  Fac = list()

  for (f in 1:numComponents) {
    u = Yres
    t = stats::rnorm(I)
    tgl = t + 2
    iter = 0

    while((norm(t-tgl, "2")/norm(t, "2") > tol) & (iter < maxIter)){
      tgl = t
      iter= iter + 1

      Z = calculateZ(unfoldedX, u, dim(X))

      # Perform SVD on Z
      svdResult = svd(Z)
      Wj = svdResult$u[, 1]
      Wk = svdResult$v[, 1]

      # Compute score vector
      t = calculateScores(unfoldedX, Wj, Wk)

      # Update Y
      Z = calculateZ(Yres, t, dim(Yres))
      Qloadings = Z / norm(Z,"2")
      Qkron = 1
      if(any(is.na(Yres))){
        for(i in 1:I){
          m = !is.na(Yres)
          u[i] = (Yres[i,mask] %*% Qkron[mask]) / (t(Qkron[mask]) %*% Qkron[mask])
        }
      } else{
        u = Yres # %*% Qkron  in NPLStoolbox
      }
    }

    # Store results
    if(f==1){
      Fac[[1]] = as.matrix(t)
      Fac[[2]] = as.matrix(Wj)
      Fac[[3]] = as.matrix(Wk)
      U = as.matrix(u)
    } else{
      Fac[[1]] = cbind(Fac[[1]], t)
      Fac[[2]] = cbind(Fac[[2]], Wj)
      Fac[[3]] = cbind(Fac[[3]], Wk)
      U = cbind(U, as.matrix(u))
    }

    # Calculate core
    Core = calculateCore(unfoldedX, Fac)

    # Calculate B
    B[1:f,f] = pracma::pinv(t(Fac[[1]]) %*% Fac[[1]]) %*% t(Fac[[1]]) %*% U[,f]

    # Calculate Xhat
    # Xhat = parafac4microbiome::reinflateTensor(Fac[[1]], Fac[[2]], Fac[[3]])
    # Xhat = rTensor::k_unfold(rTensor::as.tensor(Xhat), 1)@data
    Wkron = pracma::kron(Fac[[3]], Fac[[2]])
    Xhat = Fac[[1]] %*% Core %*% t(Wkron)

    # Calculate Yhat
    Yhat = Fac[[1]] %*% B[1:f,1:f]# %*% t(Qkron) used in PLStoolbox but always 1
    Yhat = Yhat[,f]

    # Calculate residuals
    Xres = Xres - Xhat
    Yres = Yres - Yhat
  }

  # Unify output for export
  model = list()
  model$Fac = Fac
  model$coef = B

  Wkron = pracma::kron(Fac[[3]], Fac[[2]])
  model$Xhat = Fac[[1]] %*% Core %*% t(Wkron)

  Yhat = Fac[[1]] %*% B[1:numComponents,1:numComponents] # %*% Qkron in PLStoolbox, but this is always 1
  model$Yhat = Yhat[,numComponents]

  model$Core = Core

  mask = !is.na(unfoldedX)
  Xres = unfoldedX - model$Xhat
  model$varExpX = (1 - (sum(Xres[mask]^2) / sum(unfoldedX[mask]^2))) * 100

  mask = !is.na(y)
  model$varExpY = sum(Yhat[mask]^2) / sum(y[mask]^2) * 100

  model$input = X
  model$Xhat = array(Xhat, dim=dim(X))

  Yfactors = list()
  Yfactors[[1]] = U
  Yfactors[[2]] = rep(1, numComponents)
  model$Yfactors = Yfactors

  return(model)
}
