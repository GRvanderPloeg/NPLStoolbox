#' Predict Y for new data by projecting the data onto the latent space defined by an NPLS model.
#'
#' @param model NPLS model
#' @param newX New data organized in a matrix of (Inew x J x K) with Inew new subjects
#'
#' @return Ypred: vector of the predicted value(s) of Y for the new data
#' @export
#'
#' @examples
#' Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
#' Ycnt = Y - mean(Y)
#' model = triPLS1(Cornejo2025$Tongue$data, Ycnt, numComponents=1)
#' npred(model, Cornejo2025$Tongue$data[1,,])
npred = function(model, newX){

  I = nrow(model$Fac[[1]])
  J = nrow(model$Fac[[2]])
  K = nrow(model$Fac[[3]])
  numComponents = ncol(model$Fac[[1]])
  numModes = length(model$Fac)
  Fac = model$Fac

  if(length(dim(newX)) == 2){ # Only one sample given
    newX = array(newX, dim=c(1,dim(newX)))
  }
  numSamples = dim(newX)[1]
  unfoldedX = rTensor::k_unfold(rTensor::as.tensor(newX), 1)@data

  W = matrix(0L, nrow=J*K, ncol=numComponents)
  for(f in 1:numComponents){
    W[,f] = pracma::kron(Fac[[3]][,f], Fac[[2]][,f])
  }

  T_scores = matrix(0L, nrow=numSamples, ncol=numComponents)
  for(f in 1:numComponents){
    T_scores[,f] = calculateScores(unfoldedX, Fac[[2]][,f], Fac[[3]][,f])
  }

  Ypred = T_scores %*% model$coef[1:numComponents,1:numComponents] %*% rep(1,numComponents)

  return(Ypred)
}
