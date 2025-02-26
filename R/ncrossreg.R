#' Cross-validation of NPLS by classical K-fold CV.
#'
#' This function runs ACMTF-R with cross-validation. A deterministic K–fold partition
#' is used: the subjects are split in order into `cvFolds` groups. For each fold the
#' training set consists of the other folds and the test set is the current fold.
#'
#' @inheritParams triPLS1
#' @param maxNumComponents Maximum number of components to investigate (default 5).
#' @param cvFolds Number of folds to use in the cross-validation. For example, if `cvFolds`
#'                is 5, then the subjects are deterministically partitioned into 5 groups
#'                (each CV iteration uses 4/5 for training and 1/5 for testing). Default: equal to the number of subjects (i.e. jack-knifing).
#'
#' @return A list with two elements:
#'         - **varExp**: a tibble with the variance–explained (for X and Y) per number of components.
#'         - **RMSE**: a tibble with the RMSE (computed over the unified CV prediction vector) per number of components.
#'
#' @export
#'
#' @examples
#' Y = as.numeric(as.factor(Cornejo2025$Tongue$mode1$GenderID))
#' Ycnt = Y - mean(Y)
#' ncrossreg(Cornejo2025$Tongue$data, Ycnt, cvFolds=2)
ncrossreg = function(X, y,
                     maxNumComponents = 5,
                     maxIter = 120,
                     cvFolds = dim(X)[1]) {

  numSubjects = dim(X)[1]
  X = rTensor::as.tensor(X)

  # Create CV folds
  indices = seq_len(numSubjects)
  foldsPartition = split(indices, cut(seq_along(indices), breaks = cvFolds, labels = FALSE))
  uniqueFolds = seq_len(cvFolds)

  ## --- Create Settings Data Frame ---
  settings = expand.grid(numComponents = 1:maxNumComponents,
                         fold = uniqueFolds)
  settings = settings[order(settings$numComponents, settings$fold), ]

  ## --- Run the Parallel Loop Over All Settings ---
  resultsList = vector("list", nrow(settings))
  for (i in 1:nrow(settings)) {
    currentRow = settings[i, ]
    currentComp = currentRow$numComponents
    foldID = currentRow$fold

    testIdx = foldsPartition[[foldID]]
    trainIdx = setdiff(seq_len(numSubjects), testIdx)

    ## Prepare X
    Xtrain_final = list()
    Xtest_final = list()

    Xtrain = rTensor::as.tensor(X@data[trainIdx,,])
    Xtest = rTensor::as.tensor(X@data[testIdx,,])

    # Centering Xtrain
    unfoldedXtrain = rTensor::k_unfold(Xtrain, 1)@data
    means = colMeans(unfoldedXtrain, na.rm=TRUE)
    unfoldedXtrain_cnt = sweep(unfoldedXtrain, 2, means, FUN="-")
    Xtrain_cnt = rTensor::k_fold(unfoldedXtrain_cnt, m=1, modes=Xtrain@modes)

    # Scaling Xtrain
    # unfoldedXtrain = rTensor::k_unfold(Xtrain_cnt, 2)@data
    # stds = apply(unfoldedXtrain, 1, function(x){stats::sd(x, na.rm=TRUE)})
    # unfoldedXtrain_scl = sweep(unfoldedXtrain, 1, stds, FUN="/")
    Xtrain_cnt_scl = Xtrain_cnt #rTensor::k_fold(unfoldedXtrain_cnt, m=2, modes=Xtrain@modes)

    # Use the means and stds to center and scale Xtest as well
    if(length(dim(Xtest)) == 2){ # Only one sample given
      Xtest = rTensor::as.tensor(array(Xtest@data, dim=c(1,dim(Xtest))))
    }

    unfoldedXtest = rTensor::k_unfold(Xtest, 1)@data
    unfoldedXtest_cnt = sweep(unfoldedXtest, 2, means, FUN="-")
    # unfoldedXtest_scl = sweep(unfoldedXtest_cnt, 1, stds, FUN="/")
    Xtest_cnt_scl = rTensor::k_fold(unfoldedXtest_cnt, m=1, modes=Xtest@modes)

    Xtrain = Xtrain_cnt_scl@data
    Xtest = Xtest_cnt_scl@data

    ## Prepare Y
    Ytrain = y[trainIdx]
    Ymean = mean(Ytrain)
    Ytrain = Ytrain - Ymean

    # Fit model
    model = NPLStoolbox::triPLS1(Xtrain, Ytrain, numComponents = currentComp)

    resultsList[[i]] = list(numComponents = currentComp,
                            fold = foldID,
                            testIdx = testIdx,
                            Xtrain = Xtrain,
                            Xtest = Xtest,
                            model = model,
                            Ymean = Ymean)
  }

  ## --- Group the Results by (numComponents, fold) and Select the Best Model --- ##
  # Create a grouping key for each result.
  keys = sapply(resultsList, function(x) paste(x$numComponents, x$fold, sep = "_"))
  resultsByGroup = split(resultsList, keys)

  ## --- Assemble Predictions and Compute RMSE for Each Number of Components --- ##
  RMSE_list = rep(NA, maxNumComponents)
  predictionsByComp = vector("list", maxNumComponents)

  for (comp in 1:maxNumComponents) {
    Ytest_comp = rep(NA, numSubjects)
    Ypred_comp = rep(NA, numSubjects)
    foldsToUse = uniqueFolds

    for (fold in foldsToUse) {
      key = paste(comp, fold, sep = "_")
      if (!is.null(resultsByGroup[[key]])) {
        bestEntry = resultsByGroup[[key]][[1]]

        pred = npred(bestEntry$model, bestEntry$Xtest)
        pred_original = pred + bestEntry$Ymean

        Ypred_comp[bestEntry$testIdx] = pred_original
      }
    }
    predictionsByComp[[comp]] = Ypred_comp
    RMSE_list[comp] = sqrt(mean((y - Ypred_comp)^2))
  }
  RMSEdata = dplyr::tibble(numComponents = 1:maxNumComponents,
                           RMSE = RMSE_list)

  ## --- Fit Full-Data Models to Compute Variance-Explained --- ##
  varExpX = rep(NA, maxNumComponents)
  varExpY = rep(NA, maxNumComponents)
  for (i in 1:maxNumComponents) {
    X_cnt = parafac4microbiome::multiwayCenter(X@data, mode=1)
    X_cnt_scl = parafac4microbiome::multiwayScale(X_cnt, mode=2)
    Ycnt = y - mean(y)

    full_model = NPLStoolbox::triPLS1(X_cnt_scl, Ycnt, numComponents = i)
    varExpX[i] = full_model$varExpX
    varExpY[i]   = full_model$varExpY
  }
  varExpData = dplyr::as_tibble(cbind(numComponents = 1:maxNumComponents,
                                      X = varExpX,
                                      Y = varExpY))

  return(list("varExp" = varExpData,
              "RMSE"   = RMSEdata))
}

# Ugly solution to namespace issues caused by dplyr
RMSE = NULL
Y = NULL
