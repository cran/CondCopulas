

#' Fit Kendall's regression, a GLM-type model for conditional Kendall's tau
#'
#' The function \code{CKT.kendallReg.fit} fits a regression-type model for the
#' conditional Kendall's tau between two variables \eqn{X_1} and \eqn{X_2}
#' conditionally to some predictors Z.
#' More precisely, it fits the model
#' \deqn{\Lambda(\tau_{X_1, X_2 | Z = z}) = \sum_{j=1}^{p'} \beta_j \psi_j(z),}
#' where \eqn{\tau_{X_1, X_2 | Z = z}} is the conditional Kendall's tau
#' between \eqn{X_1} and \eqn{X_2} conditionally to \eqn{Z=z},
#' \eqn{\Lambda} is a function from \eqn{]-1, 1]} to \eqn{R},
#' \eqn{(\beta_1, \dots, \beta_p)} are unknown coefficients to be estimated
#' and \eqn{\psi_1, \dots, \psi_{p'})} are a dictionary of functions.
#' To estimate \eqn{beta}, we used the penalized estimator which is defined
#' as the minimizer of the following criteria
#' \deqn{\frac{1}{2n'} \sum_{i=1}^{n'} [\Lambda(\hat\tau_{X_1, X_2 | Z = z_i})
#' - \sum_{j=1}^{p'} \beta_j \psi_j(z_i)]^2 + \lambda * |\beta|_1,}
#' where the \eqn{z_i} are a second sample (here denoted by \code{ZToEstimate}).
#'
#' @param X1 a vector of \code{n} observations of the first variable \eqn{X_1}.
#'
#' @param X2 a vector of \code{n} observations of the second variable \eqn{X_2}.
#'
#' @param Z a vector of \code{n} observations of the conditioning variable,
#' or a matrix with \code{n} rows of observations of the conditioning vector
#' (if \eqn{Z} is multivariate).
#'
#' @param ZToEstimate the intermediary dataset of observations of \eqn{Z}
#' at which the conditional Kendall's tau should be estimated.
#'
#' @param designMatrixZ the transformation of the \code{ZToEstimate} that
#' will be used as predictors. By default, no transformation is applied.
#'
#' @param h_kernel bandwidth used for the first step of kernel smoothing.
#'
#' @param Lambda the function to be applied on conditional Kendall's tau.
#' By default, the identity function is used.
#'
#' @param Lambda_inv the functional inverse of \code{Lambda}.
#' By default, the identity function is used.
#'
#' @param lambda the regularization parameter. If \code{NULL},
#' then it is chosen by K-fold cross validation.
#' Internally, cross-validation is performed by the function
#' \code{\link{CKT.KendallReg.LambdaCV}}.
#'
#' @param Kfolds_lambda the number of folds used in the cross-validation
#' procedure to choose \code{lambda}.
#'
#' @param h_lambda the smoothing bandwidth used in the cross-validation
#' procedure to choose \code{lambda}.
#'
#' @param newZ the new observations of the conditioning variable.
#'
#' @param l_norm type of norm used for selection of the optimal lambda by cross-validation.
#' \code{l_norm=1} corresponds to the sum of absolute values of differences
#' between predicted and estimated conditional Kendall's tau
#' while \code{l_norm=2} corresponds to the sum of squares of differences.
#'
#' @param ... other arguments to be passed to \code{\link{CKT.kernel}}
#' for the first step (kernel-based) estimator of conditional Kendall's tau.
#'
#'
#' @param observedX1,observedX2,observedZ old parameter names for \code{X1},
#' \code{X2}, \code{Z}. Support for this will be removed at a later version.
#'
#'
#' @return The function \code{CKT.kendallReg.fit} returns
#' a list with the following components:
#' \itemize{
#'     \item \code{estimatedCKT}: the estimated CKT at the new data points \code{newZ}.
#'
#'     \item \code{fit}: the fitted model, of S3 class glmnet
#'     (see \code{glmnet::\link[glmnet]{glmnet}} for more details).
#'
#'     \item \code{lambda}: the value of the penalized parameter used.
#'     (i.e. either the one supplied by the user or
#'     the one determined by cross-validation)
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2020).
#' On Kendall’s regression.
#' Journal of Multivariate Analysis, 178, 104610.
#' \doi{10.1016/j.jmva.2020.104610}
#'
#' @seealso See also other estimators of conditional Kendall's tau:
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.nNets}}, \code{\link{CKT.predict.kNN}},
#' \code{\link{CKT.kernel}}, \code{\link{CKT.fit.GLM}},
#' and the more general wrapper \code{\link{CKT.estimate}}.
#'
#' See also the test of the simplifying assumption that a
#' conditional copula does not depend on the value of the
#' conditioning variable using the nullity of Kendall's regression
#' coefficients: \code{\link{simpA.kendallReg}}.
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 400
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' newZ = seq(2, 10, by = 0.1)
#' estimatedCKT_kendallReg <- CKT.kendallReg.fit(
#'    X1 = X1, X2 = X2, Z = Z,
#'    ZToEstimate = newZ, h_kernel = 0.07)
#'
#' coef(estimatedCKT_kendallReg$fit,
#'      s = estimatedCKT_kendallReg$lambda)
#'
#' # Comparison between true Kendall's tau (in black)
#' # and estimated Kendall's tau (in red)
#' trueConditionalTau = -0.9 + 1.8 * pnorm(newZ, mean = 5, sd = 2)
#' plot(newZ, trueConditionalTau , col="black",
#'    type = "l", ylim = c(-1, 1))
#' lines(newZ, estimatedCKT_kendallReg$estimatedCKT, col = "red")
#'
#'
#'
#' @export
#'
CKT.kendallReg.fit <- function(
  X1 = NULL, X2 = NULL, Z = NULL, ZToEstimate,
  designMatrixZ = cbind(ZToEstimate, ZToEstimate^2, ZToEstimate^3),
  newZ = designMatrixZ,
  h_kernel,
  Lambda = identity, Lambda_inv = identity,
  lambda = NULL, Kfolds_lambda = 10, l_norm = 1, h_lambda = h_kernel, ... ,
  observedX1 = NULL, observedX2 = NULL, observedZ = NULL)
{
  # Back-compatibility code to allow users to use the old "observedX1 = ..."
  env = environment()
  .observedX1X2_to_X1X2(env)
  .observedZ_to_Z(env)

  kernelEstCKT = CKT.kernel(
    X1 = X1, X2 = X2, Z = Z,
    newZ = ZToEstimate, h = h_kernel, ...)

  stopifnot(ncol(designMatrixZ) == ncol(newZ), ncol(designMatrixZ) > 1)
  stopifnot(NROW(ZToEstimate) == nrow(designMatrixZ))

  whichFinite = which( is.finite(kernelEstCKT$estimatedCKT))
  if (is.null(whichFinite)) {
    stop("No kernel estimation successful. ",
         "Maybe h_kernel is too small?")
  }
  fit = glmnet::glmnet(x = designMatrixZ[whichFinite, ],
                       y = Lambda(kernelEstCKT$estimatedCKT[whichFinite]),
                       family = "gaussian")

  if (is.null(lambda)){
    cat("Beginning of the cross-validation\n")
    resultCV <- CKT.KendallReg.LambdaCV(
      X1 = X1, X2 = X2, Z = Z,
      ZToEstimate = ZToEstimate,
      designMatrixZ = designMatrixZ,
      typeEstCKT = 4, Lambda = Lambda, h_lambda = h_lambda, kernel.name = "Epa",
      Kfolds_lambda = Kfolds_lambda, l_norm = l_norm)

    lambda <- resultCV$lambdaCV
  }

  estimatedCKT = CKT.kendallReg.predict(fit = fit, newZ = newZ,
                                        lambda = lambda, Lambda_inv = Lambda_inv)

  return (list(estimatedCKT = estimatedCKT,
               fit = fit, lambda = lambda))
}


#' Predict conditional Kendall's tau using Kendall's regression
#'
#' The function \code{CKT.kendallReg.predict} predicts
#' the conditional Kendall's tau between two variables
#' \eqn{X_1} and \eqn{X_2} given \eqn{Z=z} for some new
#' values of \eqn{z}.
#'
#' @param fit the fitted model, obtained by a call
#' to \code{CKT.kendallReg.fit}.
#'
#' @param newZ the new observations of the conditioning variable.
#'
#' @importFrom glmnet glmnet
#'
#' @return \code{CKT.kendallReg.predict} returns
#' the predicted values of conditional Kendall's tau.
#'
#' @rdname CKT.kendallReg.fit
#' @export
CKT.kendallReg.predict <- function(fit, newZ, lambda = NULL, Lambda_inv = identity){
  return (Lambda_inv(glmnet::predict.glmnet(fit, newx = newZ,
                                            s = lambda, type = "response")))
}


#' Kendall's regression: choice of the penalization parameter by K-folds cross-validation
#'
#' @description
#' In this model, three variables \eqn{X_1}, \eqn{X_2} and \eqn{Z} are observed.
#' We try to model the conditional Kendall's tau between \eqn{X_1} and \eqn{X_2} conditionally
#' to \eqn{Z=z}, as follows:
#' \deqn{\Lambda(\tau_{X_1, X_2 | Z = z})
#' = \sum_{i=1}^{p'} \beta_i \psi_i(z),}
#' where \eqn{\tau_{X_1, X_2 | Z = z}} is the conditional Kendall's tau
#' between \eqn{X_1} and \eqn{X_2} conditionally to \eqn{Z=z},
#' \eqn{\Lambda} is a function from \eqn{]-1, 1[]} to \eqn{R},
#' \eqn{(\beta_1, \dots, \beta_p)} are unknown coefficients to be estimated
#' and \eqn{\psi_1, \dots, \psi_{p'})} are a dictionary of functions.
#' To estimate \eqn{beta}, we used the penalized estimator which is defined
#' as the minimizer of the following criteria
#' \deqn{\frac{1}{2n'} \sum_{i=1}^{n'} [\Lambda(\hat\tau_{X_1, X_2 | Z = z})
#' - \sum_{j=1}^{p'} \beta_j \psi_j(z)]^2 + \lambda * |\beta|_1.}
#' This function chooses the penalization parameter \eqn{lambda}
#' by cross-validation.
#'
#' @param X1 a vector of n observations of the first variable \eqn{X_1}.
#'
#' @param X2 a vector of n observations of the second variable \eqn{X_2}.
#'
#' @param Z a vector of n observations of the conditioning variable,
#' or a matrix with n rows of observations of the conditioning vector
#' (if \eqn{Z} is multivariate).
#'
#' @param ZToEstimate the new data of observations of Z at which
#' the conditional Kendall's tau should be estimated.
#'
#' @param designMatrixZ the transformation of the ZToEstimate that
#' will be used as predictors. By default, no transformation is applied.
#'
#' @param Lambda the function to be applied on conditional Kendall's tau.
#' By default, the identity function is used.
#'
#' @param Kfolds_lambda the number of folds used in the cross-validation
#' procedure to choose \code{lambda}.
#'
#' @param h_lambda the smoothing bandwidth used in the cross-validation
#' procedure to choose \code{lambda}.
#'
#' @param typeEstCKT type of estimation of the conditional Kendall's tau.
#'
#' @param l_norm type of norm used for selection of the optimal lambda.
#' l_norm=1 corresponds to the sum of absolute values of differences
#' between predicted and estimated conditional Kendall's tau
#' while l_norm=2 corresponds to the sum of squares of differences.
#'
#' @param kernel.name name of the kernel. Possible choices are
#' "Gaussian" (Gaussian kernel) and "Epa" (Epanechnikov kernel).
#'
#' @param matrixSignsPairs the results of a call to
#' \code{\link{computeMatrixSignPairs}} (if already computed).
#' If \code{NULL} (the default value), the \code{matrixSignsPairs}
#' will be computed again from the data.
#'
#' @param progressBars should progress bars be displayed?
#' Possible values are
#' \itemize{
#'    \item \code{"none"}: no progress bar at all.
#'
#'    \item \code{"global"}: only one global progress bar (default behavior)
#'
#'    \item \code{"eachStep"}: uses a global progress bar + one progress bar
#'    for each kernel smoothing step.
#' }
#'
#' @param observedX1,observedX2,observedZ old parameter names for \code{X1},
#' \code{X2}, \code{Z}. Support for this will be removed at a later version.
#'
#'
#' @return A list with the following components
#' \itemize{
#'   \item \code{lambdaCV}: the chosen value of the
#'   penalization parameters \code{lambda}.
#'
#'   \item \code{vectorLambda}: a vector containing the values of
#'   \code{lambda} that have been compared.
#'
#'   \item \code{vectorMSEMean}: the estimated MSE for each value of
#'   \code{lambda} in \code{vectorLambda}
#'
#'   \item \code{vectorMSESD}: the estimated standard deviation of the
#'   MSE for each \code{lambda}. It can be used to construct confidence
#'   intervals for estimates of the MSE given by \code{vectorMSEMean}.
#' }
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2020).
#' On Kendall’s regression.
#' Journal of Multivariate Analysis, 178, 104610.
#'
#' @seealso the main fitting function \code{\link{CKT.kendallReg.fit}}.
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 400
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = -0.9 + 1.8 * pnorm(Z, mean = 5, sd = 2)
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' newZ = seq(2, 10, by = 0.1)
#' result <- CKT.KendallReg.LambdaCV(X1 = X1, X2 = X2, Z = Z,
#'                                   ZToEstimate = newZ, h_lambda = 2)
#'
#' plot(x = result$vectorLambda, y = result$vectorMSEMean,
#'      type = "l", log = "x")
#'
#' @export
#'
CKT.KendallReg.LambdaCV <- function(
  X1 = NULL, X2 = NULL, Z = NULL, ZToEstimate,
  designMatrixZ = cbind(ZToEstimate, ZToEstimate^2, ZToEstimate^3),
  typeEstCKT = 4, h_lambda, Lambda = identity, kernel.name = "Epa",
  Kfolds_lambda = 10, l_norm = 1, matrixSignsPairs = NULL,
  progressBars = "global" ,
  observedX1 = NULL, observedX2 = NULL, observedZ = NULL)
{
  # Back-compatibility code to allow users to use the old "observedX1 = ..."
  env = environment()
  .observedX1X2_to_X1X2(env)
  .observedZ_to_Z(env)

  n = length(X1)
  nZprime = length(ZToEstimate)

  if (is.null(matrixSignsPairs)){
    matrixSignsPairs = computeMatrixSignPairs(
      vectorX1 = X1, vectorX2 = X2, typeEstCKT = typeEstCKT)
  }

  if (is.vector(Z)){
    estimCKTNP <- CKT.kernel.univariate
  } else {
    estimCKTNP <- CKT.kernel.multivariate
  }
  switch (progressBars,

          "none" = { globProgressBar = FALSE
          indivProgressBar = FALSE },

          "global" = { globProgressBar = TRUE
          indivProgressBar = FALSE },

          "eachStep" = { globProgressBar = FALSE
          indivProgressBar = TRUE },

          {stop("Uncorrect value for 'progressBars'.")}
  )

  foldid = sample(rep(seq(Kfolds_lambda), length = n))
  list_vectorEstimate = as.list(seq(Kfolds_lambda))
  list_vectorEstimate_comp = as.list(seq(Kfolds_lambda))
  list_resultEstimation = as.list(seq(Kfolds_lambda))
  vectorLambda = c()
  if (globProgressBar) {pb <- pbapply::startpb(min = 0, max = 3 * Kfolds_lambda)}
  for (i in seq(Kfolds_lambda)) {
    which = foldid == i
    list_vectorEstimate[[i]] = estimCKTNP(
      matrixSignsPairs = matrixSignsPairs[!which, !which],
      Z = Z[!which], typeEstCKT = typeEstCKT,
      h = h_lambda, ZToEstimate = ZToEstimate, kernel.name = kernel.name,
      progressBar = indivProgressBar)

    if (globProgressBar) {pbapply::setpb(pb, pbapply::getpb(pb) + 1)}

    list_vectorEstimate_comp[[i]] = estimCKTNP(
      matrixSignsPairs = matrixSignsPairs[which, which],
      Z = Z[which], typeEstCKT = typeEstCKT,
      h = h_lambda, ZToEstimate = ZToEstimate, kernel.name = kernel.name,
      progressBar = indivProgressBar)

    if (globProgressBar) {pbapply::setpb(pb, pbapply::getpb(pb) + 1)}

    whichna = which(!is.finite(list_vectorEstimate[[i]]))
    if (length(which(is.finite(list_vectorEstimate[[i]]))) == 0) {
      stop("Unable to estimate the conditional Kendall's tau in chooseLambdaCV. ",
           "Possible explanation: h_lambda and/or the sample size are too small, ",
           "or Kfolds_lambda is too large.")
    } else if (length(whichna) > 0) {
      list_resultEstimation[[i]] =
        glmnet::glmnet(y = Lambda(list_vectorEstimate[[i]][-whichna]),
                       x = designMatrixZ[-whichna,])
    } else {
      list_resultEstimation[[i]] =
        glmnet::glmnet(y = Lambda(list_vectorEstimate[[i]]),
                       x = designMatrixZ)
    }
    vectorLambda = c(vectorLambda, list_resultEstimation[[i]]$lambda)

    if (globProgressBar) {pbapply::setpb(pb, pbapply::getpb(pb) + 1)}
  }
  if (globProgressBar) {pbapply::closepb(pb)}

  # Choice of the lambda used
  uniqueLambda = unique(vectorLambda)
  vectorLambda = sort(uniqueLambda[which( (1:length(vectorLambda))%%3 == 0 ) ])
  nLambda = length(vectorLambda)
  matrixDiffLambda = matrix(nrow = Kfolds_lambda, ncol = nLambda)
  for (i in seq(Kfolds_lambda)) {
    whichna = which(!is.finite(list_vectorEstimate_comp[[i]]))
    if (length(whichna) > 0)
    {
      matrix_prediction_tau_CV =
        matrix(rep(list_vectorEstimate_comp[[i]][-whichna], nLambda),
               nrow = nZprime - length(whichna), ncol = nLambda)
      matrix_prediction_tau_lambda =
        glmnet::predict.glmnet(object = list_resultEstimation[[i]] ,
                               newx = designMatrixZ[-whichna,], s = vectorLambda)
      matrixDiffLambda[i,] = as.vector(
        rep(1, nZprime - length(whichna)) %*%
          abs( Lambda(matrix_prediction_tau_CV) - matrix_prediction_tau_lambda ) )
    } else {
      matrix_prediction_tau_CV = matrix(rep(list_vectorEstimate_comp[[i]], nLambda),
                                        nrow = nZprime, ncol = nLambda)
      matrix_prediction_tau_lambda =
        glmnet::predict.glmnet(object = list_resultEstimation[[i]] ,
                               newx = designMatrixZ, s = vectorLambda)
      matrixDiffLambda[i,] = as.vector(
        rep(1, nZprime) %*%
          abs( Lambda(matrix_prediction_tau_CV) - matrix_prediction_tau_lambda ) )
    }
  }
  if (l_norm == 2)
  {
    matrixDiffLambda = matrixDiffLambda^2
  }
  vectorMSEMean = apply(matrixDiffLambda, MARGIN = 2, FUN = mean, na.rm = TRUE)
  vectorMSESD = apply(matrixDiffLambda, MARGIN = 2, FUN = stats::sd, na.rm = TRUE)
  lambdaCV = vectorLambda[which.min(vectorMSEMean)]

  return(list(lambdaCV = lambdaCV, vectorLambda = vectorLambda,
              vectorMSEMean = vectorMSEMean, vectorMSESD = vectorMSESD) )
}



