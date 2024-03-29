
#' Estimation of conditional Kendall's taus by penalized GLM
#'
#' The function \code{CKT.fit.GLM} fits a regression model for the
#' conditional Kendall's tau \eqn{\tau_{1,2|Z}}
#' between two variables \eqn{X_1} and \eqn{X_2}
#' conditionally to some predictors \eqn{Z}.
#' More precisely, this function fits the model
#' \deqn{\tau_{1,2|Z} =
#' 2 * \Lambda( \beta_0 + \beta_1 \phi_1(Z) + ... + \beta_p \phi_p(Z) )}
#' for a link function \eqn{\Lambda},
#' and \eqn{p} real-valued functions \eqn{\phi_1, ..., \phi_p}.
#' The function \code{CKT.predict.GLM} predicts the values of
#' conditional Kendall's tau for some values of the conditioning variable \eqn{Z}.
#'
#' @param datasetPairs the matrix of pairs and corresponding values of the kernel
#' as provided by \code{\link{datasetPairs}}.
#'
#' @param designMatrix the matrix of predictor to be used for the fitting of the model.
#' It should have the same number of rows as the \code{datasetPairs}.
#'
#' @param link link function, can be one of
#' \code{logit}, \code{probit}, \code{cloglog}, \code{cauchit}).
#'
#' @param ... other parameters passed to
#' \code{ordinalNet::\link[ordinalNet]{ordinalNet}()}.
#'
#' @return \code{CKT.fit.GLM} returns the fitted GLM,
#' an object with S3 class \code{ordinalNet}.
#'
#'
#' @references
#' Derumigny, A., & Fermanian, J. D. (2019).
#' A classification point-of-view about conditional Kendall’s tau.
#' Computational Statistics & Data Analysis, 135, 70-94.
#' (Algorithm 2)
#' \doi{10.1016/j.csda.2019.01.013}
#'
#' @seealso See also other estimators of conditional Kendall's tau:
#' \code{\link{CKT.fit.tree}}, \code{\link{CKT.fit.randomForest}},
#' \code{\link{CKT.fit.nNets}}, \code{\link{CKT.predict.kNN}},
#' \code{\link{CKT.kernel}}, \code{\link{CKT.kendallReg.fit}},
#' and the more general wrapper \code{\link{CKT.estimate}}.
#'
#' @examples
#' # We simulate from a conditional copula
#' set.seed(1)
#' N = 400
#' Z = rnorm(n = N, mean = 5, sd = 2)
#' conditionalTau = 2*plogis(-1 + 0.8*Z - 0.1*Z^2) - 1
#' simCopula = VineCopula::BiCopSim(N=N , family = 1,
#'     par = VineCopula::BiCopTau2Par(1 , conditionalTau ))
#' X1 = qnorm(simCopula[,1])
#' X2 = qnorm(simCopula[,2])
#'
#' datasetP = datasetPairs(X1 = X1, X2 = X2, Z = Z, h = 0.07, cut = 0.9)
#' designMatrix = cbind(datasetP[,2], datasetP[,2]^2)
#' fitCKT_GLM <- CKT.fit.GLM(
#'   datasetPairs = datasetP, designMatrix = designMatrix,
#'   maxiterOut = 10, maxiterIn = 5)
#' print(coef(fitCKT_GLM))
#' # These are rather close to the true coefficients -1, 0.8, -0.1
#' # used to generate the data above.
#'
#' newZ = seq(2,10,by = 0.1)
#' estimatedCKT_GLM = CKT.predict.GLM(
#'   fit = fitCKT_GLM, newZ = cbind(newZ, newZ^2))
#'
#' # Comparison between true Kendall's tau (in red)
#' # and estimated Kendall's tau (in black)
#' trueConditionalTau = 2*plogis(-1 + 0.8*newZ - 0.1*newZ^2) - 1
#' plot(newZ, trueConditionalTau , col="red",
#'    type = "l", ylim = c(-1, 1))
#' lines(newZ, estimatedCKT_GLM)
#'
#' @export
#'
CKT.fit.GLM <- function(datasetPairs,
                        designMatrix = datasetPairs[,2:(ncol(datasetPairs)-3),drop=FALSE],
                        link = "logit", ...)
{
  dim_Z = ncol(datasetPairs) - 4

  yMatrix = cbind(datasetPairs[ ,1], 1-datasetPairs[ ,1]) *
    datasetPairs[ , 2 + dim_Z]

  glm1 = ordinalNet::ordinalNet(
    y = yMatrix, x = designMatrix, family = "cumulative", link = link, ...)

  return(glm1)
}


#' Predict the values of conditional Kendall's tau by penalized GLM
#'
#' @param fit result of a call to \code{CKT.fit.GLM}
#'
#' @param newZ new matrix of observations of the conditioning vector \eqn{Z},
#' with the same number of variables and same names as the \code{designMatrix}
#' that was used to fit the GLM.
#'
#' @return \code{CKT.predict.GLM} returns
#' a vector of (predicted) conditional Kendall's taus of the same size
#' as the number of rows of the matrix \code{newZ}.
#'
#' @rdname CKT.fit.GLM
#' @export
#'
CKT.predict.GLM <- function(fit, newZ){

  prediction = 2*stats::predict(fit, newx = newZ, type = "response" )[,1] - 1

  return(prediction)
}


