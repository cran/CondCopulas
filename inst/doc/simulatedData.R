## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7,
  fig.height=5
)

## ----setup--------------------------------------------------------------------
library(CondCopulas)
library(VineCopula)
library(ggplot2)
set.seed(1)

## -----------------------------------------------------------------------------
n = 2000

meanX3 = 5
sdX3 = 2
X3 = rnorm(n = n, mean = meanX3, sd = sdX3)

# Computation of conditional parameters
meanX1 = sin(X3)
sdX1 = abs(X3)/2

rateX2 = exp(X3/4)

ckt12_3_fun = function(x){return(1 / (1 + 0.1*(x)^2))}
ckt12_3 = ckt12_3_fun(X3)

ggplot() + geom_line(aes(X3, ckt12_3)) + 
  ggtitle("Conditional Kendall's tau between X1 and X2 conditionally to X3 = x3") +
  xlab("Value of x3") + ylab("Conditional Kendall's tau given X3 = x3")

copFamily12_3 = 3

# Simulation of X1 and X2
X1 = rep(NA, n)
X2 = rep(NA, n)

for (i in 1:n) {
  simCopula = BiCopSim(N=1 , family = copFamily12_3, par = BiCopTau2Par(copFamily12_3 , ckt12_3[i] ))
  X1[i] = qnorm(simCopula[1], mean = meanX1[i], sd = sdX1[i])
  X2[i] = qexp(simCopula[2], rate = rateX2[i])
}


## -----------------------------------------------------------------------------
ggplot() + geom_point(aes(X3, X1))
ggplot() + geom_point(aes(X3, X2))

## -----------------------------------------------------------------------------
newX3 = c(1,4,8,11)
matrixK = computeKernelMatrix(observedX = X3, newX = newX3, kernel = "Gaussian", h = 0.5)

gridXi = seq(0, 1, by = 0.01)
matrixCondQuantilesX1 = estimateCondQuantiles(observedX1 = X1, probsX1 = gridXi, matrixK)
matrixCondQuantilesX2 = estimateCondQuantiles(observedX1 = X2, probsX1 = gridXi, matrixK)
dataCondQuantiles = data.frame(
  pX = gridXi, 
  qX1 = as.numeric(matrixCondQuantilesX1),
  qX2 = as.numeric(matrixCondQuantilesX2),
  valX3 = rep(newX3, each = length(gridXi)))

ggplot(dataCondQuantiles, aes(x = pX, y = qX1, color = valX3, group = as.factor(valX3))) + geom_line() +
  ggtitle("Conditional quantiles of X1 given X3 = x3") +
  xlab("Level of the conditional quantile") +
  ylab("Value of the conditional quantile") + scale_color_continuous(name = "Conditioning\nvalue x3")

ggplot(dataCondQuantiles, aes(x = pX, y = qX2, color = valX3, group = as.factor(valX3))) + geom_line() +
  ggtitle("Conditional quantiles of X2 given X3 = x3") +
  xlab("Level of the conditional quantile") +
  ylab("Value of the conditional quantile") + scale_color_continuous(name = "Conditioning\nvalue x3")


## -----------------------------------------------------------------------------
ggplot() + geom_point(aes(x=X1, y=X2, color = X3)) + 
  scale_color_gradient(low="blue", high="yellow") +
  ggtitle("Scatterplot of (X1, X2) colored by X3")


## -----------------------------------------------------------------------------
whichX3Low = which(X3 <= 2)
whichX3High = which(X3 > 8)
ggplot() + geom_point(aes(x=X1[whichX3Low], y=X2[whichX3Low]))
ggplot() + geom_point(aes(x=X1[whichX3High], y=X2[whichX3High]))


## -----------------------------------------------------------------------------
newX3 = seq(2, 9, by = 0.1)

vecEstimatedThetas = estimateParCondCopula(
  observedX1 = X1, observedX2 = X2, observedX3 = X3, newX3 = newX3, family = 3, h = 0.08)

trueCKT12_3 = ckt12_3_fun(newX3)

comparison = data.frame(
  x3    = rep(newX3, 2),
  param = c(vecEstimatedThetas , BiCopTau2Par(copFamily12_3, trueCKT12_3)),
  CKT   = c(BiCopPar2Tau(copFamily12_3, vecEstimatedThetas) , trueCKT12_3),
  id    = rep(1:2,each=length(newX3)))

ggplot(data = comparison) +
  geom_line(aes(x = x3, y = param, color = as.factor(id), group = as.factor(id))) +
  scale_color_discrete(labels = c("Estimated parameter", "True parameter"), name="") + 
  xlab("Conditioning value x3") +
  ylab(expression(paste("Conditional parameter ",theta,"(x3)"))) + 
  ggtitle(paste0("Comparison between the estimated conditional parameter\n",
                 "of the Clayton copula and the true conditional parameter\n",
                 "as functions of the conditioning variable X3"))

## -----------------------------------------------------------------------------

ggplot(data = comparison, aes(x = x3, y = CKT, color = as.factor(id), group = as.factor(id))) +
  geom_line() +
  scale_color_discrete(labels = c("Estimated CKT","True CKT"), name="") + 
  xlab("Conditioning value x3") +
  ylab("Conditional Kendall's tau") + 
  ggtitle(paste0("Comparison between the estimated conditional Kendall's tau\n",
                 "and the true conditional Kendall's tau\n",
                 "as functions of the conditioning variable X3"))


