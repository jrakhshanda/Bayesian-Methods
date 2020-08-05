library(manipulate)
library(latex2exp)

####################################################################
## Plotting the beta pdf
####################################################################

BetaPlot <- function(a,b){
  xGrid <- seq(0.001, 0.999, by=0.001)
  prior = dbeta(xGrid, a, b)
  maxDensity <- max(prior) # Use to make the y-axis high enough
  plot(xGrid, prior, type = 'l', lwd = 3, col = "blue",
       xlim <- c(0,1), ylim <- c(0, maxDensity), xlab = TeX("$\\theta$"), 
       ylab = 'Density', main = TeX("Beta($\\alpha, \\beta$)"))
}

manipulate(
  BetaPlot(a,b),
  a = slider(1, 10, step=1, initial = 2, label = "The hyperparameter a in Beta(a,b) density"),
  b = slider(1, 10, step=1, initial = 2, label = "The hyperparameter b in Beta(a,b) density")
)

####################################################################
## Plotting the prior-to-posterior mapping for the Bernoulli model.
####################################################################

BetaPriorPostPlot <- function(a,b,n,p){
  xGrid <- seq(0.001, 0.999, by=0.001)
  normalizedLikelihood = dbeta(xGrid, n*p+1, n*(1-p)+1)
  prior = dbeta(xGrid, a, b)
  posterior = dbeta(xGrid, a+n*p, b+n*(1-p))
  maxDensity <- max(normalizedLikelihood, prior, posterior) # Use to make the y-axis high enough
  plot(xGrid, normalizedLikelihood, type = 'l', lwd = 3, col = "blue", xlim <- c(0,1), 
       ylim <- c(0, maxDensity), xlab = TeX("$\\theta$"), 
       ylab = 'Density', main = TeX("Bernoulli Model - Beta($\\alpha, \\beta$) Prior"))
  lines(xGrid, posterior, lwd = 3, col = "red")
  lines(xGrid, prior, lwd = 3, col = "green")
  legend(x = 0, y = maxDensity,bty = 'n',
         legend = c("Likelihood (normalized)", "Prior", "Posterior"), 
         col = c("blue","green","red"), lwd = c(3,3,3), cex = 0.7)
}

manipulate(
  BetaPriorPostPlot(a,b,n,p),
  a = slider(1, 100, step=1, initial = 2, label = "The hyperparameter a in Beta(a,b) prior"),
  b = slider(1, 100, step=1, initial = 2, label = "The hyperparameter b in Beta(a,b) prior"),
  n = slider(1, 1000, step=1, initial = 10, label = "The number of trials, n"),
  p = slider(0, 1, step=0.01, initial = 0.4, label = "Success proportion in the sample")
)

