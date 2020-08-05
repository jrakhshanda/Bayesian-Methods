library(manipulate)
library(latex2exp)

#______(1)______#
BetaPost = function(n, alpha=7, beta=17){
  trueMean = alpha/(alpha+beta)
  trueSD = sqrt((alpha*beta)/((alpha+beta)^2*(alpha+beta+1)))
    Mean = numeric()
    sd = numeric()
    for (i in 1:n) {
      posterior = rbeta(i, alpha, beta)
      Mean[i] = mean(posterior)
      sd[i] = sd(posterior)
    }
    par(mfrow = c(1,2))
    ts.plot(Mean, main = paste('Posterior Mean of Sample Size',n))
    abline(h = trueMean, col = 'red')
    legend('topright', legend = c('Sample Value', 'True Value'),bty = 'n',
           col = c('black','red'), lty = c(1,1),cex=0.7, lwd = 2,xjust = 0)
    
    ts.plot(sd, main = paste('Posterior SD of Sample Size',n))
    abline(h = trueSD, col = 'red')
    legend('topright', legend = c('Sample Value', 'True Value'),bty = 'n',
           col = c('black','red'), lty = c(1,1), cex = 0.7, lwd = 2,xjust = 0)
}

# a,b The hyperparameters in Beta(a,b) prior
# n  number of trials
# p Success proportion in the sample

BetaPriorPostPlot = function(a,b,n,p){
  xGrid = seq(0.001, 0.999, by=0.001)
  normalizedLikelihood = dbeta(xGrid, n*p+1, n*(1-p)+1)
  prior = dbeta(xGrid, a, b)
  posterior = dbeta(xGrid, a+n*p, b+n*(1-p))
  maxDensity = max(normalizedLikelihood, prior, posterior) # Use to make the y-axis high enough
  plot(xGrid, normalizedLikelihood, type = 'l', lwd = 3, xlim = c(0,1), ylim = c(0, maxDensity),
       xlab = expression(theta), ylab = 'Density', 
       main = paste("Bernoulli model - Beta(",a,",",b,") prior"))
  lines(xGrid, posterior, lwd = 2, col = "maroon")
  lines(xGrid, prior, lwd = 2, col = "dodgerblue3")
  legend(x = 0.4, y = maxDensity*0.95, legend = c("Likelihood (normalized)", "Prior", "Posterior"),
         col =  c("black","maroon","dodgerblue3"), lwd = c(3,3,3), cex = 0.8, bty = 'n')
}

BetaPriorPostPlot(a = 2, b = 2, n = 1000, p = 0.25)

#____(2)____#
theta = 0.3
Sample = rbeta(n = 10000, 7, 17)
posterior_prob = mean(Sample>theta)
actual_prob = pbeta(theta, 7,17, lower.tail = F)
#______(3)____#
log_odds = log(Sample/(1-Sample))
hist(log_odds,xlab = 'Log-odds', probability = T,
     main = 'Probability Distribution of log-odds', breaks = 40)
lines(density(log_odds), col = 'red', lwd = 2)
legend('topright', legend = c('Density'), lty = 1, cex = 0.8,
       col = c('red'), lwd = 2, bty = 'n')

