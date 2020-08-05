library(latex2exp)
#______(1)_________#
rInvchisq <- function(n, df, scale) (df*scale)/rchisq(n,df=df)
y = c(44,25,45,52,30,63,19,50,34,67)
n = 10
mu = 3.7
set.seed(123456)
tau2 = sum((log(y)-mu)^2)/n
varPosterior = rInvchisq(n = 100000,df = n,scale = tau2)
hist(varPosterior, breaks = 100, probability = TRUE,xlim = c(0,1),
     main = TeX('Random Posterior Distribution from $\\sigma^2$'))
lines(density(varPosterior), col = 'red', lwd = 2)
legend('topright', legend = 'Postrior Density', col = 'red', lwd = 2, bty = 'n')

#___Theroretical Inv- chi square(n,tau2) Distribution___#
xGrid <- seq(0.001, 0.999, by=0.001)
posterior = dinvchisq(xGrid, n, tau2, log=FALSE)
plot(xGrid, posterior, type = 'l', lwd = 2, xlim <- c(0,1), ylim <- c(0, maxDensity),
     xlab = expression(theta), 
     ylab = 'Density',
     main = TeX(paste("$Inv-\\chi^2(n,\\tau^2)$ Posterior Distribution")))
legend('topright', legend = 'Postrior Density', col = 'black', lwd = 2, bty = 'n')

#______(2)_______#
gini_coef = 2*pnorm(sqrt(varPosterior/2))-1
hist(gini_coef, breaks = 50, probability = TRUE,
     xlab = 'Gini Coefficient',main = '')
lines(density(gini_coef), col = 'red', lwd = 2)
legend('topright', legend = 'Density', col = 'red', lwd = 2, bty = 'n')
#______(3)_______#
CI = quantile(gini_coef, c(0.05,0.95))
hist(gini_coef, breaks = 70, probability = T,freq = F,xlim = c(0,0.6),
     xlab = 'Gini Coefficient',main = '')
abline(v = CI[1], col = 'red', lwd = 3, lty = 3)
abline(v = CI[2], col = 'red', lwd = 3, lty = 3)
legend('topright', legend = 'Quantile(5%,5%)', lwd = 3,
       col = 'red', lty = 3, bty = 'n')
#_____(HPD)______#
library(HDInterval)
HPD_region = hdi(gini_coef, credMass = 0.9)[1:2]
hist(gini_coef, breaks = 70, probability = T,freq = F,xlim = c(0,0.6),
     xlab = 'Gini Coefficient',main = '')
abline(v = HPD_region[1], col = 'red', lwd = 3, lty = 4)
abline(v = HPD_region[2], col = 'red', lwd = 3, lty = 4)
legend('topright', legend = 'HPD Region', lwd = 3,
       col = 'red', lty = 3, bty = 'n')
