library(ggplot2)
#_______(a)______#
y = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)
dVonMises = function(y, mu=2.39, kappa){
  res = 1
  for (i in y){
    res = res*(exp(kappa*cos(i-mu))/(2*pi*besselI(kappa,0)))
  }
  return(res)
}

kGrid = seq(0, 10, length.out = 1000)
prior = dexp(kappa)
likeli_hood = dVonMises(y, kappa = kGrid)
normal_likeli = (likeli_hood/sum(likeli_hood))/0.005
posterior = (prior*likeli_hood)
posterior = (posterior/sum(posterior))/0.005

plot(x = kGrid, y = posterior,type = 'l',lwd = 2, ylab = 'Density',
     main = 'Von Mises-Wind Direction')
lines(x = kGrid, y = prior, col = 'maroon', lwd = 2)
lines(x = kGrid, y = normal_likeli, col = 'dodgerblue3',lwd = 2)
legend(x = 5, y = 0.8, legend = c("Posterior", "Prior", "Likelihood (normalized)"), 
       col = c("black","maroon","dodgerblue3"), lwd = c(3,3,3), cex = 0.7, bty = 'n')

