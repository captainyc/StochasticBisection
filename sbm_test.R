source("sbm.R")

# general parameter setup
sigma = 0.1
x0 = (sqrt(5)-1)/2
n = 1000
budgets = exp(seq(log(50),log(5000),length.out=40))
kappas = seq(1.25,3,0.25)
xs = list()

for (kappa in kappas) {
  g = function(x) rnorm(1,abs(x-x0)^(kappa-1)*sign(x-x0),sigma)
  temp = matrix(0,nrow=n,ncol=length(budgets))
  pb = txtProgressBar(style=3)
  for (k in 1:length(budgets)) {
    for (i in 1:n) {
      temp[i,k] = sbm(g,0,1,budget=budgets[k],pval.plot=FALSE,trace=FALSE)  
    }
    setTxtProgressBar(pb,k/length(budgets))
  }
  close(pb)
  xs[[paste0("kappa_",kappa)]] = temp
}

par(mfrow=c(2,4))
for (kappa in kappas) {
  risk = colMeans(abs(xs[[paste0("kappa_",kappa)]]-x0))
  plot(log(budgets),kappa*log(risk),type="l",ylab="log risk",xlab="budget",main=paste0("kappa=",kappa))
  slope = -kappa/(2*(kappa-1))
  abline(a=kappa*log(risk[1])-slope*log(budgets[1]),b=slope,col="gray",lty=2)
  title(main=paste0("slope=",round(-kappa/(2*(kappa-1)),3)),line=0.5,cex.main=1)
}