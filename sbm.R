update_mean_s2 <- function(z, mean, s2, t) {
	s2 = (s2*(t-2)+(t-1)/t*z^2+(t-1)/t*mean^2-2*(t-1)/t*mean*z)/(t-1)
	mean = (mean*(t-1)+z)/t
	return(c(mean,s2))
}

sbm <- function(oracle, start.left, start.right, end.epoch=Inf, budget=Inf, alpha=0.01, pval.plot=TRUE, trace=TRUE) {
	e = 1
	t.total = 0
	while (e < end.epoch) {
		x = list()
		x$mid = (start.left+start.right)/2
		x$left = 3/4*start.left + 1/4*start.right
		x$right = 1/4*start.left + 3/4*start.right
		
		# perform two steps of preliminary 
		mean = list()
		s2 = list()
		pval = list()
		for (name in c("mid","left","right")) {
			z = c( oracle(x[[name]]), oracle(x[[name]]) )
			mean[[name]] = mean(z)
			s2[[name]] = var(z)
		}
		t = 2
		t.total = t.total+6
		if (pval.plot) pval.record = matrix(0, nrow=0, ncol=3)
		proceed = FALSE
		
		while (TRUE) {
			for (name in c("mid","left","right")) {
				temp = update_mean_s2(oracle(x[[name]]), mean[[name]], s2[[name]], t)
				mean[[name]] = temp[1]
				s2[[name]] = temp[2]
			}
			t.total = t.total+3
			pval$mid = 2*(1-pt(abs(mean$mid/sqrt(s2$mid))*sqrt(t),t-1))
			pval$left = pt(mean$left/sqrt(s2$left)*sqrt(t),t-1)
			pval$right = 1-pt(mean$right/sqrt(s2$right)*sqrt(t),t-1)
			# if total queries exceed budget, early termination
			if (t.total > budget) {
				if ( pval$mid < max(pval$left,pval$right) ) {
					if (mean$mid > 0) return(x$left)
					else return(x$right)
				}
				else return(x$mid)
			}
			if (pval$mid < alpha) {
				if (mean$mid > 0) start.right = x$mid
				else start.left = x$mid
				proceed = TRUE
			}
			else if (max(pval$left,pval$right)<alpha) {
				start.left = x$left
				start.right = x$right
				proceed = TRUE
			}
			if (pval.plot) pval.record = rbind(pval.record, c(pval$left,pval$mid,pval$right))
			t = t+1
			if (proceed) break
		}
		if (trace) {
			cat("Current interval: (",start.left,",",start.right,"). Number of steps: ",t,"\n",sep="")
		}
		if (pval.plot) {
			plot(1:nrow(pval.record), pval.record[,2], type="b", xlab="t", ylab="p value", ylim=c(0,1), pch=16)
			lines(1:nrow(pval.record), pmax(pval.record[,1],pval.record[,3]), type="b", col="red", pch=16)
			abline(h=alpha, col="gray", lty=2)
			legend("topright", legend=c("mid point","quarter points"), pch=c(16,16), col=c("black","red"), inset=0.02, cex=0.75)
			readline(prompt="[next epoch] ")
		}
		e = e+1
	}
	return((start.left+start.right)/2)
}