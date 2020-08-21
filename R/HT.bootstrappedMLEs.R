#' Maximum Likelihood estimation of the HT model parameters along with Parametric
#' Bootstrap confidence intervals.
#'
#' @param problem.counts is a matrix.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies.  One column per replicated platings.
#' @param nrepis The number of replicated platings. Must be equal to ncol(day.counts)-1
#' @param gen Number of generations per day
#' @param Dmat: A matrix whose number of rows is the number of days that the assay was ran.  The number of columns corresponds to the number of replicated platings. The entries correspond to the total number of colonies picked for testing in a particular day and replicate.
#' @return A list with the BIC, the negative log-likelihood and a matrix with the model MLEs along with 95% Parametric boostrap confidence intervals
#' @export

HT.bootstrappedMLEs <- function(problem.counts, nrepsi,Dmat, gens = 10,B = 1000){
	
	#Dmat <- problem.counts[,(nrepsi+2):(2*nrepsi +1)]
	days <- problem.counts[,1]

	# 1.  First, we define an initial value for the HT model parameters: beta0, lambda, sigma, theta, gamma.
	#     These parameters represent, respectively: initial fraction, segregation, cost, Michaelis-Menten cst, Horiz. Transf. prob.
	init.parmHT <- c((log(0.00000039)-log(1-0.00000039)),log(0.00026)-log(1-0.00026),log(0.146),log(0.5)-log(1-0.5),(log(0.25)-log(1-0.25)));


	is.there0s <- sum(Dmat[1,]==0)>0
	if(is.there0s==1){
	
		zero.ind <- which(Dmat[1,]==0, arr.ind=TRUE)
		all.ind <- 1:nrepsi;
		not.zero.ind <- all.ind[-zero.ind];
		first.counts <- problem.counts[1, 2:(nrepsi+1)];
		counts.4frac <- first.counts[not.zero.ind];
		beta0.guess <- mean(counts.4frac/Dmat[1,not.zero.ind]);
		
	}else{
		beta0.guess <- mean(problem.counts[1, 2:(nrepsi+1)]/Dmat[1,]);
	}

	if(beta0.guess==0){beta0.guess<- 0.00000001};
	
	init.parmHT[1] <- log(beta0.guess)-log(1-beta0.guess);

	#####  First the usual stuff:  ML estimation + Fisher's info CI's

	# ML estimation for the problematic data set. MLES are for beta0, lambda, sigma, theta, gamma
	fit.dc1 <- optim(par=init.parmHT, fn=negloglikeHT, method="Nelder-Mead",day.counts=problem.counts, nreps=nrepsi, gen=gens, D1=Dmat);
	for(i in 1:4){
		fit.dc1 <- optim(par=fit.dc1$par*1.001, fn=negloglikeHT, method="Nelder-Mead",day.counts=problem.counts, nreps=nrepsi, gen=gens, D1=Dmat);
	}
	
	mls1 <- c(1/(1+exp(-fit.dc1$par[1:2])), exp(fit.dc1$par[3]),1/(1+exp(-fit.dc1$par[4:5])));
	mls.orig <- mls1;
	untrsf.mles <- fit.dc1$par

	# nll and BIC
	nll1 <- fit.dc1$val;
	BIC.HT.data1 <- 2*nll1 + 5*log(nrepsi*nrow(problem.counts));

	# Now a parametric bootstrap CI calculation:
	mlesmat <- matrix(0,nrow=B,ncol=5)

	for(i in 1:B){
	
		simdat <- HTdatasim(beta.o=mls.orig[1], lam=mls.orig[2], sig=mls.orig[3], theta=mls.orig[4], gam=mls.orig[5], days=days, gen=gens, Dmat=Dmat)
		simdat.all <- cbind(days,simdat,Dmat)

		fit.dc1 <- optim(par=untrsf.mles, fn=negloglikeHT, method="Nelder-Mead",day.counts=simdat.all, nreps=nrepsi, gen=gens, D1=Dmat);
		for(j in 1:4){
			fit.dc1 <- optim(par=fit.dc1$par*1.001, fn=negloglikeHT, method="Nelder-Mead",day.counts=simdat.all, nreps=nrepsi, gen=gens, D1=Dmat);
		}

		
		mls1 <- c(1/(1+exp(-fit.dc1$par[1:2])), exp(fit.dc1$par[3]),1/(1+exp(-fit.dc1$par[4:5])));
		mlesmat[i,] <- mls1
	}

	cis.mat <- apply(mlesmat,2,FUN=function(x){quantile(x, probs=c(0.025,0.975))})
	mles.cis <- rbind(cis.mat[1,],mls.orig,cis.mat[2,])
	colnames(mles.cis) <- c("Initial fraction", "Segregation", "Cost", "M-M cst.", "Conjugation");
	row.names(mles.cis) <- c("2.5%", "MLE", "97.5%")
	mles.cis
	
	return(list(mles.cis=mles.cis, min.neg.ll = nll1, BIC.HT = BIC.HT.data1, mlesmat = mlesmat));
	
}
