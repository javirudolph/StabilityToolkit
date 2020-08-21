#' Joint estimation of the SS model parameters
#'
#' Joint estimation of the SS model parameters assuming that the parameters for two different stability 
#' assays are the same. For a description of the arguments see the 'negloglikeSS.R' documentation
#' Does not work for cases when you assume that the cost is not a true cost, but a benefit because
#' it assumes that the cost is positive.
#' @param pars a vector of the initial fraction, the segregation and the plasmid cost. Both, the initial fraction and the segregation are given in logit scale (remember that logit 'p' is logit(p) = log(p)-log(1-p)). The cost is given in the log scale
#' @param dc1 is a matrix with the days and counts of the first data set.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies.  One column per replicated platings.
#' @param dc2  is a matrix with the days and counts of the second data set.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies.  One column per replicated platings.
#' @param r1 the number of replicated platings per day, for data set 1
#' @param r2 the number of replicated platings per day, for data set 2
#' @param gens the number of generations per day
#' @param NCols1 A matrix whose number of rows is the number of days that the first assay was ran.  The number of columns corresponds to the number of replicated platings. The entries correspond to the total number of colonies picked for testing in a particular day and replicate.
#' @param NCols2 Same as NCols1, but for the second assay's data (dc2)
#' @export

pairwise.fitSS <- function(pars, dc1, dc2, r1,r2, gens, NCols1, NCols2){
	
	fit.dc1 <- optim(par=pars, fn=negloglikeSS, method="Nelder-Mead",day.counts=dc1, nreps=r1, gen=gens, D1=NCols1);
	fit.dc2 <- optim(par=pars, fn=negloglikeSS, method="Nelder-Mead",day.counts=dc2, nreps=r2, gen=gens, D1=NCols2);
	joint.fit <- optim(par=pars, fn=negloglikeSS.comb, method="Nelder-Mead", day.counts1=dc1, day.counts2=dc2, nreps1=r1,nreps2=r2, gen=gens,D1=NCols1, D2=NCols2);

	
	mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), exp(fit.dc1$par[3]));
	mls2 <- c(1/(1+exp(-fit.dc2$par[1])), 1/(1+exp(-fit.dc2$par[2])), exp(fit.dc2$par[3]));
	mls3 <- c(1/(1+exp(-joint.fit$par[1])), 1/(1+exp(-joint.fit$par[2])), exp(joint.fit$par[3]));		

	#hessian.1 <- numDeriv::hessian(func=negloglikeSS, x=fit.dc1$par, method="Richardson", day.counts=dc1, nreps=r1, gen=gens, D1=NCols1 )
	#hessian.2 <- numDeriv::hessian(func=negloglikeSS, x=fit.dc2$par, method="Richardson",day.counts=dc2, nreps=r2, gen=gens, D1=NCols2)
	#hessian.3 <- numDeriv::hessian(func=negloglikeSS.comb, x=joint.fit$par, method="Richardson", day.counts1=dc1, day.counts2=dc2, nreps1=r1,nreps2=r2, gen=gens,D1=NCols1, D2=NCols2)	

	#Fish.inv1 <- 1/diag(hessian.1);
	#Fish.inv2 <- 1/diag(hessian.2);
	#Fish.inv3 <- 1/diag(hessian.3);	

	
	#zalphahalf <- qnorm(p=0.975, mean=0,sd=1)
	#sigs1 <- zalphahalf*sqrt(Fish.inv1)
	#sigs2 <- zalphahalf*sqrt(Fish.inv2)
	#sigs3 <- zalphahalf*sqrt(Fish.inv3)
	
	#CI1.lo <- fit.dc1$par-sigs1;CI1.hi <- fit.dc1$par+sigs1;	
	#CI2.lo <- fit.dc2$par-sigs2;CI2.hi <- fit.dc2$par+sigs2;		
	#CI3.lo <- joint.fit$par-sigs3;CI3.hi <- joint.fit$par+sigs3;	
	
	#CI1.L <- c(1/(1+exp(-CI1.lo[1])), 1/(1+exp(-CI1.lo[2])), exp(CI1.lo[3]));
	#CI2.L <- c(1/(1+exp(-CI2.lo[1])), 1/(1+exp(-CI2.lo[2])), exp(CI2.lo[3]));	
	#CI3.L <- c(1/(1+exp(-CI3.lo[1])), 1/(1+exp(-CI3.lo[2])), exp(CI3.lo[3]));

	#CI1.H <- c(1/(1+exp(-CI1.hi[1])), 1/(1+exp(-CI1.hi[2])), exp(CI1.hi[3]))
	#CI2.H <- c(1/(1+exp(-CI2.hi[1])), 1/(1+exp(-CI2.hi[2])), exp(CI2.hi[3]))	
	#CI3.H <- c(1/(1+exp(-CI3.hi[1])), 1/(1+exp(-CI3.hi[2])), exp(CI3.hi[3]))

	#CI1 <- cbind(CI1.L,CI1.H);CI2 <- cbind(CI2.L,CI2.H);CI3 <- cbind(CI3.L,CI3.H);
	
	nll1 <- fit.dc1$val;		
	nll2 <- fit.dc2$val;		
	nll3 <- joint.fit$val;
	
	BIC.SS.data1 <- 2*nll1 + 3*log(r1*nrow(dc1));
	BIC.SS.data2 <- 2*nll2 + 3*log(r2*nrow(dc2));
	
	BIC.sep   <- 2*nll1 + 2*nll2 + 6*log(r1*nrow(dc1) + r2*nrow(dc2));
	BIC.joint <- 2*nll3 + 3*log(r1*nrow(dc1) + r2*nrow(dc2));		

	Best.model <- "The 2 data sets have different dynamics"
	if(BIC.joint<BIC.sep){Best.model <- "The 2 data sets have the same dynamics"}

	names.pars <- c("Beta0 (initial frac.)", "lambda (segregation)", "sigma (cost)")

	#out <- list(names.pars=names.pars, mls1=mls1, mls2=mls2, joint.mls =mls3, nll1=nll1, nll2=nll2, joint.nll=nll3, BIC.SS.data1=BIC.SS.data1, BIC.SS.data2 = BIC.SS.data2, BIC.sep = BIC.sep, BIC.joint=BIC.joint, Best.model=Best.model,Vcov1 = Fish.inv1, Vcov2=Fish.inv2, Vcov.joint = Fish.inv3, CI.1=CI1, CI.2=CI2, CI.3=CI3);

		out <- list(names.pars=names.pars, mls1=mls1, mls2=mls2, joint.mls =mls3, nll1=nll1, nll2=nll2, joint.nll=nll3, BIC.SS.data1=BIC.SS.data1, BIC.SS.data2 = BIC.SS.data2, BIC.sep = BIC.sep, BIC.joint=BIC.joint, Best.model=Best.model)


	return(out)
					
}
