#' Joint estimation of the HT model parameters 
#'
#' Joint estimation of the HT model parameters assuming that the parameters for two different stability 
#' assays are the same. For a description of the arguments see the 'negloglikeHT.R' documentation
#' @param pars a vector of the initial fraction, the segregation, the plasmid cost, the HT constant and the HT probability. All the parameters, except for the cost, have to be specified in logit scale (remember that logit 'p' is logit(p) = log(p)-log(1-p)). The cost is given in the log scale
#' @param dc1 is a matrix with the days and counts of the first data set.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies.  One column per replicated platings.
#' @param dc2  is a matrix with the days and counts of the second data set.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies.  One column per replicated platings.
#' @param r1 the number of replicated platings per day, for data set 1
#' @param r2 the number of replicated platings per day, for data set 2
#' @param gens the number of generations per day
#' @param NCols1 A matrix whose number of rows is the number of days that the first assay was ran.  The number of columns corresponds to the number of replicated platings. The entries correspond to the total number of colonies picked for testing in a particular day and replicate.
#' @param NCols2 Same as NCols1, but for the second assay's data (dc2)
#' @export
pairwise.fitHT <- function(pars, dc1, dc2, r1,r2, gens, NCols1, NCols2){
	
	fit.dc1 <- optim(par=pars, fn=negloglikeHT, method="Nelder-Mead",day.counts=dc1, nreps=r1, gen=gens, D1=NCols1);
	fit.dc2 <- optim(par=pars, fn=negloglikeHT, method="Nelder-Mead",day.counts=dc2, nreps=r2, gen=gens, D1=NCols2);
	joint.fit <- optim(par=pars, fn=negloglikeHT.comb, method="Nelder-Mead", day.counts1=dc1, day.counts2=dc2, nreps1=r1,nreps2=r2, gen=gens,D1=NCols1, D2=NCols2);

	
	mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), exp(fit.dc1$par[3]), 1/(1+exp(-fit.dc1$par[4])),1/(1+exp(-fit.dc1$par[5])));
	mls2 <- c(1/(1+exp(-fit.dc2$par[1])), 1/(1+exp(-fit.dc2$par[2])), exp(fit.dc2$par[3]), 1/(1+exp(-fit.dc2$par[4])),1/(1+exp(-fit.dc2$par[5])));
	mls3 <- c(1/(1+exp(-joint.fit$par[1])), 1/(1+exp(-joint.fit$par[2])), exp(joint.fit$par[3]), 1/(1+exp(-joint.fit$par[4])),1/(1+exp(-joint.fit$par[5])));		

	nll1 <- fit.dc1$val;		
	nll2 <- fit.dc2$val;		
	nll3 <- joint.fit$val;
	
	BIC.HT.data1 <- 2*nll1 + 3*log(r1*nrow(dc1));
	BIC.HT.data2 <- 2*nll2 + 3*log(r2*nrow(dc2));

	
	BIC.sep   <- 2*nll1 + 2*nll2 + 10*log(r1*nrow(dc1) + r2*nrow(dc2));
	BIC.joint <- 2*nll3 + 5*log(r1*nrow(dc1) + r2*nrow(dc2));		

	Best.model <- "The 2 data sets have different dynamics";
	if(BIC.joint<BIC.sep){Best.model <- "The 2 data sets have the same dynamics"}

	names.pars <- c("Beta0 (initial frac.)", "lambda (segregation)", "sigma (cost)", "theta (HT-cst)", "gamma (HT prob.)")

	out <- list(names.pars=names.pars, mls1=mls1, mls2=mls2, joint.mls =mls3, nll1=nll1, nll2=nll2, joint.nll=nll3, BIC.HT.data1=BIC.HT.data1, BIC.HT.data2 = BIC.HT.data2, BIC.sep = BIC.sep, BIC.joint=BIC.joint, Best.model=Best.model);
	
	return(out)
					
}

