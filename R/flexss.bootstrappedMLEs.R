#' Maximum Likelihood estimation of the SS model parameters along with Parametric
#' Bootstrap confidence intervals.
#'
#' @param problem.counts is a matrix.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies.  One column per replicated platings.
#' @param nrepis The number of replicated platings. Must be equal to ncol(day.counts)-1
#' @param gen Number of generations per day
#' @param Dmat: A matrix whose number of rows is the number of days that the assay was ran.  The number of columns corresponds to the number of replicated platings. The entries correspond to the total number of colonies picked for testing in a particular day and replicate. 
#' @return A list with the BIC, the negative log-likelihood and a matrix with the model MLEs along with 95% Parametric boostrap confidence intervals
#' @export
#' @examples
#' flexss.bootstrappedMLEs(problem.counts, nrepsi, Dmat, gens, B=1000, neg.cost=TRUE)

flexss.bootstrappedMLEs <- function(problem.counts, nrepsi, Dmat, gens, B=1000, neg.cost=TRUE){
  
  days <- problem.counts[,1]
  init.parmSS <- c(log(0.4)-log(1-0.4),log(0.026)-log(1-0.026),0);
  if(neg.cost==FALSE){init.parmSS <- c(log(0.4)-log(1-0.4),log(0.026)-log(1-0.026), log(0.04)- log(1-0.04));}
  
  is.there0s <- sum(Dmat[1,]==0)>0
  
  if(is.there0s==1){
    
    zero.ind <- which(Dmat[1,]==0, arr.ind=TRUE)
    all.ind <- 1:nrepsi;
    not.zero.ind <- all.ind[-zero.ind];
    first.counts <- problem.counts[1, 2:(nrepsi+1)];
    counts.4frac <- first.counts[not.zero.ind];
    beta0.guess <- mean(counts.4frac/Dmat[1,not.zero.ind]);
    
  }else{
    beta0.guess <- sum(problem.counts[1, 2:(nrepsi+1)]/Dmat[1,])/nrepsi;
  }
  
  if(beta0.guess==0){beta0.guess<- 0.00000001};
  
  init.parmSS[1] <- log(beta0.guess)-log(1-beta0.guess);
  
  #####  First the usual stuff:  ML estimation

  fit.dc1 <- optim(par=init.parmSS, fn=negllSS.flex, method="Nelder-Mead", day.counts=problem.counts, nreps=nrepsi, gen=gens, D1=Dmat,neg.cost=neg.cost);
  for(i in 1:4){
    fit.dc1 <- optim(par=fit.dc1$par*1.0001, fn=negllSS.flex, method="Nelder-Mead",day.counts=problem.counts, nreps=nrepsi, gen=gens, D1=Dmat,neg.cost=neg.cost);
    #print(fit.dc1$val)
  }
  
  if(neg.cost==TRUE){  
  		mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), -1/(1+exp(-fit.dc1$par[3])));
  	}else{
  		mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), exp(fit.dc1$par[3]));  		
  }
  mls.orig <- mls1;
  nll1 <- fit.dc1$val;  
  	
  
  
  # Now a parametric bootstrap CI calculation:
  mlesmat <- matrix(0,nrow=B,ncol=3)
  
  for(i in 1:B){
    
    simdat <- SSdatasim(beta.o=mls.orig[1], lam=mls.orig[2], sig=mls.orig[3], days=days, gen=gens, Dmat=Dmat)
    simdat.all <- cbind(days,simdat,Dmat)
    
    fit.dc1 <- optim(par=init.parmSS, fn=negllSS.flex, method="Nelder-Mead",day.counts=simdat.all, nreps=nrepsi, gen=gens, D1=Dmat, neg.cost=neg.cost);
    
	if(neg.cost==TRUE){  
  		mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), -1/(1+exp(-fit.dc1$par[3])));
  	}else{
  		mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), exp(fit.dc1$par[3]));  		
  	}

    mlesmat[i,] <- mls1
    
  }
  
  cis.mat <- apply(mlesmat,2,FUN=function(x){quantile(x, probs=c(0.025,0.975))})
  mles.cis <- rbind(cis.mat[1,],mls.orig,cis.mat[2,])
  colnames(mles.cis) <- c("Initial fraction", "Segregation", "Cost");
  row.names(mles.cis) <- c("2.5%", "MLE", "97.5%")
  #return(mles.cis)
  
  
  is.low.ok <- sum(mles.cis[2,]>mles.cis[1,])==3
  is.high.ok <- sum(mles.cis[2,]<mles.cis[3,])==3
   
  ci.problem.present <- 2-(is.low.ok + is.high.ok)
  
  if(ci.problem.present==TRUE){
    # Re-doing the parametric bootstrap if CI problem is present
    mlesmat <- matrix(0,nrow=B,ncol=3)
    
    for(i in 1:B){
      
      simdat <- SSdatasim(beta.o=mls.orig[1], lam=mls.orig[2], sig=mls.orig[3], days=days, gen=gens, Dmat=Dmat)
      simdat.all <- cbind(days,simdat,Dmat)
      
	  fit.dc1 <- optim(par=init.parmSS, fn=negllSS.flex, method="Nelder-Mead",day.counts=simdat.all, nreps=nrepsi, gen=gens, D1=Dmat, neg.cost=neg.cost);
      
      	for(j in 1:20){
        	fit.dc1 <- optim(par=fit.dc1$par*1.000001, fn=negllSS.flex, method="Nelder-Mead",day.counts=simdat.all, nreps=nrepsi, gen=gens, D1=Dmat, neg.cost=neg.cost);
        	#print(c(fit.dc1$par,fit.dc1$val))
      	}
      
      
		if(neg.cost==TRUE){  
  			mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), -1/(1+exp(-fit.dc1$par[3])));
  		}else{
  			mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), exp(fit.dc1$par[3]));  		
  		}
      mlesmat[i,] <- mls1
  	}
    
    cis.mat <- apply(mlesmat,2,FUN=function(x){quantile(x, probs=c(0.025,0.975))})
    mles.cis <- rbind(cis.mat[1,],mls.orig,cis.mat[2,])
    colnames(mles.cis) <- c("Initial fraction", "Segregation", "Cost");
    row.names(mles.cis) <- c("2.5%", "MLE", "97.5%")
    
  }
  
  # neg.loglike, BIC and mles+CI's mat
  BIC.SS.data1 <- 2*nll1 + 3*log(nrepsi*nrow(problem.counts));
  
  return(list(nll=nll1, BIC=BIC.SS.data1, mles.cis=mles.cis))
}
