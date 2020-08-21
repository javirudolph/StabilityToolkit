#' Joint Maximum Likelihood estimation of the SS model parameters along with Parametric
#' Bootstrap confidence intervals for competition and stability assay data.
#'
#' @param problem.counts is a matrix.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies.  One column per replicated platings.
#' @param data.compet is a matrix for the competition experiment.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies (one column per replicated platings), followed by the same number of columns with the total number of colonies picked per replicate, per day. 
#' @param nrc The number of replicated platings for the competition experiment. 
#' @param data.stabexp is a matrix for the stability experiment.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies (one column per replicated platings), followed by the same number of columns with the total number of colonies picked per replicate, per day.  
#' @param nrepexp The number of replicated platings for the stability experiment. 
#' @param gen Number of generations per day
#' @return A list with the BIC, the negative log-likelihood and a matrix with the model MLEs along with 95% Parametric boostrap confidence intervals
#' @export
#' @examples
#' jointflexssMLEs(data.compet=compet.data, data.stabexp= exp.data, nrc=3,nrepexp=3, gen=10, B=1000, neg.cost=TRUE)

jointflexssMLEs <- function(data.compet, data.stabexp, nrc,nrepexp, gen, B=1000, neg.cost=TRUE){
  
  	init.parmSS <- c(log(0.4)-log(1-0.4),log(0.026)-log(1-0.026),0);
  	if(neg.cost==FALSE){
  		init.parmSS <- c(log(0.4)-log(1-0.4),log(0.026)-log(1-0.026), log(0.04)- log(1-0.04));
  	}

	day.countsexp <- data.stabexp[,1:(nrepexp+1)];
	D1.1  <- data.stabexp[,((nrepexp+2):(2*nrepexp +1))];
  	days1 <- day.countsexp[,1]
	
	day.countsc <- data.compet[,1:(nrc+1)];
	D1.2  <- data.compet[,((nrc+2):(2*nrc +1))];
	days2 <- data.compet[,1];
  
   
  	is.there0s <- sum(D1.1[1,]==0)>0
  
  	if(is.there0s==1){
    
    	zero.ind <- which(D1.1[1,]==0, arr.ind=TRUE)
    	all.ind <- 1:nrepexp;
    	not.zero.ind <- all.ind[-zero.ind];
    	first.counts <- day.countsexp[1, 2:(nrepexp+1)];
    	counts.4frac <- first.counts[not.zero.ind];
    	beta0.guess <- mean(counts.4frac/D1.1[1,not.zero.ind]);
    
  	}else{
    	beta0.guess <- sum(day.countsexp[1, 2:(nrepexp+1)]/D1.1[1,])/nrepexp;
  	}
  
  	if(beta0.guess==0){beta0.guess<- 0.00000001};
  	init.parmSS[1] <- log(beta0.guess)-log(1-beta0.guess);
  

	# Optimization  
  	fit.dc1 <- optim(par=init.parmSS, fn=joint.negllSS.flex, method="Nelder-Mead", data.compet=data.compet, 
  	data.stabexp= data.stabexp, nrc=nrc,nrepexp=nrepexp, gen=gen,neg.cost=neg.cost);
  	for(i in 1:4){
    	fit.dc1 <- optim(par=fit.dc1$par*1.01, fn=joint.negllSS.flex, method="Nelder-Mead", data.compet=data.compet, 
    	data.stabexp= data.stabexp, nrc=nrc,nrepexp=nrepexp, gen=gen,neg.cost=neg.cost);  
    }
  
  	if(neg.cost==TRUE){  
  			mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), -1/(1+exp(-fit.dc1$par[3])));
  		}else{
  			mls1 <- c(1/(1+exp(-fit.dc1$par[1])), 1/(1+exp(-fit.dc1$par[2])), exp(fit.dc1$par[3]));  		
  	}
  	untrsf.mls <- fit.dc1$par
  	mls.orig <- mls1;
  	nll1 <- fit.dc1$val;  
  	
  
  
  	# Now a parametric bootstrap CI calculation:
  	mlesmat <- matrix(0,nrow=B,ncol=3)
  
  	for(i in 1:B){
    
    	simdat1 <- SSdatasim(beta.o=mls.orig[1], lam=mls.orig[2], sig=mls.orig[3], days=days1, gen=gen, Dmat=D1.1)
    	simdat.all1 <- cbind(days1,simdat1,D1.1)

    	simdat2 <- SSdatasim(beta.o=0.5, lam=mls.orig[2], sig=mls.orig[3], days=days2, gen=gen, Dmat=D1.2)
    	simdat.all2 <- cbind(days2,simdat2,D1.2)

    	fit.dc1 <- optim(par=untrsf.mls, fn=joint.negllSS.flex, method="Nelder-Mead", data.compet=simdat.all2, 
    	data.stabexp= simdat.all1, nrc=nrc,nrepexp=nrepexp, gen=gen,neg.cost=neg.cost);
    
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

  	is.low.ok <- sum(mles.cis[2,]>mles.cis[1,])==3
  	is.high.ok <- sum(mles.cis[2,]<mles.cis[3,])==3
   
  	ci.problem.present <- min((2-(is.low.ok + is.high.ok)),1)
  
  	if(ci.problem.present==TRUE){
    	# Re-doing the parametric bootstrap if CI problem is present
    	mlesmat <- matrix(0,nrow=B,ncol=3)
    
    	for(i in 1:B){
      
    		simdat1 <- SSdatasim(beta.o=mls.orig[1], lam=mls.orig[2], sig=mls.orig[3], days=days1, gen=gen, Dmat=D1.1)
    		simdat.all1 <- cbind(days1,simdat1,D1.1)

    		simdat2 <- SSdatasim(beta.o=0.5, lam=mls.orig[2], sig=mls.orig[3], days=days2, gen=gen, Dmat=D1.2)
    		simdat.all2 <- cbind(days2,simdat2,D1.2)

    		fit.dc1 <- optim(par=untrsf.mls, fn=joint.negllSS.flex, method="Nelder-Mead", data.compet=simdat.all2, 
    		data.stabexp= simdat.all1, nrc=nrc,nrepexp=nrepexp, gen=gen,neg.cost=neg.cost);
            
      		for(j in 1:20){
        		fit.dc1 <- optim(par=fit.dc1$par*1.000001, fn=joint.negllSS.flex, method="Nelder-Mead", 
        		data.compet=simdat.all2, data.stabexp= simdat.all1, nrc=nrc,nrepexp=nrepexp, gen=gen,neg.cost=neg.cost);
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
    
  	}# End if(ci.problem.present==TRUE)
  
  	# neg.loglike, BIC and mles+CI's mat

  	BIC.SS.data1 <- 2*nll1 + 3*log(nrepexp*length(days1) + nrc*length(days2));
  
  	return(list(nll=nll1, BIC=BIC.SS.data1, mles.cis=mles.cis))
}
