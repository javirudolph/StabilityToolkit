#' Simulates picking colonies data from a serial transfer experiment
#'
#' Simulates picking colonies data from a serial transfer experiment
#' using the solution of the difference equation for the fraction 
#' of plasmid-free cells for the Segregation Selection model.
#' Calls the SS.predict() function. 
#' See De Gelder et al 2004, 2008 and Ponciano et al 2007.
#' @param beta.o Initial fraction of plasmid free cells
#' @param lam Segregation rate (a number between 0 and 1)
#' @param sig plasmid cost. Usually sig>0 but if -1<sig<0, then this parameter represents the benefit of having a plasmid   
#' @param days vector
#' @param gen Number of generations per day
#' @param Dmat A matrix whose number of rows is the number of days that the assay was ran.  The number of columns corresponds to the number of replicated platings. The entries correspond to the total number of colonies picked for testing in a particular day and replicate.
#' @return a vector of same length as 'says' with the fraction of plasmid free cells. The output is a matrix of the number of plasmid-free colonies picked over time. The number of columns of this matrix is equal to the number of specified replicates. The number of replicates is specified implicitly by the number of columns of Dmat.
#' @export
#' @examples
#' mydmat <- matrix(52,nrow=8,ncol=3)
#' mydays <- 0:(nrow(mydmat)-1)
#' SSdatasim(beta.o=0.01, lam=0.003,sig=0.4, days=7, gen=mydays, Dmat=mydmat)
SSdatasim <- function(beta.o, lam,sig, days, gen, Dmat){
	
	beta.lk1 <- SS.predict(beta.o=beta.o, lam=lam,sig=sig,days=days,gen=gen);
	ndays <- length(beta.lk1);
	nreps <- ncol(Dmat);
	sims <- matrix(0,nrow=ndays, ncol=nreps)
	for(i in 1:nreps){
		
		for(j in 1:ndays){
			
			sims[j,i] <- rbinom(n=1, size=Dmat[j,i], prob=beta.lk1[j]);
		}
	}
	return(sims)
}
