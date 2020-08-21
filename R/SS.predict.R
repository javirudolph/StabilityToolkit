#' Computes the solution of the difference equation for the fraction 
#' of plasmid-free cells for the Segregation Selection model. 
#'
#' Computes the solution of the difference equation for the fraction 
#' of plasmid-free cells for the Segregation Selection model. 
#' The output is a time series of the fraction of plasmid-free cells
#' See De Gelder et al 2004, 2008 and Ponciano et al 2007.
#' @param beta.o Initial fraction of plasmid free cells
#' @param lam Segregation rate (a number between 0 and 1)
#' @param sig plasmid cost. Usually sig>0 but if -1<sig<0, then this parameter represents the benefit of having a plasmid   
#' @param days vector
#' @param gen Number of generations per day
#' @return a vector of same length as 'says' with the fraction of plasmid free cells
#' @export
#' @examples
#' SS.predict(beta.o=0.01, lam=0.003,sig=0.4, days=7, gen=10)

SS.predict <- function(beta.o, lam,sig, days, gen){
	
	k1<-days;
	lk1 <-length(k1);
	npoints<-(k1[lk1])*gen + 1;
	beta.lk <- rep(0,npoints);
	beta.lk[1]<- beta.o;

	for(i in 2:npoints){
		beta.im1 <- beta.lk[(i-1)];			
		beta.lk[i]<-(beta.im1*2^(1+sig) + 2*lam*(1-beta.im1))/(beta.im1*2^(1+sig) + 2*(1-beta.im1));
    }

	indexvals1<-k1*gen +1;
	beta.lk1<- beta.lk[indexvals1];

	return(beta.lk1)	
	
}
