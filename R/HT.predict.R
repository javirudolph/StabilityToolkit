#'  Computing the fraction over time of the plasmid-free cells in a stability assay
#'
#' Function to compute the fraction over time of the plasmid-free cells in a stability assay
#' assuming that Horizontal transfer is taking place
#' See Ponciano et al 2007.
#' @param beta.o Initial fraction of plasmid free cells
#' @param lam Segregation rate (a number between 0 and 1)
#' @param sig plasmid cost. Usually sig>0 but if -1<sig<0, then this parameter represents the benefit of having a plasmid   
#' @param theta fraction of plasmid free cells at which the horizontal transfer attains half its maximum
#' @param gam Horizontal transfer probability
#' @param days vector
#' @param gen Number of generations per day
#' @export
HT.predict <- function(beta.o, lam,sig, theta, gam, days, gen){
	
	k1<-days;
	lk1 <-length(k1);
	npoints<-(k1[lk1])*gen + 1;
	beta.lk <- rep(0,npoints);
	beta.lk[1]<- beta.o;

	for(i in 2:npoints){
		beta.im1 <- beta.lk[(i-1)];			
		beta.lk[i]<-((1-(gam*(1-beta.im1))/(theta+1-beta.im1))*beta.im1*2^(1+sig) + 2*lam*(1-beta.im1))/(beta.im1*2^(1+sig) + 2*(1-beta.im1));
    }

	indexvals1<-k1*gen +1;
	beta.lk1<- beta.lk[indexvals1];

	return(beta.lk1)	
	
}
