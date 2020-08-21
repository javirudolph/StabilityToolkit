#' Computes the negative log-likelihood for the Horizontal Transfer model
#'
#' Computes the negative log-likelihood for the Horizontal Transfer model, under binomial sampling 
#' See Ponciano et al 2007.
#' @param parm=c(beta.o,lam,sig, theta, gamma) which are the initial fraction of plasmid free cells, segregation rate and the plasmid cost. All three parameters must be specified according to the inverse transform of either the 1/(1+e-x) or the exp(x) functions. See example below
#' @param day.counts is a matrix.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies.  One column per replicated platings.
#' @param nreps The number of replicated platings. Must be equal to ncol(day.counts)-1
#' @param gen Number of generations per day
#' @param D1: Either a matrix or a vector.  If a matrix, its number of rows is the number of days that the assay was ran.  The number of columns corresponds to the number of replicated platings. The entries correspond to the total number of colonies picked for testing in a particular day and replicate. If a vector, it is assumed that the daily number of colonies picked for each replicate is the same.
#' @return the negative log-likelihood
#' @export

negloglikeHT <- function(parm, day.counts, nreps, gen,D1){
	
	beta.o <-1/(1+exp(-parm[1]));#to transform between 0 and 1.
	lam	<- 1/(1+exp(-parm[2]));
	sig	<- exp(parm[3]);
	theta<- 1/(1+exp(-parm[4]));
	gamma<- 1/(1+exp(-parm[5]));
	
	k1<-day.counts[,1];lk1 <-length(k1);
	npoints1<-(k1[lk1])*gen + 1;
	npoints<- npoints1;   #max(c(npoints1,npoints2));	

	if(length(D1)==1){D1 <- rep(D1,lk1);}

	beta.lk <- rep(0,npoints);
	beta.lk[1]<- beta.o;

	for(i in 2:npoints){
		beta.im1 <- beta.lk[(i-1)];			
		beta.lk[i]<-((1-(gamma*(1-beta.im1))/(theta+1-beta.im1))*beta.im1*2^(1+sig) + 2*lam*(1-beta.im1))/(beta.im1*2^(1+sig) + 2*(1-beta.im1));
    }

	indexvals1<-k1*gen +1;
	beta.lk1<- beta.lk[indexvals1];
	nloglike1<- matrix(0,nrow=lk1, ncol=nreps);
	
	ismatD <- is.matrix(D1)
	if(ismatD==TRUE){

			for(i in 1:lk1){
				for(j in 1:nreps){
					nloglike1[i,j]<- dbinom(x=day.counts[i,(j+1)],size=D1[i,j],prob=beta.lk1[i],log=T);			
				}
			}	
		
	}else{
			for(i in 1:lk1){
				for(j in 1:nreps){
					nloglike1[i,j]<- dbinom(x=day.counts[i,(j+1)],size=D1[i],prob=beta.lk1[i],log=T);			
				}	
			}
	}	

	nloglike<- (-1)*(sum(nloglike1))#+sum(nloglike2)+sum(nloglike3));
	return(nloglike);
	
}
