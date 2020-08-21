#' Joint likelihood for two different stability assays
#'
#' Joint likelihood assuming that the parameters for two different stability assays are the same
#' For a description of the arguments see the 'negloglikeSS.R' documentation
#' Mostly used internally by the function 'parwise.fitSS.R'
#' @export
negloglikeSS.comb <- function(parm, day.counts1, day.counts2, nreps1,nreps2, gen,D1, D2){
	
	#Variable defs
	beta.o <-1/(1+exp(-parm[1]));#to transform between 0 and 1.
	lam	<- 1/(1+exp(-parm[2]));
	sig	<- exp(parm[3]);
	
	k1<-day.counts1[,1];lk1 <-length(k1);
	npoints1<-(k1[lk1])*gen + 1;
	k2<-day.counts2[,1];lk2 <-length(k2);
	npoints2<-(k2[lk2])*gen + 1;
	npoints<- max(c(npoints1,npoints2));	

	if(length(D1)==1){D1 <- rep(D1,lk1);}
	if(length(D2)==1){D2 <- rep(D2,lk2);}
	
	beta.lk <- rep(0,npoints);
	beta.lk[1]<- beta.o;

	for(i in 2:npoints){
		beta.im1 <- beta.lk[(i-1)];			
		beta.lk[i]<-(beta.im1*2^(1+sig) + 2*lam*(1-beta.im1))/(beta.im1*2^(1+sig) + 2*(1-beta.im1));
    }

	indexvals1<-k1*gen +1;
	beta.lk1<- beta.lk[indexvals1];
	nloglike1<- matrix(0,nrow=lk1, ncol=nreps1);
	
	indexvals2<-k2*gen +1;
	beta.lk2<- beta.lk[indexvals2];
	nloglike2<- matrix(0,nrow=lk2, ncol=nreps2);
	
	ismatD1 <- is.matrix(D1)
	if(ismatD1==TRUE){

			for(i in 1:lk1){
				for(j in 1:nreps1){
					nloglike1[i,j]<- dbinom(x=day.counts1[i,(j+1)],size=D1[i,j],prob=beta.lk1[i],log=T);			
				}
			}	
		
	}else{
			for(i in 1:lk1){
				for(j in 1:nreps1){
					nloglike1[i,j]<- dbinom(x=day.counts1[i,(j+1)],size=D1[i],prob=beta.lk1[i],log=T);			
				}	
			}
	}	

	ismatD2 <- is.matrix(D2)
	if(ismatD2==TRUE){

			for(i in 1:lk2){
				for(j in 1:nreps2){
					nloglike2[i,j]<- dbinom(x=day.counts2[i,(j+1)],size=D2[i,j],prob=beta.lk2[i],log=T);			
				}
			}	
		
	}else{
			for(i in 1:lk2){
				for(j in 1:nreps2){
					nloglike2[i,j]<- dbinom(x=day.counts2[i,(j+1)],size=D2[i],prob=beta.lk2[i],log=T);			
				}	
			}
	}	

	nloglike<- (-1)*(sum(nloglike1)+sum(nloglike2));
	return(nloglike);
	
}
