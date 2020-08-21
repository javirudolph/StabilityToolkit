#' joint likelihod from a competition experiment and initial fraction of 0.5 and from a stability assay
#'
#' Uses data from a competition experiment starting with an initial fraction of 0.5 and from a stability assay
#' to compute the joint negative log-likelihood for the Segregation Selection model, under binomial sampling 
#' @param parm=c(beta.o,lam,sig) which are the initial fraction of plasmid free cells, segregation rate and the plasmid cost. All three parameters must be specified according to the inverse transform of either the 1/(1+e^(-x)) or the -1/(1+e^(-x)) functions.
#' @param data.compet is a matrix for the competition experiment.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies (one column per replicated platings), followed by the same number of columns with the total number of colonies picked per replicate, per day. 
#' @param nrc The number of replicated platings for the competition experiment. 
#' @param data.stabexp is a matrix for the stability experiment.  The first column are the days, the next columns are the counts over time of the plasmid-free colonies (one column per replicated platings), followed by the same number of columns with the total number of colonies picked per replicate, per day.  
#' @param nrepexp The number of replicated platings for the stability experiment. 
#' @param gen Number of generations per day
#' @param neg.cost Logical. If TRUE, then assumes cost is between -1 and 0.
#' @return the negative log-likelihood
#' @export
#' @examples
#' joint.negllSS.flex(parm, data.compet, data.stabexp,nrc,nrepexp, gen,neg.cost=TRUE)

joint.negllSS.flex <- function(parm, data.compet, data.stabexp,nrc,nrepexp, gen,neg.cost=TRUE){
	
	day.countsexp <- data.stabexp[,1:(nrepexp+1)];
	D1.1  <- data.stabexp[,((nrepexp+2):(2*nrepexp +1))];
	
	negllSS.1 <- negllSS.flex(parm=parm, day.counts=day.countsexp, nreps=nrepexp, gen=gen,D1=D1.1,neg.cost=neg.cost)
	
	day.countsc <- data.compet[,1:(nrc+1)];
	D1.2  <- data.compet[,((nrc+2):(2*nrc +1))];
	
	
	negllSS.2 <- negllSS.flex.nobeta0(parm=parm[2:3], day.counts=day.countsc, nreps=nrc, gen=gen,D1=D1.2,neg.cost=neg.cost)
	
	return(negllSS.1 + negllSS.2)	
	
}
