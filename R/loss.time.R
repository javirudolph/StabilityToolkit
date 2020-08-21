#' Function to compute the time (in generations) it takes for a plasmid to be lost
#' or quasi-extinct.  
#'
#' @param beta.o is the initial fraction of plasmid-free cells
#' @param lam is the segregation rate
#' @param sig is the cost
#' @param crit.pct is the quasi-extinction threshold, typically very close to 1
#' @param stop.when is a safe-guard stopping limit:  if the critical threshold has not been reached by 'stop.when' generations, the function stops accumulating generations
#' @export
loss.time <- function(beta.o, lam, sig, crit.pct,stop.when=300){
  
  pct.loss <- beta.o
  gen <- 0 
  while(pct.loss < crit.pct){
    beta.im1 <- pct.loss;
    beta.lk <- (beta.im1*2^(1+sig) + 2*lam*(1-beta.im1))/(beta.im1*2^(1+sig) + 2*(1-beta.im1));
    pct.loss <- beta.lk;
    if(pct.loss >= crit.pct){gen<- gen}else{gen <- gen+1}
    if(gen >= stop.when){
    	print(paste("the critical percent has not been reached by generation ",stop.when,sep=""));
    	break();
    	}
    
  }
  
  return(gen)
}
