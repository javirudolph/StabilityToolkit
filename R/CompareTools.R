
#' Log likelihood for beta binomial model
#'
#' @param guess
#' @param datamat
#' @param nreps
#'
#' @return
#' @export
#'
#' @examples
#'
negll.bblogist<- function(guess, datamat,nreps=3){

  days <- datamat[,1] # First column should be the 'Days'
  repdays <- rep(days,nreps)
  counts  <- as.vector(datamat[,2:(nreps+1)])
  Ntrials <- as.vector(datamat[,(nreps+2):(nreps*2+1)])

  beta0 <- guess[1];
  beta1 <- guess[2];
  beta2 <- guess[3];
  beta  <- exp(guess[4]);
  p.x   <- 1/(1+ exp(-(beta0+beta1*repdays  + beta2*repdays^2)))
  alpha.x <- (beta*p.x)/(1-p.x)
  llike <- dbbinom(x=counts, size=Ntrials, alpha=alpha.x, beta=beta, log=TRUE)


  #beta  <- exp(guess[3]);
  #p.x   <- 1/(1+ exp(-(beta0+beta1*repdays)))
  #llike <- dbinom(x=counts, size=Ntrials, prob=p.x, log=TRUE)


  negll <- -sum(llike)
  return(negll)
}

#' Simulates data picking from a beta binomial
#'
#' @param mles
#' @param days
#' @param Dmat
#'
#' @return
#' @export
#'
#' @examples
bblogist.sim <- function(mles, days, Dmat){

  beta0 <- mles[1];beta1 <- mles[2]; beta2 <- mles[3]
  beta <- mles[4];
  p.x   <- 1/(1+ exp(-(beta0+beta1*days  + beta2*days^2)))
  alpha.x <- (beta*p.x)/(1-p.x)
  ndays <- lenght(days)
  nreps <- ncol(Dmat)
  data.out <- matrix(0,nrow=ndays,ncol=nreps)
  for(i in 1:ndays){
    data.out[i,] <- rbbinom(n=nreps, size=Dmat[i,], alpha=alpha.x[i], beta=beta)
  }

  out <- cbind(days, data.out,Dmat)
  colnames(out) <- c("Day", paste0("rep",1:nreps), paste0("D", 1:nreps))
  return(out)
}


#' Comparisson between fits
#'
#' @param datamatAnc
#' @param datamatEv
#' @param nrepsAnc
#' @param nrepsEv
#' @param my.guess
#' @param my.method
#' @param plot.Anc
#' @param plot.Ev
#'
#' @return
#' @export
#'
#' @examples
joint.vs.sep.fit <- function(datamatAnc, datamatEv,nrepsAnc=3,nrepsEv=3, my.guess = c(0.5, 0.6, 0.05,log(1.5)),my.method="BFGS",plot.Anc=FALSE, plot.Ev=FALSE){

  #### Function only works if Ancestral and Evolved have the SAME number of replicates
  #### This restriction can be changed but that extra coding is not needed now (Jan 2019)
  # Old my.guess <- c(0.5, 0.6, 0.05,log(3.5))

  # Ancestral data fit
  nAnc <- nrow(datamatAnc)*nrepsAnc
  Anc.mlest <- optim(par=my.guess, fn=negll.bblogist, method=my.method, hessian=TRUE, datamat=datamatAnc,nreps=nrepsAnc)
  Anc.mles <- Anc.mlest$par
  my.hess<- Anc.mlest$hessian
  Fish.Inv <- ginv(my.hess)
  zalphahalf <- qnorm(p=0.975, mean=0, sd=1)
  st.errs <-zalphahalf*sqrt(diag(Fish.Inv))
  low.cis <- Anc.mles - st.errs
  hi.cis <- Anc.mles + st.errs
  AncCIs.mat <- cbind(low.cis, Anc.mles, hi.cis)
  colnames(AncCIs.mat) <- c("2.5%", "MLE", "97.5%")
  #row.names(AncCIs.mat) <- c("Beta0","Beta1")
  #row.names(AncCIs.mat) <- c("Beta0","Beta1","log(Beta)")
  row.names(AncCIs.mat) <- c("Beta0","Beta1", "Beta2","log(Beta)")

  if(plot.Anc==TRUE){
    days <- datamatAnc[,1] # First column should be the 'Days'
    repdays <- rep(days,nrepsAnc)
    counts  <- as.vector(datamatAnc[,2:(nrepsAnc+1)])
    Ntrials <- as.vector(datamatAnc[,(nrepsAnc+2):(nrepsAnc*2+1)])
    cont.days <- seq(0,max(days),by=0.01)
    logp.posx <- Anc.mles[1] + Anc.mles[2]*cont.days + Anc.mles[3]*cont.days^2
    p.x <- 1/(1+exp(-logp.posx))
    plot(repdays, counts/Ntrials, pch=16)
    points(cont.days,p.x,type="l", col="red", lwd=2)
    title("Ancestor data fit")
  }

  # Evolved data fit
  nEv <- nrow(datamatEv)*nrepsEv
  Ev.mlest <- optim(par=my.guess, fn=negll.bblogist, method=my.method, hessian=TRUE, datamat=datamatEv,nreps=nrepsEv)
  Ev.mles <- Ev.mlest$par
  my.hess<- Ev.mlest$hessian
  Fish.Inv <- ginv(my.hess)
  zalphahalf <- qnorm(p=0.975, mean=0, sd=1)
  st.errs <-zalphahalf*sqrt(diag(Fish.Inv))
  low.cis <- Ev.mles - st.errs
  hi.cis <- Ev.mles + st.errs
  EvCIs.mat <- cbind(low.cis, Ev.mles, hi.cis)
  colnames(EvCIs.mat) <- c("2.5%", "MLE", "97.5%")
  #row.names(EvCIs.mat) <- c("Beta0","Beta1")
  #row.names(EvCIs.mat) <- c("Beta0","Beta1","log(Beta)")
  row.names(EvCIs.mat) <- c("Beta0","Beta1", "Beta2","log(Beta)")

  # joint data fit
  datamatJoint <- rbind(datamatAnc,datamatEv)

  Joint.mlest <- optim(par=my.guess, fn=negll.bblogist, method=my.method, hessian=TRUE, datamat=datamatJoint)
  Joint.mles <- Joint.mlest$par
  my.hess<- Joint.mlest$hessian
  Fish.Inv <- ginv(my.hess)
  zalphahalf <- qnorm(p=0.975, mean=0, sd=1)
  st.errs <-zalphahalf*sqrt(diag(Fish.Inv))
  low.cis <- Joint.mles - st.errs
  hi.cis <- Joint.mles + st.errs
  JointCIs.mat <- cbind(low.cis, Joint.mles, hi.cis)
  colnames(JointCIs.mat) <- c("2.5%", "MLE", "97.5%")
  #row.names(JointCIs.mat) <- c("Beta0","Beta1")
  #row.names(JointCIs.mat) <- c("Beta0","Beta1", "log(Beta)")
  row.names(JointCIs.mat) <- c("Beta0","Beta1", "Beta2","log(Beta)")

  # BICs

  BIC.sep <- (2*Anc.mlest$val) + (2*Ev.mlest$val) + 2*length(my.guess)*log(nEv+nAnc)
  BIC.joint <- (2*Joint.mlest$val) + length(my.guess)*log(nEv+nAnc)

  if(BIC.sep< BIC.joint){Best.model="Separate dynamics is best"}else if(BIC.sep>BIC.joint){Best.model="Joint dynamics is best"}

  return(list(AncCIs.mat = AncCIs.mat, EvCIs.mat = EvCIs.mat, JointCIs.mat =JointCIs.mat,BIC.sep=BIC.sep,BIC.joint=BIC.joint,Best.model=Best.model, negll.anc =Anc.mlest$val, negll.ev = Ev.mlest$val, negll.joint = Joint.mlest$val ))

}
