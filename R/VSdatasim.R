#' Data simulation under the variable selection model
#'
#' @param lam lambda parameter in the model
#' @param mu mu parameter in the model
#' @param tausq tau square parameter
#' @param days number of days
#' @param x0 day zero
#' @param gen number of generations per day
#' @param Dmat starting point
#'
#' @return
#' @export
#'
#' @examples
#'
#' ndays <-22
#' days.vec <- 0:(ndays-1)
#' nreps <- 3
#' D <- matrix(50, ncol=nreps, nrow=ndays)
#' varselsim <- vs.sim(lam=0.001533, mu=0.011666, tausq=0.158333, days=days.vec, x0 = 0.0066, gen=10, Dmat=D)
#'
#'
vs.sim <- function(lam=0.001533, mu=0.011666, tausq=0.158333,
                   days, x0 = 0.0066, gen=10, Dmat = matrix(50, ncol=nreps, nrow=ndays)){

  # 'days' have to match 'Dmat'
  gens.vec <- days*gen
  lengens  <- length(gens.vec)
  last.gen <- gens.vec[lengens]
  all.gens <- 0:last.gen
  num.gens <- length(all.gens)
  nreps    <- ncol(Dmat)

  xt.mat <- matrix(0, nrow=num.gens,ncol=nreps)

  xt.mat[1,] <- rep(x0,nreps)

  for(i in 2:num.gens){

    im1 <- i-1
    xtm1.vec <- xt.mat[im1,]
    St.vec   <- rnorm(n=nreps,mean=mu, sd=sqrt(tausq))

    numer  <- xtm1.vec*2^(1+St.vec) + 2*lam*(1-xtm1.vec)
    denom  <- xtm1.vec*2^(1+St.vec) + 2*(1-xtm1.vec)
    xt.mat[i,] <- numer/denom

  }

  x4samps <- xt.mat[match(x=gens.vec, table=all.gens),]
  out.samps <- matrix(0,nrow=lengens, ncol=nreps)

  for(i in 1:lengens){

    out.samps[i,] <- rbinom(n=3, size=Dmat[i,], prob=x4samps[i,])

  }
  return(out.samps)
}
