#' General function to fit either the SS model or the HT model
#'
#' General function to fit either the SS model or the HT model to the data
#' coming from either one or two combined data sets (in which case each data set
#' is allowed to have its own number of replicates). An option is given to   
#' use experimental data coming from a competition assay to inform the cost estimates
#' Parametric bootstrap confidence intervals, the BIC and the minimized negative
#' log-likelihood score is given.  If two data sets are used, the function 
#' computes the BIC of the model that states that a single set of parameters
#' is needed to describe the data and that of the model that states that rather,
#' two distinct sets of parameters are needed to describe both data sets.
#' The number of parametric bootstrap replicates is set to 1000 by default.
#' Finally, an option is given to indicate whether you want to assume that the cost
#' is a true cost, in the sense that in the absence of the plasmid, an individual
#' reproduces faster.  Alternatively, the function can operate under the assumption 
#' that by carrying the plasmid, an individual reproduces faster.
#' @param data1 character string with the object name of the stability series data
#' @param model Which model do you want to fit? model="SS" fits the segregation selection model, model = "HT" fits the horizontal transfer model
#' @param gen The number of generations per day, defaults to 10.
#' @param true.cost Logical.  If TRUE (default), then the sigma is assumed >0, If FALSE, -1<sigma<0
#' @param add.compet.data Logical.  Defaults to FALSE. If TRUE, then you must specify a second object name. This object should be a matrix of data from a competition assay, specified using the argument 'data2'.  The information in this second data set will supplement the stability assay data set (data1) to jointly estimate the model parameters.  Currently only works if model="SS"
#' @param comb Logical.  Defaults to FALSE.  If TRUE, then you must specify a second object name.The data from this second data matrix will be combined with the data from the first data matrix (data1) to jointly estimate the model parameters.  This option works for both, model="SS" and model="HT".
#' @param data2 Character string.  Corresponds to the name of a second data matrix that should contain data from a second stability assay.   
#' @param B  Number of parametric bootstrap iterations.  Defaults to 1000.
#' @export
dynamic.fit <- function(data1,model="SS",gen=10,true.cost=TRUE, add.compet.data=FALSE, comb=FALSE, data2=NULL, B=1000){

	# Reading the data dimensions
	nreps1 <- (ncol(data1)-1)/2;
	day.counts1 <- data1[,1:(nreps1+1)]
	Dmat1 <- data1[,((nreps1+2):ncol(data1))]
	if((add.compet.data==TRUE)|(comb==TRUE)){
		nreps2 <- (ncol(data2)-1)/2;
		day.counts2 <- data2[,1:(nreps2+1)]
		Dmat2 <- data2[,((nreps2+2):ncol(data2))]
	}
	
	
	# Decision tree following the specified options
	if(model=="SS"){
				
		simplefit <- flexss.bootstrappedMLEs(problem.counts=day.counts1, nrepsi=nreps1,
		Dmat=Dmat1, gens=gen,B=1,neg.cost=1-true.cost)
		mles <- simplefit$mles.cis[2,]
		if(true.cost==TRUE){init.parms <- log(mles) - log(1-mles)}else{
			init.parms <- c(log(mles[1:2]) - log(1-mles[1:2]),log(-mles[3])-log(mles[3]+1))
		}
		
		if(comb==TRUE){
			
			print("I'm fittting the SS model by combining two stability assays data sets")
			fit <- pairwise.fitSS(pars=init.parms, dc1=day.counts1, 			
			dc2=day.counts2,r1=nreps1,r2=nreps2,gens=gen,NCols1=Dmat1, NCols2=Dmat2)
			return(fit)
			
			}else if(add.compet.data==TRUE){
				
				print("I am using data from a competition assay to inform the cost estimate")	
				fit <- jointflexssMLEs(data.compet=data2, data.stabexp=data1, 
				nrc=nreps2, nrepexp=nreps1,gen=gen, B=B, neg.cost=1-true.cost)			
				return(fit)
						
			}else{

				print("I am fitting the SS model to a single data set")
				fit <- flexss.bootstrappedMLEs(problem.counts=day.counts1, nrepsi=nreps1,
				Dmat=Dmat1, gens=gen,B=B,neg.cost=1-true.cost)
				return(fit)
			
			}
		
	}else if(model=="HT"){

		simplefit <- HT.bootstrappedMLEs(problem.counts=day.counts1, nrepsi=nreps1,
		Dmat=Dmat1, gens=gen,B=1)
		mles <- simplefit$mles.cis[2,]
		is.any1 <- sum(mles==1)
		if(is.any1>0){
			index <- which(mles==1, arr.ind=TRUE)
			mles[index] <- 0.99
			
		}
		init.parms <- c(log(mles[1:2]) - log(1 - mles[1:2]),log(mles[3]), log(mles[4:5]) - log(1 - mles[4:5]))
		
		if(comb==TRUE){
			
			print("I'm fittting the HT model by combining two stability assays data sets")
			fit <- pairwise.fitHT(pars=init.parms,dc1=day.counts1, 			
			dc2=day.counts2,r1=nreps1,r2=nreps2,gens=gen,NCols1=Dmat1, NCols2=Dmat2)
			return(fit)
			
		}else{

			print("I am fitting the HT model to a single data set")			
			fit <- HT.bootstrappedMLEs(problem.counts=day.counts1, nrepsi=nreps1,
			Dmat=Dmat1, gens=gen,B=B)
			return(fit)
		}
		
	}else{print("You specified a non-existent model, can only do SS or HT");}
	
}