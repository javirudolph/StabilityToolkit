% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dynamic.fit.R
\name{dynamic.fit}
\alias{dynamic.fit}
\title{General function to fit either the SS model or the HT model}
\usage{
dynamic.fit(
  data1,
  model = "SS",
  gen = 10,
  true.cost = TRUE,
  add.compet.data = FALSE,
  comb = FALSE,
  data2 = NULL,
  B = 1000
)
}
\arguments{
\item{data1}{character string with the object name of the stability series data}

\item{model}{Which model do you want to fit? model="SS" fits the segregation selection model, model = "HT" fits the horizontal transfer model}

\item{gen}{The number of generations per day, defaults to 10.}

\item{true.cost}{Logical.  If TRUE (default), then the sigma is assumed >0, If FALSE, -1<sigma<0}

\item{add.compet.data}{Logical.  Defaults to FALSE. If TRUE, then you must specify a second object name. This object should be a matrix of data from a competition assay, specified using the argument 'data2'.  The information in this second data set will supplement the stability assay data set (data1) to jointly estimate the model parameters.  Currently only works if model="SS"}

\item{comb}{Logical.  Defaults to FALSE.  If TRUE, then you must specify a second object name.The data from this second data matrix will be combined with the data from the first data matrix (data1) to jointly estimate the model parameters.  This option works for both, model="SS" and model="HT".}

\item{data2}{Character string.  Corresponds to the name of a second data matrix that should contain data from a second stability assay.}

\item{B}{Number of parametric bootstrap iterations.  Defaults to 1000.}
}
\description{
General function to fit either the SS model or the HT model to the data
coming from either one or two combined data sets (in which case each data set
is allowed to have its own number of replicates). An option is given to
use experimental data coming from a competition assay to inform the cost estimates
Parametric bootstrap confidence intervals, the BIC and the minimized negative
log-likelihood score is given.  If two data sets are used, the function
computes the BIC of the model that states that a single set of parameters
is needed to describe the data and that of the model that states that rather,
two distinct sets of parameters are needed to describe both data sets.
The number of parametric bootstrap replicates is set to 1000 by default.
Finally, an option is given to indicate whether you want to assume that the cost
is a true cost, in the sense that in the absence of the plasmid, an individual
reproduces faster.  Alternatively, the function can operate under the assumption
that by carrying the plasmid, an individual reproduces faster.
}
