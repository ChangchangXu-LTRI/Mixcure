\name{mixcure.penal.profile.CI.nested}
\alias{mixcure.penal.profile.CI.nested}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mixcure.penal.profile.CI.nested
}
\description{
Profile likelihood confidence interval (PLCI) at alpha level (alpha=0.05 as default) for estimated parameter (MLE or FT-PLE). The PLCI endpoints of single parameter theta, are identified as intersects of profile likelihood function of theta and maximum loglikelihood value subtract 1/2 times chi square distribution (df=1) at level alpha. For FT-PLE the reduced model, nested deviance method was applied to determine the dimension of information matrix dim(I) in the penalty.
}
\usage{
mixcure.penal.profile.CI.nested(formula, data, init, pl, apct = 0.05, LRT.pval = F, iterlim=200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
MC model specification in the form of a conventional survival model, i.e., Surv()~. Note the whatever x variables specified in the latency part of the model are assumed to be also regressed on in the incidence part, and vice versa.
}
  \item{data}{
  a data object in the form of dataframe.
  }
  \item{init}{
  a vector of initial values input for the optimization of the vector of parameters specified in the likelihood. Usually, specified in the order of 1) regression parameters for the incidence part, 2) regression parameters for the latency part, and 3) the shape parameter of Weibull hazard function for the latency part.
  }
  \item{pl}{
  if TRUE, uses FT-PL; otherwise uses the usual likelihood for the model.
  }
  \item{iterlim}{
  specifies the number of maximum iterations for model parameter optimization. DEFAULT set to 200.
  }
  \item{apct}{
  specifies the alpha level of confidence interval bounds. DEFAULT set to 0.05.
  }
  \item{LRT.pval}{
  specifies if needs to calculate LRT p-values at the same time. DEFAULT set to FALSE.
  }
}
\details{
see mixcure.penal.est
}
\value{
The mixcure.penal.profile.CI.nested function returns an object of class 'mixcure.plci' that encompasses a list of the followings:
  \item{coefficients}{contains 3 tables of parameter estimates with lower bounds and upper bounds of PLCIs: i) 'CURE' table for regression parameters of the incidence part; ii) 'SURVIVAL' table for regression parameters of the latency part; iii) 'ALPHA' table for shape parameter of the Weibull baseline hazard function in the latency part.}
}
\references{
Bull (2007), Heinze (2013)
}
\author{
Changchang Xu
}

\examples{
# Begin Example

data(ANNbcBMdat1)
# Fit the MC model using maximum likelihoodlizards.
mc.mle.plci1 <- mixcure.penal.profile.CI.nested(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1,init=c(1,-0.1,-10,1,1), pl=F)
# Now the bias-reduced fit:
mc.ple.plci1 <- mixcure.penal.profile.CI.nested(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1,init=c(1,-0.1,-10,1,1), pl=T)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
