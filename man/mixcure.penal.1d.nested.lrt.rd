\name{mixcure.penal.1d.nested.lrt}
\alias{mixcure.penal.1d.nested.lrt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
mixcure.penal.1d.nested.lrt
}
\description{
Likelihood ratio test (LRT) or penalized likelihood ratio test (PLRT) for estimated parameter, that outputs p-values based on asymptotic chi square distribution with df=1. Under nested deviance method, to perform LRT for a single parameter theta of a MC model with information martix of dim(I)=k, for MLE the reduced model assumes theta=0 and dim(I)=k-1; for FT-PLE the reduced model assumes theta=0 and dim(I)=k for the penalty term.
}
\usage{
mixcure.penal.1d.nested.lrt(formula, data, init, pl, loglik, iterlim = 200)
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
  \item{loglik}{
  specifies the (penalized) loglikelihood value of the full model or non-restrected model.
  }
}
\details{
see mixcure.penal.est
}
\value{
The mixcure.penal.1d.nested.lrt function returns an object of class 'mixcure.1d.nested.lrt' that encompasses a list of the followings:
  \item{coefficients}{returns 3 tables, each contains the (penalized) loglikelihood value of the reduced model, the (penalized) loglikelihood ratio in compared to the (penalized) loglikelihood value of the full model and the p-value under LRT or PLRT: i) 'CURE' table for regression parameters of the incidence part; ii) 'SURVIVAL' table for regression parameters of the latency part; iii) 'ALPHA' table for shape parameter of the Weibull baseline hazard function in the latency part.}
}
\references{
Heinze(2020), Firth(1993)
}
\author{
Changchang Xu
}

\examples{
# Begin example

data(ANNbcBMdat1)

mc.mle1 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1,init=c(1,-0.1,-10,1,1), pl=F)
mc.ple1 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1,init=c(1,-0.1,-10,1,1), pl=T)

# Perform LRT for the MC model using maximum likelihoods:
mc.mle.lrt1 <- mixcure.penal.1d.nested.lrt(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1, mc.mle1$loglikelihood, init=c(1,-0.1,-10,1,1), pl=F)

# Now for the bias-reduced model using penalized maximum likelihood:
mc.ple.lrt1 <- mixcure.penal.1d.nested.lrt(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1, mc.ple1$loglikelihood,init=c(1,-0.1,-10,1,1), pl=T)

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
