\name{mixcure.penal.est}
\alias{mixcure.penal.est}
\title{
mixcure.penal.est
}
\description{
Parameter estimation of MC model via a nonlinear minimization function ('nlm'), which utilizes a newton-type algorithm for optimization of all parameters. Optimization is based on either the usual likelihood or the Firth-type penalized likelihood (FT-PL). In the FT-PL method, Jeffry's invariant prior was used for constructing the penalized likelihood.  For hypothesis testing, a simple Wald-type test statistic was employed here, assuming normality of distribution of the estimates. The standard errors are obtained by extracting the diagonals of variance-covariancematrix, calculated by inverting the numeric hessian matrix using nlm.}
\usage{
mixcure.penal.est(formula, data, init, pl, iterlim = 200)
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
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The mixcure.penal.est function returns an object of class 'mixcure' that encompasses a list of the followings:
  \item{coefficients}{contains 3 tables of parameter estimates with standard errors, as well as Wald-type test statistics, p values and interval estimations: i) 'CURE' table for regression parameters of the incidence part; ii) 'SURVIVAL' table for regression parameters of the latency part; iii) 'ALPHA' table for shape parameter of the Weibull baseline hazard function in the latency part.}
  \item{cov}{variance-covariance matrix of all the parameters in the order of incidence part, latency part and shape parameter.}
}
\references{
Firth (1993)
}
\author{
Changchang Xu
}
\note{
1. Supported methods for objects of class of 'mixcure' include:
i) print.mixcure
ii) coef.mixcure
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Begin Example
data(ANNbcBMdat2)
# Fit the MC model using maximum likelihoodlizards.
mc.mle <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat2,init=c(0.5,-0.1,-5,1,0.1), pl=F)
# Now the bias-reduced fit:
mc.ple <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat2,init=c(0.5,-0.1,-5,1,0.1), pl=T)

mc.mle
mc.ple
}
