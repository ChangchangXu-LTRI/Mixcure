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
1) For categorical variables of 3 or more categories, pre-treatment for transforming it to binary dummy variables are needed, when being specified in the formula.
}
\value{
The mixcure.penal.est function returns an object of class 'mixcure' that encompasses a list of the followings:
  \item{coefficients}{contains 3 tables of parameter estimates with standard errors, as well as Wald-type test statistics, p values and interval estimations: i) 'CURE' table for regression parameters of the incidence part; ii) 'SURVIVAL' table for regression parameters of the latency part; iii) 'ALPHA' table for shape parameter of the Weibull baseline hazard function in the latency part.}
  \item{cov}{variance-covariance matrix of all the parameters in the order of incidence part, latency part and shape parameter.}
  \item{formula} {formula of fitted model.}
  \item{init} {initial input values for the parameters.}
  \item{data} {model fitted data.}
  \item{loglikelihood} {returns the (penalized) maximum loglikelihood value.}
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

# High event rate, univariate
data(ANNbcBMdat1)
# Fit the MC model using maximum likelihoodlizards.
mc.mle1 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1,init=c(1,-0.1,-10,1,1), pl=F)
# Now the bias-reduced fit:
mc.ple1 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1,init=c(1,-0.1,-10,1,1), pl=T)

mc.mle1$coefficients
mc.ple1$coefficients

# Relatively low event rate, univariate
data(ANNbcBMdat2)
# Fit the MC model using maximum likelihoodlizards.
mc.mle2 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat2,init=c(1,-0.1,-10,1,1), pl=F)
# Now the bias-reduced fit:
mc.ple2 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat2,init=c(1,-0.1,-10,1,1), pl=T)

mc.mle2$coefficients
mc.ple2$coefficients

# Low event rate, 5 variable
data(ANNbcBMdat5)
# Fit the MC model using maximum likelihoodlizards.
mc.mle5 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,data=ANNbcBMdat5,init=c(5,rep(0,5),-10,rep(0,5),0.1), pl=F)
# Now the bias-reduced fit:
mc.ple5 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2 + LuminalA +TN +MENS0 + TUMCT,data=ANNbcBMdat5,init=c(5,rep(0,5),-10,rep(0,5),0.1), pl=T)

mc.mle5$coefficients
mc.ple5$coefficients



}
