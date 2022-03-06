\name{separation.detection}
\alias{separation.detection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
separation.detection
}
\description{
A function for identifying whether or not separation has occurred under MC model.
}
\usage{
separation.detection(fit, nsteps = 30, pl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{mixcure object as the result of mixcure.penal.est function fit.
}
  \item{pl}{if TRUE, uses FT-PL; otherwise uses the usual likelihood for the model.
}
  \item{nsteps}{Starting from iterlim = 1, the mixcure.penal.est is refitted for iterlim = 2, iterlim = 3, . . . , maxit = nsteps. Default value is 30}
}
\details{
Identifies separated cases for object returned from mixcure.penal.est(), by refitting the model. At each iteration
the maximum number of allowed IWLS iterations is fixed starting from 1 to nsteps. For each value
of iterlim (under mixcure.penal.est function), the estimated asymptotic standard errors are divided to the corresponding ones resulted for iterlim=1. Based on the results in Lesaffre & Albert (1989), if the sequence of ratios in any column of the resultant matrix diverges, then separation occurs and the maximum likelihood estimate for the corresponding parameter has value minus or plus infinity.}
\value{
A matrix of dimension nsteps by length of number of parameters, that contains the ratios of the estimated asymptotic standard errors.
}
\references{
Lesaffre (1989)
}
\author{
Changchang Xu
}

\examples{
data(ANNbcBMdat1)
# Fit the MC model using maximum likelihoodlizards.
mc.mle1 <- mixcure.penal.est(Surv(Time, CENS == 1) ~ Her2,data=ANNbcBMdat1,init=c(1,-0.1,-10,1,1), pl=F)
separation.detection(mc.mle1, pl=F, nsteps =30)
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
