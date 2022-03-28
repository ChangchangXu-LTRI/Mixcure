#####################################################################
### Function for ikelihood ratio test (df=1) mixcure model         ##
#### via the nested deviance method under penalized loglikelihoods ##
#####################################################################
###################################################
#### previously 'mixcure.penal.ESTV.nested1d.r'  ##
###################################################

mixcure.penal.1d.nested.lrt <- function(formula, data, init, pl, loglik, iterlim = 200) {
require(splines)
require(survival)
require(abind)

  #########################################################################################
  mat.inv <- function(matx) {

    detm = det(matx)
    #2x2 matrix inverse;
    if (ncol(matx) == 2) {
      inv.matx = (1/detm) * matrix(c(matx[2,2],-matx[2,1],-matx[1,2],matx[1,1]), nrow = 2)
    }

    else {
      #For any n>2 dimension square matrix;
      adjug.matx <- matrix(rep(0, ncol(matx)^2), nrow = nrow(matx))
      for (i in 1:nrow(matx)) {
        for (j in 1:ncol(matx)) {
          adjug.matx[i,j] <- (-1)^(i+j)*det(matx[-i,][,-j])
        }
      }
      inv.matx <- t(adjug.matx/detm)
    }

    return(inv.matx)
  }

  #########################################################################################

  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];

  design.matrix <- model.matrix(formula, data = design.matrix);

  # index ranges of coefficients of glm and cox models
  index.cure.v <- 1 : ncol(design.matrix);
  index.surv.v <- (ncol(design.matrix) + 1) : (2*length(index.cure.v))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.v)+1;

  samp.s <- nrow(design.matrix)


#############################################################
#############################################################
#############################################################

  #create CI for profile likelihood, this option only outputs estimates and PL under specified model;


  ##loglik function for testing parameters of cure or surv part;
  loglik.mixture.part <- function(p, survt, design.matrix1, design.matrix0,
                                  index.cure.var=index.cure.v,
                                  index.surv.var=index.surv.v, pl) {  #design.matrix1-surv, design.matrix0-cure

    design.mtx.comb = cbind(design.matrix0,design.matrix1)

    #parameter and variable dep parameters;
    if (k <= length(index.cure.v)) {
    theta = 1/(1+exp(-design.mtx.comb[,index.cure.var]%*%as.matrix(p[index.cure.v[-length(index.cure.v)]])))
    eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.var-1]))
    } else {
      theta = 1/(1+exp(-design.mtx.comb[,index.cure.var]%*%as.matrix(p[index.cure.var])))
      eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.v[-length(index.surv.v)]]))
    }
     eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta)
    pi = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi


    max.len = max(length(index.cure.var),length(index.surv.var))
    n.elema = max.len^2
    a.sub1 <- matrix(rep(0,n.elema), nrow = max.len)
    a.sub2 <- matrix(rep(0,n.elema), nrow = max.len)

    for (i in c(index.cure.v)) {
      for (j in c(index.cure.v)) {
        a.sub1[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*theta*(1-theta))[survt[, 2] == 1])
        a.sub2[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*kap)[survt[, 2] == 0])
      }
    }
    info.a = (a.sub1 + a.sub2)


    ##info matrix block B
    design.xt0 <- cbind(design.matrix0, log(survt[,1]))
    n.elemb <- max.len*(max.len+1)
    b.sub <- matrix(rep(0,n.elemb), nrow = max.len)

    for (i in c(index.cure.v)) {
      for (j in c(1:length(index.surv.v), max.len+1)) {
        b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])
      }
    }
    info.b = b.sub

    design.xt1 <- cbind(design.matrix1, log(survt[,1]))

    n.elemd <- (max.len+1)^2
    d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
    d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))

    for (i in c(index.cure.v, max.len +1)) {
      for (j in c(index.cure.v, max.len +1)) {
        d.sub1[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*eps)[survt[, 2] == 1])
        d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
      }
    }
    d.sub = d.sub1 + d.sub2 +
      matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-1]^2)),
             nrow = (max.len + 1))

    info.d = d.sub


    info.d.inv = mat.inv(info.d)

    #fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
    #hessian.mat = -fisher.info

    # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
    info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)

    #determinant of hessian matrix;
    det.info = det(info.set0)*det(info.d)

    #calculate loglikelihood for the unpenalized;
    cure.par <- p[index.cure.var];
    surv.par <- p[index.surv.var];
    p.gamma <- p[index.gamma-1];  #use original shape parameter instead of exp();

    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) -
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );


    if (pl == FALSE)
    {
     loglik.part = loglikelihood
    }
    else if (pl == TRUE)
    {
     loglik.part = loglikelihood - 0.5*log(det.info)
    }

    return(loglik.part)
  }


  #################################################################

  #### parameter estimation under H0 for individual parameter
  #### loglikelihood ratio test statistics for each cure part variable;

  dim.v <- ncol(design.matrix)
  ll.cure <- rep(0,dim.v)
  llr.cure <- rep(0,dim.v)
  pval.cure <- rep(0,dim.v)
# index.cure.v[-1] for no intercept calculation
  for (k in index.cure.v) {
    maximizer <- nlm(
      f = loglik.mixture.part, p = init[-k],
      survt = survt, design.matrix0 = design.matrix,
      design.matrix1=design.matrix,
      index.cure.var=index.cure.v[-k], pl=pl,
      iterlim = iterlim, hessian=T
    );
    loglik.part = -maximizer$minimum;
    dif.ll = -2*(loglik.part-loglik);  #loglik is ll under Ha;
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.cure[k]<- loglik.part
    llr.cure[k]<- dif.ll
    pval.cure[k]<- pval
    if (det(maximizer$hessian) < 1e-05)
      diag(maximizer$hessian) <- diag(maximizer$hessian) + 1e-06

  }



  ### loglikelihood calculation for each surv part variable;

  ll.surv <- rep(0,ncol(design.matrix))
  llr.surv <- rep(0,ncol(design.matrix))
  pval.surv <- rep(0,ncol(design.matrix))


  for (k in index.surv.v) {
    is=k-length(index.cure.v)
    maximizer <- nlm(
      f = loglik.mixture.part, p = init[-k],
      survt = survt, design.matrix1 = design.matrix,
      design.matrix0=design.matrix,
      index.surv.var=index.surv.v[-is], pl=pl,
      iterlim = iterlim, hessian=FALSE
    );

    loglik.part = -maximizer$minimum;
    dif.ll = -2*(loglik.part-loglik);
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.surv[is]<- loglik.part
    llr.surv[is]<-dif.ll
    pval.surv[is]<-pval
  }



  coef.table.cure <- cbind(
    'LL.cure' = ll.cure,
    'LLR'         = llr.cure,
    'Pr(>chisq)'  = pval.cure
  );
  rownames(coef.table.cure) <- colnames(design.matrix);

  coef.table.surv <- cbind(
      'LL.surv' = ll.surv,
    'LLR'         = llr.surv,
    'Pr(>chisq)'  = pval.surv
  );
  rownames(coef.table.surv) <- colnames(design.matrix);

  coef.table.alpha <- "NA";


#run.time = proc.time() - init.time

out <- list(
  coefficients = list(
    cure = coef.table.cure,
    surv = coef.table.surv,
    alpha = coef.table.alpha
 #   run.time
  )
);
class(out) <- c('mixcure.1d.nested.lrt', 'list');

return(out);

}

