#####################################################################
####      Function for LRT of mixcure model for              ########
####      FT-PL or ML via direct likelihood ratio comparison ########
#####################################################################

mixcure.penal.lrt.test <- function(formula, data, init, pl, iterlim = 200) {
require(splines)
require(survival)
require(abind)
  matrix.det<- function(matx) {

    #determinant for 2x2;
    TwobyTwo <- function(matfunc) {
      detm = matfunc[1,1]*matfunc[2,2]-matfunc[1,2]*matfunc[2,1]
      return(detm)
    }

    #determinant for 3x3;
    ThebyThe <- function(matfunc) {
      detm =  matfunc[1,1] * TwobyTwo(matfunc[-1,][,-1]) - matfunc[1,2] * TwobyTwo(matfunc[-1,][,-2]) + matfunc[1,3] * TwobyTwo(matfunc[-1,][,-3])
      return(detm)
    }

    #determinant for 4x4;
    FourbyFour <- function(matfunc){
      detm = matfunc[1,1] * ThebyThe(matfunc[-1,][,-1]) - matfunc[1,2] * ThebyThe(matfunc[-1,][,-2]) + matfunc[1,3] * ThebyThe(matfunc[-1,][,-3]) - matfunc[1,4] * ThebyThe(matfunc[-1,][,-4])
      return(detm)
    }

    #determinant for 5x5;
    FivebyFive <- function(matfunc){
      detm = matfunc[1,1] * FourbyFour(matfunc[-1,][,-1])- matfunc[1,2] * FourbyFour(matfunc[-1,][,-2]) + matfunc[1,3] * FourbyFour(matfunc[-1,][,-3])- matfunc[1,4] * FourbyFour(matfunc[-1,][,-4])+ matfunc[1,5] * FourbyFour(matfunc[-1,][,-5])
      return(detm)
    }

    #determinant for 6x6;
    SixbySix <- function(matfunc){
      detm = matfunc[1,1] * FivebyFive(matfunc[-1,][,-1]) - matfunc[1,2] * FivebyFive(matfunc[-1,][,-2]) + matfunc[1,3] * FivebyFive(matfunc[-1,][,-3]) - matfunc[1,4] * FivebyFive(matfunc[-1,][,-4]) + matfunc[1,5] * FivebyFive(matfunc[-1,][,-5]) - matfunc[1,6] * FivebyFive(matfunc[-1,][,-6])
      return(detm)
    }

    #determinant for 7x7;
    SevbySev <- function(matfunc){
      detm = matfunc[1,1] * SixbySix(matfunc[-1,][,-1]) - matfunc[1,2] * SixbySix(matfunc[-1,][,-2]) + matfunc[1,3] * SixbySix(matfunc[-1,][,-3]) - matfunc[1,4] * SixbySix(matfunc[-1,][,-4]) + matfunc[1,5] * SixbySix(matfunc[-1,][,-5]) - matfunc[1,6] * SixbySix(matfunc[-1,][,-6]) + matfunc[1,7] * SixbySix(matfunc[-1,][,-7])
      return(detm)
    }

    #determinant for 8x8;
    EigbyEig <- function(matfunc){
      detm = matfunc[1,1] * SevbySev(matfunc[-1,][,-1]) - matfunc[1,2] * SevbySev(matfunc[-1,][,-2]) + matfunc[1,3] * SevbySev(matfunc[-1,][,-3]) - matfunc[1,4] * SevbySev(matfunc[-1,][,-4]) + matfunc[1,5] * SevbySev(matfunc[-1,][,-5]) - matfunc[1,6] * SevbySev(matfunc[-1,][,-6]) + matfunc[1,7] * SevbySev(matfunc[-1,][,-7]) - matfunc[1,8] * SevbySev(matfunc[-1,][,-8])
      return(detm)
    }

    #determinant for 9x9;
    NinbyNin <- function(matfunc){
      detm = matfunc[1,1] * EigbyEig(matfunc[-1,][,-1]) - matfunc[1,2] * EigbyEig(matfunc[-1,][,-2]) + matfunc[1,3] * EigbyEig(matfunc[-1,][,-3]) - matfunc[1,4] * EigbyEig(matfunc[-1,][,-4]) + matfunc[1,5] * EigbyEig(matfunc[-1,][,-5]) - matfunc[1,6] * EigbyEig(matfunc[-1,][,-6]) + matfunc[1,7] * EigbyEig(matfunc[-1,][,-7]) - matfunc[1,8] * EigbyEig(matfunc[-1,][,-8]) + matfunc[1,9] * EigbyEig(matfunc[-1,][,-9])
      return(detm)
    }

    #determinant for 10x10;
    TenbyTen <- function(matfunc){
      detm = matfunc[1,1] * NinbyNin(matfunc[-1,][,-1]) - matfunc[1,2] * NinbyNin(matfunc[-1,][,-2]) + matfunc[1,3] * NinbyNin(matfunc[-1,][,-3]) - matfunc[1,4] * NinbyNin(matfunc[-1,][,-4]) + matfunc[1,5] * NinbyNin(matfunc[-1,][,-5]) - matfunc[1,6] * NinbyNin(matfunc[-1,][,-6]) + matfunc[1,7] * NinbyNin(matfunc[-1,][,-7]) - matfunc[1,8] * NinbyNin(matfunc[-1,][,-8]) + matfunc[1,9] * NinbyNin(matfunc[-1,][,-9]) - matfunc[1,10] * NinbyNin(matfunc[-1,][,-10])
      return(detm)
    }

    #determinant for 11x11;

    ElebyEle <- function(matfunc){
      detm = matfunc[1,1] * TenbyTen(matfunc[-1,][,-1]) - matfunc[1,2] * TenbyTen(matfunc[-1,][,-2]) + matfunc[1,3] * TenbyTen(matfunc[-1,][,-3]) - matfunc[1,4] * TenbyTen(matfunc[-1,][,-4]) + matfunc[1,5] * TenbyTen(matfunc[-1,][,-5]) - matfunc[1,6] * TenbyTen(matfunc[-1,][,-6]) + matfunc[1,7] * TenbyTen(matfunc[-1,][,-7]) - matfunc[1,8] * TenbyTen(matfunc[-1,][,-8]) + matfunc[1,9] * TenbyTen(matfunc[-1,][,-9]) - matfunc[1,10] * TenbyTen(matfunc[-1,][,-10]) + matfunc[1,11] * TenbyTen(matfunc[-1,][,-11])
      return(detm)
    }


    if (ncol(matx) == 1){
      determt = matx
    }

    if (ncol(matx) == 2){
      determt = TwobyTwo(matx)
    }

    if (ncol(matx) == 3){
      determt = ThebyThe(matx)
    }

    if (ncol(matx) == 4){
      determt = FourbyFour(matx)
    }

    if (ncol(matx) == 5){
      determt = FivebyFive(matx)
    }

    if (ncol(matx) == 6){
      determt = SixbySix(matx)
    }

    if (ncol(matx) == 7){
      determt = SevbySev(matx)
    }

    if (ncol(matx) == 8){
      determt = EigbyEig(matx)
    }

    if (ncol(matx) == 9){
      determt = NinbyNin(matx)
    }

    if (ncol(matx) == 10){
      determt = TenbyTen(matx)
    }

    if (ncol(matx) == 11){
      determt = ElebyEle(matx)
    }
    return(determt)
  }

  #########################################################################################
  mat.inv <- function(matx) {

    detm = matrix.det(matx)
    #2x2 matrix inverse;
    if (ncol(matx) == 2) {
      inv.matx = (1/detm) * matrix(c(matx[2,2],-matx[2,1],-matx[1,2],matx[1,1]), nrow = 2)
    }

    else {
      #For any n>2 dimension square matrix;
      adjug.matx <- matrix(rep(0, ncol(matx)^2), nrow = nrow(matx))
      for (i in 1:nrow(matx)) {
        for (j in 1:ncol(matx)) {
          adjug.matx[i,j] <- (-1)^(i+j)*matrix.det(matx[-i,][,-j])
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


#############################################################
#### loglik function of full model

  loglik.mixture <- function(p, survt, design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {

    ####  parameter and variable dependent parameters;
    #####
    theta = 1/(1+exp(-design.matrix%*%p[index.cure.var]))
    eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.var])
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap = -theta*(1-theta)*(1-eta)+(1-theta)^2*eta*(1-eta)
    pi = exp(eps)*eps*eta^2

    #calculate loglikelihood for the unpenalized;
    cure.par <- p[1 : ncol(design.matrix) ];
    surv.par <- p[ (ncol(design.matrix) + 1) : (2*length(cure.par)) ];
    p.gamma <- p[ 2*length(cure.par) + 1 ];  #use original shape parameter instead of exp();

    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) -
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );

    if (pl==T) {
      ####calculate inverse of info matrix by block matrix;

      n.elema = length(index.cure.var)^2
      a.sub1 <- matrix(rep(0,n.elema), nrow = length(index.cure.var))
      a.sub2 <- matrix(rep(0,n.elema), nrow = length(index.cure.var))

      for (i in c(index.cure.var)) {
        for (j in c(index.cure.var)) {
          a.sub1[i,j] <- sum((design.matrix[,i]*design.matrix[,j]*theta*(1-theta))[survt[, 2] == 1])
          a.sub2[i,j] <- sum((design.matrix[,i]*design.matrix[,j]*kap)[survt[, 2] == 0])
        }
      }
      info.a = a.sub1 + a.sub2

      design.xt <- cbind(design.matrix, log(survt[,1]))
      n.elemb <- length(index.cure.var)*(length(index.cure.var)+1)
      b.sub <- matrix(rep(0,n.elemb), nrow = length(index.surv.var))

      for (i in c(index.cure.var)) {
        for (j in c(index.cure.var,length(index.surv.var)+1)) {
          b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*eps*(1-delta)*delta)[survt[, 2] == 0])

        }
      }
      info.b = b.sub  #Upper right block of fisher.info;


      n.elemd <- (length(index.surv.var)+1)^2
      d.sub1 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
      d.sub2 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))

      for (i in c(index.cure.var,length(index.surv.var)+1)) {
        for (j in c(index.cure.var,length(index.surv.var)+1)) {
          d.sub1[i,j] <- sum((design.xt[,i]*design.xt[,j]*eps)[survt[, 2] == 1])
          d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*(delta*(1-delta))))[survt[, 2] == 0])

        }
      }
      info.d = d.sub1 + d.sub2 +
        matrix(c(rep(0, (n.elemd-1)),sum(survt[, 2] == 1)/(p[index.gamma]^2)),nrow = (length(index.surv.var)+1))


      info.d.inv = mat.inv(info.d)

       info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)

      #determinant of hessian matrix;
      det.info = matrix.det(info.set0)*matrix.det(info.d)

       loglik = loglikelihood - 0.5*log(det.info)
    } else if (pl == FALSE) {
      loglik = loglikelihood
    }

        return(loglik)

  }

  ######END of loglik.mixture####################################


  # Parameter estimation under Ha (non-restricted likelihood)
  # maximize penalized or unpenalized loglikelihood by nlm;
  maximizer0 <- nlm(
    f = loglik.mixture, p = init, survt=survt, design.matrix=design.matrix,
    pl = pl,
    iterlim = iterlim, hessian=F);

  loglik <- -maximizer0$minimum  #in loglik function loglik was calculated as minus of actual loglik value


#############################################################


  ##loglik function for testing parameters of cure or surv part;
  #design.matrix1-surv part, design.matrix0-cure part;

  loglik.mixture.part <- function(p, survt, design.matrix1, design.matrix0, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {

    design.mtx.comb = cbind(design.matrix0,design.matrix1)


    #parameter and variable dep parameters;
    if (k > length(index.cure.v)) {
    theta = 1/(1+exp(-design.matrix[,index.cure.var]%*%as.matrix(p[index.cure.var])))
    eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[-c(index.cure.var,index.gamma-1)]))
    } else {
    theta = 1/(1+exp(-design.matrix[,index.cure.var]%*%as.matrix(p[-c(index.surv.var-1,index.gamma-1)])))
    eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.var-1]))
    }

    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap = -theta*(1-theta)*(1-eta)+(1-theta)^2*eta*(1-eta)
    pi = exp(eps)*eps*eta^2

    ####################################################################################################
    # Note: below constructs fisher info matrix; steps are divide into 4 blocks, 2 square blocks (A&D) #
    # on upper left and lower right, 2 identical transposed blocks (B) on upper right and lower left;  #
    # the identical B blocks are not identical in reduced models unless it's a global LRT, needs to be C#
    ####################################################################################################

    #calculate loglikelihood for the unpenalized;
    p.gamma <- p[index.gamma-1];  #use original shape parameter instead of exp();

    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer;
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             + log(eps)-eps )[survt[, 2] == 1] ) -
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );

    if (pl == F) {loglik.part = loglikelihood} else {


    max.len = max(length(index.cure.var),length(index.surv.var))
    n.elema = max.len^2
    a.sub1 <- matrix(rep(0,n.elema), nrow = max.len)
    a.sub2 <- matrix(rep(0,n.elema), nrow = max.len)

    for (i in c(index.cure.var)) {
      for (j in c(index.cure.var)) {
        a.sub1[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*theta*(1-theta))[survt[, 2] == 1])
        a.sub2[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*kap)[survt[, 2] == 0])
      }
    }
    info.a = (a.sub1 + a.sub2)[index.cure.var,index.cure.var]

    ##info matrix block B
    design.xt0 <- cbind(design.matrix0, log(survt[,1]))
    n.elemb <- max.len*(max.len+1)
    b.sub <- matrix(rep(0,n.elemb), nrow = max.len)

    for (i in c(index.cure.var)) {
      for (j in c((index.surv.var-max.len), max.len+1)) {
        b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*eps*(1-delta)*delta)[survt[, 2] == 0])  #equivalent to expression below
      }
    }
    info.b = b.sub[index.cure.var,c(index.surv.var-max.len,index.gamma-max.len)]

        design.xt1 <- cbind(design.matrix1, log(survt[,1]))

    n.elemd <- (max.len+1)^2
    d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
    d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))

    for (i in c(index.surv.var-max.len, max.len +1)) {
      for (j in c(index.surv.var-max.len, max.len +1)) {
        d.sub1[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*eps)[survt[, 2] == 1])
        d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*(delta*(1-delta))))[survt[, 2] == 0])

      }
    }
    d.sub = d.sub1 + d.sub2 +
      matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-1]^2)),
             nrow = (max.len + 1))

    info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]

    info.d.inv = mat.inv(info.d)

       # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
    info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)

    det.info = matrix.det(info.set0)*matrix.det(info.d)

    loglik.part = loglikelihood - 0.5*log(det.info)

      }

    return(loglik.part)
  }


  #################################################################

  #### parameter estimation under H0 for individual parameter
  #### loglikelihood ratio test statistics for each of cure part variables;

  dim.v <- ncol(design.matrix)
  ll.cure <- rep(0,dim.v)
  llr.cure <- rep(0,dim.v)
  pval.cure <- rep(0,dim.v)

  for (k in index.cure.v[-1]) {
      maximizer <- nlm(
      f = loglik.mixture.part,
      p = init[-k],
      survt = survt, design.matrix0 = design.matrix,
      design.matrix1=design.matrix,
      index.cure.var=index.cure.v[-k],
      pl=pl,
      iterlim = iterlim, hessian=F
    );

    loglik.part = -maximizer$minimum;
    dif.ll = -2*(loglik.part-loglik);  #loglik is ll under non-restricted model;
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.cure[k]<- loglik.part
    llr.cure[k]<- dif.ll
    pval.cure[k]<- pval

  }


  ### loglikelihood calculation for each surv part variable;

  ll.surv <- rep(0,ncol(design.matrix))
  llr.surv <- rep(0,ncol(design.matrix))
  pval.surv <- rep(0,ncol(design.matrix))


  for (k in index.surv.v[-1]) {
    is=k-length(index.cure.v)
    maximizer <- nlm(
      f = loglik.mixture.part, p =  init[-k],
      survt = survt, design.matrix1 = design.matrix,
      design.matrix0=design.matrix,
      index.surv.var=index.surv.v[-is],
      pl=pl,
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


#######################################
## Output tables from FT-PLE or ML; ##
#######################################

out <- list(
  coefficients = list(
    cure = coef.table.cure,
    surv = coef.table.surv,
    alpha = coef.table.alpha
 #  run.time
  )
);
class(out) <- c('mixcure', 'list');

return(out);

}


