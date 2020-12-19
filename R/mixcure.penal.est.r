################################
### Function for mixcure model##
#### penalized loglikelihoods ##
##########################################################
#### Last modified Dec31 2018 for 3 or more x variables ##
##########################################################

mixcure.penal.est <- function(formula, data, init, pl, iterlim = 200) { 
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
     detm = matfunc[1,1] * NinbyNin(matfunc[-1,][,-1]) - matfunc[1,2] * NinbyNin(matfunc[-1,][,-2]) + 
            matfunc[1,3] * NinbyNin(matfunc[-1,][,-3]) - matfunc[1,4] * NinbyNin(matfunc[-1,][,-4]) + 
            matfunc[1,5] * NinbyNin(matfunc[-1,][,-5]) - matfunc[1,6] * NinbyNin(matfunc[-1,][,-6]) + 
            matfunc[1,7] * NinbyNin(matfunc[-1,][,-7]) - matfunc[1,8] * NinbyNin(matfunc[-1,][,-8]) + 
            matfunc[1,9] * NinbyNin(matfunc[-1,][,-9]) - matfunc[1,10] * NinbyNin(matfunc[-1,][,-10])
     return(detm)
    }
    #determinant for 11x11;
    ElebyEle <- function(matfunc){
      detm = matfunc[1,1] * TenbyTen(matfunc[-1,][,-1]) - matfunc[1,2] * TenbyTen(matfunc[-1,][,-2]) + 
             matfunc[1,3] * TenbyTen(matfunc[-1,][,-3]) - matfunc[1,4] * TenbyTen(matfunc[-1,][,-4]) + 
             matfunc[1,5] * TenbyTen(matfunc[-1,][,-5]) - matfunc[1,6] * TenbyTen(matfunc[-1,][,-6]) + 
             matfunc[1,7] * TenbyTen(matfunc[-1,][,-7]) - matfunc[1,8] * TenbyTen(matfunc[-1,][,-8]) + 
             matfunc[1,9] * TenbyTen(matfunc[-1,][,-9]) - matfunc[1,10] * TenbyTen(matfunc[-1,][,-10]) + 
             matfunc[1,11] * TenbyTen(matfunc[-1,][,-11])
      return(detm)
    }
    
    TwlbyTwl <- function(matfunc){
      detm = matfunc[1,1] * ElebyEle(matfunc[-1,][,-1]) - matfunc[1,2] * ElebyEle(matfunc[-1,][,-2]) + 
        matfunc[1,3] * ElebyEle(matfunc[-1,][,-3]) - matfunc[1,4] * ElebyEle(matfunc[-1,][,-4]) + 
        matfunc[1,5] * ElebyEle(matfunc[-1,][,-5]) - matfunc[1,6] * ElebyEle(matfunc[-1,][,-6]) + 
        matfunc[1,7] * ElebyEle(matfunc[-1,][,-7]) - matfunc[1,8] * ElebyEle(matfunc[-1,][,-8]) + 
        matfunc[1,9] * ElebyEle(matfunc[-1,][,-9]) - matfunc[1,10] * ElebyEle(matfunc[-1,][,-10]) + 
        matfunc[1,11] * ElebyEle(matfunc[-1,][,-11]) - matfunc[1,12] * ElebyEle(matfunc[-1,][,-12])
      return(detm)
    }
    
    TitbyTit <- function(matfunc){
      detm = matfunc[1,1] * TwlbyTwl(matfunc[-1,][,-1]) - matfunc[1,2] * TwlbyTwl(matfunc[-1,][,-2]) + 
        matfunc[1,3] * TwlbyTwl(matfunc[-1,][,-3]) - matfunc[1,4] * TwlbyTwl(matfunc[-1,][,-4]) + 
        matfunc[1,5] * TwlbyTwl(matfunc[-1,][,-5]) - matfunc[1,6] * TwlbyTwl(matfunc[-1,][,-6]) + 
        matfunc[1,7] * TwlbyTwl(matfunc[-1,][,-7]) - matfunc[1,8] * TwlbyTwl(matfunc[-1,][,-8]) + 
        matfunc[1,9] * TwlbyTwl(matfunc[-1,][,-9]) - matfunc[1,10] * TwlbyTwl(matfunc[-1,][,-10]) + 
        matfunc[1,11] * TwlbyTwl(matfunc[-1,][,-11]) - matfunc[1,12] * TwlbyTwl(matfunc[-1,][,-12]) +
        matfunc[1,13] * TwlbyTwl(matfunc[-1,][,-13])
      return(detm)
    }
    
    FotbyFot <- function(matfunc){
      detm = matfunc[1,1] * TitbyTit(matfunc[-1,][,-1]) - matfunc[1,2] * TitbyTit(matfunc[-1,][,-2]) + 
        matfunc[1,3] * TitbyTit(matfunc[-1,][,-3]) - matfunc[1,4] * TitbyTit(matfunc[-1,][,-4]) + 
        matfunc[1,5] * TitbyTit(matfunc[-1,][,-5]) - matfunc[1,6] * TitbyTit(matfunc[-1,][,-6]) + 
        matfunc[1,7] * TitbyTit(matfunc[-1,][,-7]) - matfunc[1,8] * TitbyTit(matfunc[-1,][,-8]) + 
        matfunc[1,9] * TitbyTit(matfunc[-1,][,-9]) - matfunc[1,10] * TitbyTit(matfunc[-1,][,-10]) + 
        matfunc[1,11] * TitbyTit(matfunc[-1,][,-11]) - matfunc[1,12] * TitbyTit(matfunc[-1,][,-12]) +
        matfunc[1,13] * TitbyTit(matfunc[-1,][,-13]) - matfunc[1,14] * TitbyTit(matfunc[-1,][,-14])
      return(detm)
    }
    
    FifbyFif <- function(matfunc){
      detm = matfunc[1,1] * FotbyFot(matfunc[-1,][,-1]) - matfunc[1,2] * FotbyFot(matfunc[-1,][,-2]) + 
        matfunc[1,3] * FotbyFot(matfunc[-1,][,-3]) - matfunc[1,4] * FotbyFot(matfunc[-1,][,-4]) + 
        matfunc[1,5] * FotbyFot(matfunc[-1,][,-5]) - matfunc[1,6] * FotbyFot(matfunc[-1,][,-6]) + 
        matfunc[1,7] * FotbyFot(matfunc[-1,][,-7]) - matfunc[1,8] * FotbyFot(matfunc[-1,][,-8]) + 
        matfunc[1,9] * FotbyFot(matfunc[-1,][,-9]) - matfunc[1,10] * FotbyFot(matfunc[-1,][,-10]) + 
        matfunc[1,11] * FotbyFot(matfunc[-1,][,-11]) - matfunc[1,12] * FotbyFot(matfunc[-1,][,-12]) +
        matfunc[1,13] * FotbyFot(matfunc[-1,][,-13]) - matfunc[1,14] * FotbyFot(matfunc[-1,][,-14]) +
        matfunc[1,15] * FotbyFot(matfunc[-1,][,-15])
        return(detm)
    }
    
    SitbySit <- function(matfunc){
      detm = matfunc[1,1] * FifbyFif(matfunc[-1,][,-1]) - matfunc[1,2] * FifbyFif(matfunc[-1,][,-2]) + 
        matfunc[1,3] * FifbyFif(matfunc[-1,][,-3]) - matfunc[1,4] * FifbyFif(matfunc[-1,][,-4]) + 
        matfunc[1,5] * FifbyFif(matfunc[-1,][,-5]) - matfunc[1,6] * FifbyFif(matfunc[-1,][,-6]) + 
        matfunc[1,7] * FifbyFif(matfunc[-1,][,-7]) - matfunc[1,8] * FifbyFif(matfunc[-1,][,-8]) + 
        matfunc[1,9] * FifbyFif(matfunc[-1,][,-9]) - matfunc[1,10] * FifbyFif(matfunc[-1,][,-10]) + 
        matfunc[1,11] * FifbyFif(matfunc[-1,][,-11]) - matfunc[1,12] * FifbyFif(matfunc[-1,][,-12]) +
        matfunc[1,13] * FifbyFif(matfunc[-1,][,-13]) - matfunc[1,14] * FifbyFif(matfunc[-1,][,-14]) +
        matfunc[1,15] * FifbyFif(matfunc[-1,][,-15]) - matfunc[1,16] * FifbyFif(matfunc[-1,][,-16])
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
    
    if (ncol(matx) == 12){
      determt = TwlbyTwl(matx)
    }
    
    if (ncol(matx) == 13){
      determt = TitbyTit(matx)
    }
    
    if (ncol(matx) == 14){
      determt = FotbyFot(matx)
    }

    if (ncol(matx) == 15){
      determt = FifbyFif(matx)
    }
    if (ncol(matx) == 16){
      determt = SitbySit(matx)
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
  
  #samp.s <- nrow(design.matrix)
  
  
  ####################################################
  ## nonlinear minimization algoritm to solve       ##
  ## penalized mixture cure loglikelihood functions ##
  ####################################################      
  
  loglik.mixture <- function(p, survt, design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl, PLCI=F) {
    
   ####  parameter and variable dep parameters;
   #####
    theta = 1/(1+exp(-design.matrix%*%p[index.cure.var]))
    eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.var])
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap = -theta*(1-theta)*(1-eta)+(1-theta)^2*eta*(1-eta)
    #kap= (1-eta)*(1-theta)*(theta + eta)
    pi = exp(eps)*eps*eta^2
    # lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    # phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi

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
        b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])
        #b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])
        
              }
    }
    info.b = b.sub  #Upper right block of fisher.info;
    

    n.elemd <- (length(index.surv.var)+1)^2
    d.sub1 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
    d.sub2 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
    
    for (i in c(index.cure.var,length(index.surv.var)+1)) {
      for (j in c(index.cure.var,length(index.surv.var)+1)) {
        d.sub1[i,j] <- sum((design.xt[,i]*design.xt[,j]*eps)[survt[, 2] == 1])
        #d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
        #d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta^2))[survt[, 2] == 0])
        d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*(delta*(1-delta))))[survt[, 2] == 0])
        
        }
    }
    info.d = d.sub1 + d.sub2 + 
           matrix(c(rep(0, (n.elemd-1)),sum(survt[, 2] == 1)/(p[index.gamma]^2)),nrow = (length(index.surv.var)+1))
    
    
    info.d.inv = mat.inv(info.d)

#    fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
    #hessian.mat = -fisher.info
    
    # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
    info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)

    #determinant of hessian matrix;
    det.info = matrix.det(info.set0)*matrix.det(info.d)
 #   det.info = matrix.det(fisher.info)
      
    if (PLCI==T)  
    {loglik = loglikelihood - 0.5*log(abs(det.info))} else {
      
      loglik = loglikelihood - 0.5*log(det.info)
    }
 }
    
   else if (pl == FALSE)
    {
      loglik = loglikelihood
    }

    #loglik = loglikelihood
    return(loglik)
    
  }
  
  ######END of loglik.mixture####################################
  
  
  # Parameter estimation under Ha (non-restricted likelihood)
  # maximize penalized or unpenalized loglikelihood by nlm; 
  maximizer0 <- nlm(
    f = loglik.mixture, p = init, survt=survt, design.matrix=design.matrix, 
    pl = pl, 
    iterlim = iterlim, hessian=TRUE);

  hessmat <- maximizer0$hessian
  if (det(hessmat) < 1e-05) 
    diag(hessmat) <- diag(hessmat) + 1e-06
  
var.mat <- solve(hessmat)

alpha.hat <- maximizer0$estimate[index.gamma];

loglik <- -maximizer0$minimum  #in loglik function loglik was calculated as minus of actual loglik value

######### Wald statistic inference ##########
  # confidence intervals for estimated coefficients
  z.score <- maximizer0$estimate / sqrt(diag(var.mat));
  
coef.table.cure <- cbind(
  'coef'        = maximizer0$estimate[index.cure.v],
  'exp(coef)'   = exp(maximizer0$estimate[index.cure.v]),
  'se(coef)'    = sqrt(diag(var.mat)[index.cure.v]),
  'z'           = z.score[index.cure.v],
  'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.cure.v]))),
  'LCI.95%' = maximizer0$estimate[index.cure.v] - 1.96 * sqrt(diag(var.mat)[index.cure.v]),
  'UCI.95%' = maximizer0$estimate[index.cure.v] + 1.96 * sqrt(diag(var.mat)[index.cure.v])
);
rownames(coef.table.cure) <- colnames(design.matrix);

coef.table.surv <- cbind(
  'coef'        = maximizer0$estimate[index.surv.v],
  'exp(coef)'   = exp(maximizer0$estimate[index.surv.v]),
  'se(coef)'    = sqrt(diag(var.mat)[index.surv.v]),
  'z'           = z.score[index.surv.v],
  'Pr(>|z|)'    = 2 * (1 - pnorm(abs(z.score[index.surv.v]))),
  'LCI.95%' = maximizer0$estimate[index.surv.v] - 1.96 * sqrt(diag(var.mat)[index.surv.v]),
  'UCI.95%' = maximizer0$estimate[index.surv.v] + 1.96 * sqrt(diag(var.mat)[index.surv.v])
);
rownames(coef.table.surv) <- colnames(design.matrix);

coef.table.alpha <- cbind(
  'coef'     = alpha.hat,
  'se(coef)' = sqrt(diag(var.mat)[index.gamma]),
  'z'        = z.score[index.gamma],
  'Pr(>|z|)' = 2 * (1 - pnorm(abs(z.score[index.gamma]))),
  'LCI.95%'  = maximizer0$estimate[index.gamma] - 1.96 * sqrt(diag(var.mat)[index.gamma]),
  'UCI.95%'  = maximizer0$estimate[index.gamma] + 1.96 * sqrt(diag(var.mat)[index.gamma]),
  'loglik' = -maximizer0$minimum 
);
rownames(coef.table.alpha) <- 'alpha';

#######################################
## Output tables from either method; ##
#######################################

colnames(var.mat) <- c(
  paste('cure.', colnames(design.matrix)), 
  paste('surv.', colnames(design.matrix)),
  'alpha'
);
rownames(var.mat) <- colnames(var.mat);

out <- list(
  coefficients = list(
    cure = coef.table.cure, 
    surv = coef.table.surv, 
    alpha = coef.table.alpha
 #   run.time
  ),
  cov = var.mat
);
class(out) <- c('mixcure', 'list');

return(out);

}


#### print.mixcure #############################################################
# DESCRIPTION
#   To print a mixcure object.
# INPUT
#   object : a mixcure object, which is an outcome of function mixcure.
#   digits : number of digits for printing, passed to print.default.
#   ...    : other parameters passed to print.default.
# OUTPUT
#   NULL.   

print.mixcure <- function(object, digits = 3, ...) {
  sep.line.cure   <- paste(c(rep('-', 37),  ' CURE ' , rep('-', 37)), collapse = '');
  sep.line.surv   <- paste(c(rep('-', 37),  ' SURVIVAL ' , rep('-', 37)), collapse = '');
  sep.line.alpha <- paste(c(rep('-', 36), ' ALPHA ', rep('-', 36)), collapse = '');
  
  message(sep.line.cure);
  print.default(object$coefficients$cure,   digits = digits, ...);
  
  message(sep.line.surv);
  print.default(object$coefficients$surv,   digits = digits, ...);
  
  message(sep.line.alpha);
  print.default(object$coefficients$alpha, digits = digits, ...);
  
  return(NULL);
};


#### coef.mixcure ##############################################################
coef.mixcure <- function(object) {
  coefs <- c(
    object$coefficients$cure[, 'coef'], 
    object$coefficients$surv[, 'coef'], 
    object$coefficients$alpha[, 'coef']
  );
  names(coefs) <- c( paste('cure.', rownames(object$coefficients$cure)), 
                     paste('surv.', rownames(object$coefficients$surv)),
                     rownames(object$coefficients$alpha) );
  
  return(coefs);
}


vcov.mixcure <- function(object) {
  return(object$cov);
}


