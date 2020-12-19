################################
### Function for mixcure model##
#### penalized loglikelihoods ##
###################################################
#### Last modified Dec31 2018 for one x variable ##
###################################################

########### NOTE on Oct22:
########## Error checking for cure uppper points iter2;
########## maximizer$temp1 is not consistent in and out of iter2 loop;

mixcure.penal.profile.CI.test <- function(formula, data, init, pl , method = "LR", apct = 0.05, LRT.pval = F, iterlim = 200) { 
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
  

  design.matrix <- model.frame(formula, data = data, na.action = na.omit);
  survt <- design.matrix[,1];
  
  design.matrix <- model.matrix(formula, data = design.matrix);
  
  # index ranges of coefficients of glm and cox models
  index.cure.v <- 1 : ncol(design.matrix); 
  index.surv.v <- (ncol(design.matrix) + 1) : (2*length(index.cure.v))
  # index of alpha,the shape parameter
  index.gamma <- 2*length(index.cure.v)+1;
  
  samp.s <- nrow(design.matrix)
  
  
  ####################################################
  ## nonlinear minimization algoritm to solve       ##
  ## penalized mixture cure loglikelihood functions ##
  ####################################################      
  
  # > maximizer0$minimum
  # [1] 600.1339
  # > maximizer0$estimate
  # [1]  1.71643676 -0.17877823  0.65814299  0.61905752  0.09867225 -0.77065440 -7.63529053  0.85847947 -0.20121321  1.00244978  0.71073735 -0.22989348
  # [13]  1.81956226
  # op.est = c(1.71643676, 0,  0.65814299,  0.61905752,  0.09867225, -0.77065440, -7.63529053,  0.85847947, -0.20121321,  1.00244978,  0.71073735, -0.22989348, 1.81956226)
  # > loglik.part
  # [1] 600.3221
  # 
  # op.est = c(1.71643676, -0.17877823,  0,  0.61905752,  0.09867225, -0.77065440, -7.63529053,  0.85847947, -0.20121321,  1.00244978,  0.71073735, -0.22989348, 1.81956226)
  # > loglik.part
  # [1] 604.1643
  # 
  # op.est = c(1.71643676, -0.17877823,  0.65814299,  0,  0.09867225, -0.77065440, -7.63529053,  0.85847947, -0.20121321,  1.00244978,  0.71073735, -0.22989348, 1.81956226)
  # > loglik.part
  # [1] 603.4383
  # 
  # p = op.est;
  # design.matrix1 = design.matrix; design.matrix0 = design.matrix; index.cure.var=index.cure.v;index.surv.var=index.surv.v;
  
  loglik.mixture <- function(p, survt, design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {
    
   ####  parameter and variable dep parameters;
   #####
    theta = 1/(1+exp(-design.matrix%*%p[index.cure.var]))
    eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.var])
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    #kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta)
    kap= (1-eta)*(1-theta)*(theta + eta)
    pi = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi

    

 
  ####calculate inverse of info matrix by block matrix;
    ################
    # a.sub1xsq = sum((design.matrix[,-1]^2*theta*(1-theta))[survt[, 2] == 1])
    # a.sub2xsq = sum((design.matrix[,-1]^2*kap)[survt[, 2] == 0])
    # a.sub1x = sum((design.matrix[,-1]*theta*(1-theta))[survt[, 2] == 1])
    # a.sub2x = sum((design.matrix[,-1]*kap)[survt[, 2] == 0])
    # a.sub1 = sum((theta*(1-theta))[survt[, 2] == 1])
    # a.sub2 = sum(kap[survt[, 2] == 0])
    # info.a = rbind(cbind((a.sub1+a.sub2),(a.sub1x+a.sub2x)),
    #                cbind((a.sub1x+a.sub2x),(a.sub1xsq+a.sub2xsq)))

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
    
    
    # det.a = (a.sub1xsq+a.sub2xsq)*(a.sub1+a.sub2)-(a.sub1x+a.sub2x)^2
    # info.a.inv <- (1/det.a)*matrix(c((a.sub1xsq+a.sub2xsq),((a.sub1x+a.sub2x)),
    #                                  ((a.sub1x+a.sub2x)),(a.sub1+a.sub2)),
    #                                nrow=length(index.cure.var))
    
    # b.subxsq = -sum((design.matrix[,-1]^2*theta*(1-theta)*pi)[survt[, 2] == 0])
    # b.subx = -sum((design.matrix[,-1]*theta*(1-theta)*pi)[survt[, 2] == 0])
    # b.sub = -sum((theta*(1-theta)*pi)[survt[, 2] == 0])
    # b.subxt = -sum((design.matrix[,-1]*log(survt[,1])*theta*(1-theta)*pi)[survt[, 2] == 0])
    # b.subt = -sum((log(survt[,1])*theta*(1-theta)*pi)[survt[, 2] == 0])
    # info.b <- matrix(c(b.sub,b.subx,b.subx,b.subxsq,b.subt,b.subxt),nrow=length(index.cure.var))

    
    design.xt <- cbind(design.matrix, log(survt[,1]))
    n.elemb <- length(index.cure.var)*(length(index.cure.var)+1)
    b.sub <- matrix(rep(0,n.elemb), nrow = length(index.surv.var))
    
    for (i in c(index.cure.var)) {
      for (j in c(index.cure.var,length(index.surv.var)+1)) {
        b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])
      }
    }
    info.b = b.sub  #Upper right block of fisher.info;
    
    # d.sub1 = sum(eps[survt[, 2] == 1])
    # d.sub2 = sum((eps*delta-eps^2*delta+eps^2*delta^2)[survt[, 2] == 0])
    # d.sub1x = sum((design.matrix[,-1]*eps)[survt[, 2] == 1])
    # d.sub2x = sum((design.matrix[,-1]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
    # d.sub1xsq = sum((design.matrix[,-1]^2*eps)[survt[, 2] == 1])
    # d.sub2xsq = sum((design.matrix[,-1]^2*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
    # d.sub1t = sum((log(survt[,1])*eps)[survt[, 2] == 1])
    # d.sub2t = sum((log(survt[,1])*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
    # d.sub1xt = sum((design.matrix[,-1]*log(survt[,1])*eps)[survt[, 2] == 1])
    # d.sub2xt = sum((design.matrix[,-1]*log(survt[,1])*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
    # d.sub1tsq = sum((log(survt[,1])^2*eps)[survt[, 2] == 1])
    # d.sub2tsq = sum((log(survt[,1])^2*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
    # 
    # info.d <- matrix(c((d.sub1+d.sub2),(d.sub1x+d.sub2x),(d.sub1t+d.sub2t),
    #                    (d.sub1x+d.sub2x),(d.sub1xsq+d.sub2xsq),(d.sub1xt+d.sub2xt),
    #                    (d.sub1t+d.sub2t),(d.sub1xt+d.sub2xt),
    #                    (sum(survt[, 2] == 1)/(p[index.gamma]^2)+d.sub1tsq+d.sub2tsq)), nrow = 3)
    
    
    n.elemd <- (length(index.surv.var)+1)^2
    d.sub1 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
    d.sub2 <- matrix(rep(0,n.elemd), nrow = (length(index.surv.var)+1))
     for (i in c(index.cure.var,length(index.surv.var)+1)) {
      for (j in c(index.cure.var,length(index.surv.var)+1)) {
        d.sub1[i,j] <- sum((design.xt[,i]*design.xt[,j]*eps)[survt[, 2] == 1])
        d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
       # d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta^2))[survt[, 2] == 0])
        
        }
    }
    info.d = d.sub1 + d.sub2 + 
           matrix(c(rep(0, (n.elemd-1)),sum(survt[, 2] == 1)/(p[index.gamma]^2)),nrow = (length(index.surv.var)+1))
   
    
    
    info.d.inv = mat.inv(info.d)

    fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
    #hessian.mat = -fisher.info
    
    # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
    info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)

    #determinant of hessian matrix;
    det.info = matrix.det(info.set0)*matrix.det(info.d)
    #det.info = matrix.det(fisher.info)
      
    #calculate loglikelihood for the unpenalized;
    cure.par <- p[1 : ncol(design.matrix) ];
    surv.par <- p[ (ncol(design.matrix) + 1) : (2*length(cure.par)) ];
    p.gamma <- p[ 2*length(cure.par) + 1 ];  #use original shape parameter instead of exp();
    
    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer; 
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                            +log(eps)-eps )[survt[, 2] == 1] ) - 
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
    
    if (pl == FALSE)
    {
      loglik = loglikelihood
    }
    else if (pl == TRUE)
    {
          loglik = loglikelihood - 0.5*log(det.info)
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




if (method == "Wald") {

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

#############################################################
#############################################################
#############################################################

} else if (method == "LR") {
  #create CI for profile likelihood, this option only outputs estimates and PL under specified model;

  
  ## loglik function for testing parameters of cure or surv part ##
  loglik.mixture.part <- function(p, survt, design.matrix1, design.matrix0, 
                                  index.cure.var=index.cure.v,
                                  index.surv.var=index.surv.v, pl) {  #design.matrix1-surv, design.matrix0-cure
   
    design.mtx.comb = cbind(design.matrix0,design.matrix1)

    #parameter and variable dep parameters;
    theta = 1/(1+exp(-design.mtx.comb[,index.cure.var]%*%as.matrix(p[index.cure.var])))
    eps = survt[,1]^(p[index.gamma])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.var]))
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    #kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta)
    kap= (1-eta)*(1-theta)*(theta + eta)
    pi = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
    
    ####################################################################################################
    # Note: below constructs fisher info matrix; steps are divide into 4 blocks, 2 square blocks (A&D) #
    # on upper left and lower right, 2 identical transposed blocks (B) on upper right and lower left;  #
    # the idential B blocks are not identical in reduced models unless it's a global LRT, needs to be C#
    ####################################################################################################
    
    
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
    
    ####### For error check of info.a############
    # for (i in c(1:max.len)) {
    #   for (j in c(1:max.len)) {
    #     a.sub1[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*theta*(1-theta))[survt[, 2] == 1])
    #     a.sub2[i,j] <- sum((as.matrix(design.matrix0)[,i]*as.matrix(design.matrix0)[,j]*kap)[survt[, 2] == 0])
    #   }
    # }
    # info.a = (a.sub1 + a.sub2)[index.cure.var,index.cure.var]
    
    
       ##info matrix block B 
    design.xt0 <- cbind(design.matrix0, log(survt[,1]))
    n.elemb <- max.len*(max.len+1)
    b.sub <- matrix(rep(0,n.elemb), nrow = max.len)
    
    for (i in c(index.cure.var)) {
      for (j in c(1:length(index.surv.var), max.len+1)) {
        b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])
       # b.sub[i,j] <- -sum((design.matrix1[,i]*design.xt0[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])
        
        }
    }
    info.b = b.sub[index.cure.var,c(index.surv.var-max.len,index.gamma-max.len)]
    
    ###### For error checking of info.b#########
    # for (i in c(1:max.len)) {
    #   for (j in c(1:max.len, max.len + 1)) {
    #     b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])
    #   }
    # }
    # info.b = b.sub[index.cure.var,c(index.surv.var-max.len,index.gamma-max.len)]
    
  # ##info matrix block C  
     design.xt1 <- cbind(design.matrix1, log(survt[,1]))
  #   n.elemc <- max.len*(max.len+1)
  #   c.sub <- matrix(rep(0,n.elemc), ncol = max.len)
  #   
  #   for (i in c(index.cure.var)) {
  #     for (j in c(1:length(index.surv.var), length(index.surv.var)+1)) {
  #       c.sub[j,i] <- -sum((design.matrix1[,i]*design.xt1[,j]*delta*(1-delta)*eps)[survt[, 2] == 0])
  #     }
  #   }
  #   info.c = c.sub[c(index.surv.var-max.len,index.gamma-max.len),index.cure.var]
    

     n.elemd <- (max.len+1)^2
    d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
    d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
    
    for (i in c(index.surv.var-max.len, max.len +1)) {
      for (j in c(index.surv.var-max.len, max.len +1)) {
        d.sub1[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*eps)[survt[, 2] == 1])
        d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
       # d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta^2))[survt[, 2] == 0])
      }
    }
    d.sub = d.sub1 + d.sub2 + 
      matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma]^2)), 
             nrow = (max.len + 1))
    
    info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]
    
    
    info.d.inv = mat.inv(info.d)
    
    fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
    #hessian.mat = -fisher.info
    
    # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
    info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)
    
    #determinant of hessian matrix;
    det.info = matrix.det(info.set0)*matrix.det(info.d)
    #det.info = matrix.det(fisher.info)
    
    #calculate loglikelihood for the unpenalized;
    cure.par <- p[index.cure.var];
    surv.par <- p[index.surv.var];
    p.gamma <- p[index.gamma];  #use original shape parameter instead of exp();
    
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
  
 
  
  ## loglik function for constructing profile likelihood of cure or surv part ##
  loglik.mixture.profile <- function(p, survt, k=k, design.matrix1=design.matrix, design.matrix0=design.matrix, param.est, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl) {  
    

    design.mtx.comb = cbind(design.matrix0,design.matrix1)
    ik = k-length(index.cure.var);
    
    #parameter and variable dep parameters;
    if (k > length(index.cure.v)) {
      theta = 1/(1+exp(-design.matrix[,index.cure.var]%*%as.matrix(p[index.cure.var])))
    } else {
      theta = 1/(1+exp(-design.matrix[,index.cure.var[-k]]%*%as.matrix(p[-c(index.surv.var-1,index.gamma-1)])-design.mtx.comb[,k]*param.est))
    }
    if (k > length(index.cure.v)) {
      eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var[-ik]]%*%as.matrix(p[-c(index.cure.var,index.gamma-1)])+design.mtx.comb[,k]*param.est)
    } else {
      eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.var-1]))
    }
    
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    #kap = theta*(1-theta)*(1-eta)-(1-theta)^2*eta*(1-eta)
    kap= (1-eta)*(1-theta)*(theta + eta)
    pi = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
    
    ####################################################################################################
    # Note: below constructs fisher info matrix; steps are divide into 4 blocks, 2 square blocks (A&D) #
    # on upper left and lower right, 2 identical transposed blocks (B) on upper right and lower left;  #
    # the idential B blocks are not identical in reduced models unless it's a global LRT, needs to be C#
    ####################################################################################################
    
    #calculate loglikelihood for the unpenalized;
    p.gamma <- p[index.gamma-1];  #use original shape parameter instead of exp();
    
    # loglikelihood is defined as the negative of the actual loglikelihood for feeding nlm() minimizer; 
    loglikelihood <- -sum( ( log(1-theta) + log(p.gamma)-log(survt[,1])
                             +log(eps)-eps )[survt[, 2] == 1] ) - 
      sum( (log(theta + (1-theta)*exp(-eps)))[survt[, 2] == 0] );
    
    if (pl==T) {
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
   # if (k <= length(index.cure.var)) {info.a = (a.sub1 + a.sub2)[index.cure.var[-k],index.cure.var[-k]]} else
      info.a = (a.sub1 + a.sub2)[index.cure.var,index.cure.var]
  
    ##info matrix block B 
    design.xt0 <- cbind(design.matrix0, log(survt[,1]))
    n.elemb <- max.len*(max.len+1)
    b.sub <- matrix(rep(0,n.elemb), nrow = max.len)
    
    for (i in c(index.cure.var)) {
      for (j in c(1:length(index.surv.var), max.len+1)) {
        b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])
        #b.sub[i,j] <- -sum((design.matrix1[,i]*design.xt0[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])
        
      }
    }
   # if (k <= length(index.cure.var)) {info.b = b.sub[index.cure.var[-k],c(index.surv.var-max.len,index.gamma-max.len)]} else
   # {info.b = b.sub[index.cure.var,c(index.surv.var[-ik]-max.len,index.gamma-max.len)]}
    info.b = b.sub[index.cure.var,c(index.surv.var-max.len,index.gamma-max.len)]
    
        ###info matrix block d  
    design.xt1 <- cbind(design.matrix1, log(survt[,1]))
    
    n.elemd <- (max.len+1)^2
    d.sub1 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
    d.sub2 <- matrix(rep(0,n.elemd), nrow = (max.len+1))
    
    for (i in c(index.surv.var-max.len, max.len +1)) {
      for (j in c(index.surv.var-max.len, max.len +1)) {
        d.sub1[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*eps)[survt[, 2] == 1])
        d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
        #d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta^2))[survt[, 2] == 0])
      }
    }
    d.sub = d.sub1 + d.sub2 + 
      matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-1]^2)), 
             nrow = (max.len + 1))
    
    # if (k <= length(index.cure.var)) 
    #   {info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]} else
    #   {info.d = d.sub[c(index.surv.var[-ik]-max.len,index.gamma-max.len),c(index.surv.var[-ik]-max.len,index.gamma-max.len)]}
    
    info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]
    
    info.d.inv = mat.inv(info.d)
    
    fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
    #hessian.mat = -fisher.info
    
    # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
    info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)
    
    #determinant of hessian matrix;
    det.info = matrix.det(info.set0)*matrix.det(info.d)
    #det.info = matrix.det(fisher.info)
    
    loglik.part = loglikelihood - 0.5*log(det.info)
    } else
    if (pl == FALSE)
    {
      loglik.part = loglikelihood
    }
    
    return(loglik.part)
  }
  
  
  #################################################################
  #### parameter estimation under H0 for individual parameter
  #### loglikelihood ratio test statistics for each cure part variable;
  
  dim.v <- ncol(design.matrix)

  if (LRT.pval == T) {
    ll.cure <- rep(0,dim.v)
    llr.cure <- rep(0,dim.v)
    pval.cure <- rep(0,dim.v)
    #ll.est.cure <- array(0,c((2*dim.v+1),(2*dim.v+1)))
    #varmat.A.cure <-sample(0,3*dim.v,replace = TRUE)
    #varmat.A.cure <-array(0,c((2*dim.v+1),(2*dim.v+1),(2*dim.v+1))) #vcov matrix of each single parameter (cure.var1,cure.var2) likelihood ratio test;
    
    # dim(varmat.A.cure) = c(dim.v,dim.v,dim.v)
   # score.U.cure <- array(0,c((2*dim.v+1),(2*dim.v+1))) #gradient vector for all variables of each single parameter(cure.var1,cure.var2) likelihood ratio test;
    # init = c(0.866409603,	-0.334994711,	0.004421239,	-0.070689979,	-8.11466538,	-0.349162116,	0.143310052,	-0.049019853,	0.1)
    
    for (k in index.cure.v[-1]) {
   # mle under the reduced (null) model for cure parameter;
    maximizer <- nlm(
      f = loglik.mixture.part, p = init, 
      survt = survt, design.matrix0 = design.matrix, 
      design.matrix1=design.matrix,
      index.cure.var=index.cure.v[-k], pl=pl,
      iterlim = iterlim, hessian=F
    );
   # ll.est.cure[,k] <- maximizer$estimate
    loglik.part = -maximizer$minimum; 
    dif.ll = -2*(loglik.part-loglik);  #loglik is ll under Ha;
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.cure[k]<- loglik.part
    llr.cure[k]<- dif.ll
    pval.cure[k]<- pval
    if (det(maximizer$hessian) < 1e-05) 
      diag(maximizer$hessian) <- diag(maximizer$hessian) + 1e-06
   # varmat.A.cure[,,k] <- solve(maximizer$hessian)  #if hessian matrix contains row/col of zeros, add small values;
   # score.U.cure [,k] <- maximizer$gradient
  
  }
}
    ###################################
    # Profile likelihood CI endpoint  #
    ###################################
  
  # Note:
  # loglik     -- loglikelihood of all the parameters under MLE of full likelihood;
  # l.up       -- loglik corresponds to upper or lower CI bounds B0=B+delta or B0=B-delta;
  # tol        -- tolerance level for defining loglik difference of estimated and actual parameters are converged;
  # lambda     -- lambda quantity for profile likelhood calculation as referenced;
  # delta.up   -- delta quantity for upper endpoint of profile likelihood calculation as referenced;
  # delta.lo   -- delta quantity for lower endpoint of profile likelihood calculation as referenced;
  # l.temp/l0.b-- loglik of estimated endpoint parameter values under profile likelihood for the corresponding parameter;
  
  ######################################################
  ## By parameter upper or lower endpoint calculation ##
  ######################################################
   apct=0.05
  
  l.null = loglik - 0.5 * qchisq(1-apct,df=1,ncp = 0,lower.tail=T)
  ni = 1
  
  ######## Cure part variable CI endpoints ########
  #################################################

  upper.cure <- rep(0,dim.v)
  lower.cure <- rep(0,dim.v)
  for (k in index.cure.v) {
    
  ##################upper endpoint##########################
   # tol = 0.2
    tol = 0.1
 l.temp <- loglik
 
 n=ni+1
 # iter0 <- 1; converge = FALSE;l0.b.up = 0;
 # while(l.temp > l.null & (!converge| is.nan(l0.b.up)) & iter0<=30) {
        
   #assign initial values to parameter estimates
   param.est.up <- maximizer0$estimate
   
   converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; l0.b.up = 0;
    while ((!converge| is.nan(l0.b.up)) & iter1 <= 25 & !EXIT1) {
  
             # calculate log-lik, score and hessian under l0.b;
      maximizer.temp <-  nlm(
        f = loglik.mixture, p = param.est.up, survt=survt, design.matrix=design.matrix,
        pl = pl, iterlim = 1, hessian=TRUE)
      score.temp = maximizer.temp$gradient
      hessian.temp = maximizer.temp$hessian
      if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06
      
      #### Approach 1: lambda = (2*(l0.b-l.up+e*A^-1*U)/(e*A^-1*e))^0.5
      # l0.b.up <- -loglik.mixture(p=param.est.up, survt, design.matrix,
      #                            index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl)
      l0.b.up <- -maximizer.temp$minimum
      inv.hessian.temp <- solve(hessian.temp)
      if (l0.b.up < l.null + 0.5 * score.temp %*% inv.hessian.temp %*% score.temp) {
        lambda <- ((-inv.hessian.temp %*% score.temp)[k])/inv.hessian.temp[k,k]
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda);
        } else{
      lambda <- (2*(l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
      #define increment for estimated value: delta=-A^-1(U-lambda*e)
      delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)};
      
      # maximizing loop for unpenalized estimates;
      #if (pl == F) {
      inside <- FALSE; iter2 <- 1; 
      while (!inside & iter2 <= 100 ) {
        
        # add increment to stepwise parameter value;
        param.est.temp.up <- param.est.up
        param.est.temp.up[k] <- param.est.temp.up[k]+delta.up
       # if (k==2) {param.est.temp.up[1] <- -1;param.est.temp.up[9] <- 0.1}
        
        #compute loglikelihood function using updated parameter values;
       
        maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.up[-k], survt=survt, 
                                param.est = param.est.temp.up[k], k = k,
                                pl = pl, iterlim = iterlim, hessian=TRUE)
        l.temp.up = -maximizer.temp1$minimum
        
        #if (!is.nan(l.temp.up))
          
         #compare to see if updated l is still 
          inside <- (l.temp.up > (l.null - 0.05)) #l.null - 0.05 for all others, 0.2 for k=3 of high rate H0
          #diff.up = l.temp.up - l.null 
          #converge0 <- (abs(diff.up) <= tol)
          alevel.up <- pchisq(2*(l.temp-l.temp.up),df=1,ncp=0,lower.tail = T)
         # print(c(delta.up, alevel.up, n,l.temp.up,k,iter1,iter2))
          if (!inside) {delta.up <- delta.up/((n+1)/n);iter2 <- iter2 + 1}  #(n+0.1)/n for low rate H0;
       } #for iter2
        
        #}
      #Using converged increment for parameter to get corresponding score and variance expressions;
      if (k==1) {
      param.est.up <- c(param.est.temp.up[k], maximizer.temp1$estimate)
      } else
      {
        param.est.up <- c(maximizer.temp1$estimate[1:(k-1)], param.est.temp.up[k], maximizer.temp1$estimate[-c(1:(k-1))])
      }
      l0.b.up = l.temp.up
      
      diff.up = l0.b.up - l.null 
      converge <- (abs(diff.up) <= tol)
      if (!converge| is.nan(l0.b.up)) {iter1 <- iter1 + 1; n = n + 1} else {EXIT1 = T;}
    } #for iter1
  #} #for iter0
    upper.cure[k] <- param.est.up[k]
    
 
  ###############lower endpoint#####################
 
    n=ni
 #iter0 <- 1; converge = FALSE; l0.b.lo=0
 #while(l.temp > l.null & (!converge| is.nan(l0.b.lo)) & iter0<=30) {
   
   #assign initial values to parameter estimates
   param.est.lo <- maximizer0$estimate
   
   converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; l0.b.lo=0
   while ((!converge| is.nan(l0.b.lo)) & iter1 <= 25 & !EXIT1) {
     
     # calculate log-lik, score and hessian under l0.b;
     maximizer.temp <-  nlm(
       f = loglik.mixture, p = param.est.lo, survt=survt, design.matrix=design.matrix,
       pl = pl, iterlim = 1, hessian=TRUE)
     score.temp = maximizer.temp$gradient
     hessian.temp = maximizer.temp$hessian
     if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06
     
     #### Approach 1: lambda = (2*(l0.b-l.null+e*A^-1*U)/(e*A^-1*e))^0.5
     # l0.b.lo <- -loglik.mixture(p=param.est.lo, survt, design.matrix,
     #                            index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl)
     l0.b.lo <- -maximizer.temp$minimum
     inv.hessian.temp <- solve(hessian.temp)
     if ((l0.b.lo < l.null + 0.5 * score.temp %*% inv.hessian.temp %*% score.temp)|(l0.b.lo - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k]<0) {
       lambda <- ((-inv.hessian.temp %*% score.temp)[k])/inv.hessian.temp[k,k]
       #define increment for estimated value: delta=-A^-1(U-lambda*e)
     } else{
       lambda <- -(2*(l0.b.lo - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5}
       #define increment for estimated value: delta=-A^-1(U-lambda*e)
       delta.lo <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda);
       
       # maximizing loop for unpenalized estimates;
       #if (pl == F) {
       inside <- FALSE; iter2 <- 1; 
       while (!inside & iter2 <= 100) {
         
         # add increment to stepwise parameter value;
         param.est.temp.lo <- param.est.lo
         param.est.temp.lo[k] <- param.est.temp.lo[k] + delta.lo
        # param.est.temp.lo[9] <- 0.1
         
         #compute loglikelihood function using lodated parameter values;
         maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.lo[-k], survt=survt, 
                                 param.est = param.est.temp.lo[k], k = k,
                                 pl = pl, iterlim = iterlim, hessian=TRUE)
         l.temp.lo = -maximizer.temp1$minimum

       
           inside <- (l.temp.lo > l.null - 0.1)
          # diff.lo = l.temp.lo - l.null 
          # converge0 <- (abs(diff.lo) <= tol)
           alevel.lo <- pchisq(2*(l.temp-l.temp.lo),df=1,ncp=0,lower.tail = T)
          # print(c(delta.lo, alevel.lo, n,l.temp.lo,k,iter1,iter2))
           if (!inside) {delta.lo <- delta.lo/((n+0.5)/n);iter2 <- iter2 + 1} #for ones not lowrate+H0 (n+0.5)/n;
        
       } # for iter2;
     
     #}
     #Using converged increment for parameter to get corresponding score and variance expressions;
     if (k==1) {
     param.est.lo <- c(param.est.temp.lo[k], maximizer.temp1$estimate)
     } else {
       param.est.lo <- c(maximizer.temp1$estimate[1:(k-1)], param.est.temp.lo[k], maximizer.temp1$estimate[-c(1:(k-1))])
     }
     l0.b.lo = l.temp.lo 
       
       diff.lo = l0.b.lo - l.null 
       converge <- (abs(diff.lo) <= tol)
       if (!converge | is.nan(l0.b.lo)) {iter1 <- iter1 + 1; n = n + 2} else {EXIT1 = T}
   } #for iter1
 
   lower.cure[k] <- param.est.lo[k]
   
 }
 
  

  ### loglikelihood calculation for each surv part variable;
  
 
  if (LRT.pval == T) {
    ll.surv <- rep(0,ncol(design.matrix))
    llr.surv <- rep(0,ncol(design.matrix))
    pval.surv <- rep(0,ncol(design.matrix))
    
     for (k in index.surv.v[-1]) {
    # mle under the reduced (null) model for surv parameter;
    is=k-length(index.cure.v)
    maximizer <- nlm(
      f = loglik.mixture.part, p = init, 
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
  }
  
  ######## Surv part variable CI endpoints ########
  #################################################
 
  upper.surv <- rep(0,dim.v)
  lower.surv <- rep(0,dim.v)
  for (k in index.surv.v) {
    
    is=k-length(index.cure.v)
    ##################upper endpoint##########################
    l.temp <- loglik
    tol = 0.2
    
    n=ni
    # iter0 <- 1; converge = FALSE;l0.b.up = 0;
    # while(l.temp > l.null & (!converge| is.nan(l0.b.up)) & iter0<=30) {
    
    #assign initial values to parameter estimates
    param.est.up <- maximizer0$estimate
    
    converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; l0.b.up = 0;
    while ((!converge| is.nan(l0.b.up)) & iter1 <= 25 & !EXIT1) {
      
      # calculate log-lik, score and hessian under l0.b;
      maximizer.temp <-  nlm(
        f = loglik.mixture, p = param.est.up, survt=survt, design.matrix=design.matrix,
        pl = pl, iterlim = 1, hessian=TRUE)
      score.temp = maximizer.temp$gradient
      hessian.temp = maximizer.temp$hessian
      if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06
      
      #### Approach 1: lambda = (2*(l0.b-l.up+e*A^-1*U)/(e*A^-1*e))^0.5
      # l0.b.up <- -loglik.mixture(p=param.est.up, survt, design.matrix,
      #                            index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl)
      l0.b.up <- -maximizer.temp$minimum
      inv.hessian.temp <- solve(hessian.temp)
      if (l0.b.up < l.null + 0.5 * score.temp %*% inv.hessian.temp %*% score.temp) {
        lambda <- ((-inv.hessian.temp %*% score.temp)[k])/inv.hessian.temp[k,k]
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda);
      } else{
        lambda <- (2*(l0.b.up - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.up <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)};
      
      # maximizing loop for unpenalized estimates;
      #if (pl == F) {
      inside <- FALSE; iter2 <- 1; 
      while (!inside & iter2 <= 100 ) {
        
        # add increment to stepwise parameter value;
        param.est.temp.up <- param.est.up
        param.est.temp.up[k] <- param.est.temp.up[k]+delta.up
        
        #compute loglikelihood function using updated parameter values;
        
        maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.up[-k], survt=survt, 
                                param.est = param.est.temp.up[k], k = k,
                                pl = pl, iterlim = iterlim, hessian=TRUE)
        l.temp.up = -maximizer.temp1$minimum
        
        #if (!is.nan(l.temp.up))
        
        #compare to see if updated l is still 
        inside <- (l.temp.up > (l.null - 0.2))
        #diff.up = l.temp.up - l.null 
        #converge0 <- (abs(diff.up) <= tol)
        alevel.up <- pchisq(2*(l.temp-l.temp.up),df=1,ncp=0,lower.tail = T)
        #print(c(delta.up, alevel.up, n,l.temp.up,k,iter1,iter2))
        if (!inside) {delta.up <- delta.up/((n+1)/n);iter2 <- iter2 + 1} 
      } #for iter2
      
      #}
      #Using converged increment for parameter to get corresponding score and variance expressions;
      if (k==1) {
        param.est.up <- c(param.est.temp.up[k], maximizer.temp1$estimate)
      } else
      {
        param.est.up <- c(maximizer.temp1$estimate[1:(k-1)], param.est.temp.up[k], maximizer.temp1$estimate[-c(1:(k-1))])
      }
      l0.b.up = l.temp.up
      
      diff.up = l0.b.up - l.null 
      converge <- (abs(diff.up) <= tol)
      if (!converge| is.nan(l0.b.up)) {iter1 <- iter1 + 1; n = n + 1} else {EXIT1 = T;}
    } #for iter1
    
    upper.surv[is] <- param.est.up[k]
    
    
    ###############lower endpoint#####################
    n=ni
    #iter0 <- 1; converge = FALSE; l0.b.lo=0
    #while(l.temp > l.null & (!converge| is.nan(l0.b.lo)) & iter0<=30) {
    
    #assign initial values to parameter estimates
    param.est.lo <- maximizer0$estimate
    
    converge <- FALSE; iter1 <- 1; EXIT1 <-FALSE; l0.b.lo=0
    while ((!converge| is.nan(l0.b.lo)) & iter1 <= 25 & !EXIT1) {
      
      # calculate log-lik, score and hessian under l0.b;
      maximizer.temp <-  nlm(
        f = loglik.mixture, p = param.est.lo, survt=survt, design.matrix=design.matrix,
        pl = pl, iterlim = 1, hessian=TRUE)
      score.temp = maximizer.temp$gradient
      hessian.temp = maximizer.temp$hessian
      if (det(hessian.temp) < 1e-05) diag(hessian.temp) <- diag(hessian.temp) + 1e-06
      
      #### Approach 1: lambda = (2*(l0.b-l.null+e*A^-1*U)/(e*A^-1*e))^0.5
      # l0.b.lo <- -loglik.mixture(p=param.est.lo, survt, design.matrix,
      #                            index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl)
      l0.b.lo <- -maximizer.temp$minimum
      inv.hessian.temp <- solve(hessian.temp)
      if (l0.b.lo < l.null + 0.5 * score.temp %*% inv.hessian.temp %*% score.temp) {
        lambda <- ((-inv.hessian.temp %*% score.temp)[k])/inv.hessian.temp[k,k]
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.lo <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda);
      } else{
        lambda <- -(2*(l0.b.lo - l.null + score.temp %*% inv.hessian.temp %*% score.temp)/inv.hessian.temp[k,k])^0.5
        #define increment for estimated value: delta=-A^-1(U-lambda*e)
        delta.lo <- -inv.hessian.temp[k,k] %*% (score.temp[k] - lambda)};
      
      # maximizing loop for unpenalized estimates;
      #if (pl == F) {
      inside <- FALSE; iter2 <- 1; 
      while (!inside & iter2 <= 100) {
        
        # add increment to stepwise parameter value;
        param.est.temp.lo <- param.est.lo
        param.est.temp.lo[k] <- param.est.temp.lo[k] + delta.lo
        
        #compute loglikelihood function using lodated parameter values;
        maximizer.temp1 <- nlm( f = loglik.mixture.profile, p = param.est.temp.lo[-k], survt=survt, 
                                param.est = param.est.temp.lo[k], k = k,
                                pl = pl, iterlim = iterlim, hessian=TRUE)
        l.temp.lo = -maximizer.temp1$minimum
        
        
        inside <- (l.temp.lo > l.null - 0.2)
        # diff.lo = l.temp.lo - l.null 
        # converge0 <- (abs(diff.lo) <= tol)
        alevel.lo <- pchisq(2*(l.temp-l.temp.lo),df=1,ncp=0,lower.tail = T)
       # print(c(delta.lo, alevel.lo, n,l.temp.lo,k,iter1,iter2))
        if (!inside) {delta.lo <- delta.lo/((n+1)/n);iter2 <- iter2 + 1} 
        
      } # for iter2;
      
      #}
      #Using converged increment for parameter to get corresponding score and variance expressions;
      if (k==1) {
        param.est.lo <- c(param.est.temp.lo[k], maximizer.temp1$estimate)
      } else {
        param.est.lo <- c(maximizer.temp1$estimate[1:(k-1)], param.est.temp.lo[k], maximizer.temp1$estimate[-c(1:(k-1))])
      }
      l0.b.lo = l.temp.lo 
      
      diff.lo = l0.b.lo - l.null 
      converge <- (abs(diff.lo) <= tol)
      if (!converge | is.nan(l0.b.lo)) {iter1 <- iter1 + 1; n = n + 1} else {EXIT1 = T}
    } #for iter1
    
    lower.surv[is] <- param.est.lo[k]
  }
  

  coef.table.cure <- cbind(
    'coef'        = maximizer0$estimate[index.cure.v],
    'exp(coef)'   = exp(maximizer0$estimate[index.cure.v]),
    # 'LL.cure' = ll.cure,
    # 'LLR'         = llr.cure,
    # 'Pr(>chisq)'  = pval.cure,
    'LCI.95%' = lower.cure,
    'UCI.95%' = upper.cure
  );
  rownames(coef.table.cure) <- colnames(design.matrix);
  
  coef.table.surv <- cbind(
    'coef'        = maximizer0$estimate[index.surv.v],
    'exp(coef)'   = exp(maximizer0$estimate[index.surv.v]),
    # 'LL.surv' = ll.surv,
    # 'LLR'         = llr.surv,
    # 'Pr(>chisq)'  = pval.surv,
    'LCI.95%' = lower.surv,
    'UCI.95%' = upper.surv
  );
  rownames(coef.table.surv) <- colnames(design.matrix);
  
  coef.table.alpha <- 0;
}

#run.time = proc.time() - init.time


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


