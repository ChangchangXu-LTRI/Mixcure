#######################################
### Function for LRT of mixcure model##
#### penalized loglikelihoods #########
###################################################
#### Last modified Dec31 2018 for one x variable ##
###################################################

mixcure.penal.lrt.test.uni <- function(formula, data, init, pl, iterlim = 200) { 
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
  
  ## Test for univariate model LRT ##
  ###################################
  
  # setwd('C:\\Biostats_PhD\\Missing imputation studies\\simulation data')
  # data1<-read.csv("junk2.csv", header=TRUE)
  # data2<-read.csv("ki67 data.csv", header=TRUE)
  # data3<-read.csv("EGFR data.csv", header=TRUE)
  # data4<-merge(data1,data2,by="siteid",all.x=TRUE)
  # primdata<-merge(data3,data4,by="siteid")
  # 
  # #primdata <- primdata[complete.cases(primdata),]
  # 
  # primdata$labER<-factor(with(primdata,ifelse((NESTLAB==1),1,
  #                                             ifelse((NESTLAB==2 | NESTLAB==3), 2, 
  #                                                    ifelse((NESTLAB==9),3,".")))))
  # 
  # primdata$labPR<-factor(with(primdata,ifelse((NPROLAB==1),1,
  #                                             ifelse((NPROLAB==2 | NPROLAB==3), 2, 
  #                                                    ifelse((NPROLAB==9),3,".")))))
  # 
  # primdata$MENS0<-factor(with(primdata,ifelse((MENSTAT==3),0,1)))##pre&peri vs. post
  # primdata$TUMCT2<-factor(with(primdata,ifelse((TUMCAT==5),1,0)))##>5 vs. rest
  # primdata$TUMCT5<-factor(with(primdata,ifelse((TUMCAT==4),1,0)))##2-5 vs. rest
  # primdata$TUMCT <-factor(with(primdata,ifelse((TUMCT2==1|TUMCT5==1),1,0)))##>2 vs. <2
  # 
  # primdata$Her2<-factor(with(primdata,ifelse((Cockt==1),1,
  #                                            ifelse((Cockt==0), 0, "."))))
  # 
  # primdata$TN<-factor(with(primdata,
  #                          ifelse((Cockt==0 & labER == 2 & labPR == 2),1,0)))
  # 
  # primdata$Basal<-factor(with(primdata,
  #                             ifelse((Cockt==0 & labER == 2 & labPR == 2 &(CK_5 == 1 | EGFR_NWD == 1)),1,0)))
  # 
  # primdata$TNPnonbasal<-factor(with(primdata,
  #                                   ifelse((Cockt==0 & labER == 2 & labPR == 2 &(CK_5 == 0 & EGFR_NWD == 0)),1,0)))
  # 
  # primdata$LuminalA<-factor(with(primdata,
  #                                ifelse((Cockt==0 & (labER == 1 | labPR == 1) & ki67D == 0),1,0)))
  # 
  # primdata$LuminalB<-factor(with(primdata,
  #                                ifelse((Cockt==0 & (labER == 1 | labPR == 1) & ki67D == 1),1,0)))
  # 
  # primdata$Subtype<-factor(with(primdata,
  #                               ifelse((LuminalA=='1'),"LuminalA",
  #                                      ifelse((LuminalB=='1'),"LuminalB",
  #                                             ifelse((Her2=='1'), "Her2",
  #                                                    ifelse((TN=='1'), "Triple Negative", "Missing"))))))
  # missing1 <- primdata[!primdata$Subtype=="Missing",]
  # missing3 <- missing1[,c('TSURV2','CENS','TN','LuminalA', 'MENS0','TUMCT')]
  # missing3$Her2 <- factor(with(missing1,ifelse((Her2==1),1,0)))
  # 
  # formula=Surv(TSURV2, CENS == 0) ~ Her2; data=missing3; loglik=-625; init=c(-1,rep(0,1),-1,rep(0,1),0.1); pl=T; iterlim=200
  
  
  
  #############################################################################
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
# loglik of full model
  
  loglik.mixture <- function(p, survt, design.matrix, index.cure.var=index.cure.v, index.surv.var=index.surv.v, pl, PLCI=F) {
    
    ####  parameter and variable dep parameters;
    #####
    theta = 1/(1+exp(-design.matrix%*%p[index.cure.var]))
    eps = survt[,1]^(p[index.gamma])*exp(design.matrix%*%p[index.surv.var])
    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap = -theta*(1-theta)*(1-eta)+(1-theta)^2*eta*(1-eta) #for LRT, p-large
    #kap= (1-eta)*(1-theta)*(theta + eta)   # exp for est and PLCI ; LRT p-small
    pi = exp(eps)*eps*eta^2
    lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
    
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
          #b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*theta*(1-theta)*pi)[survt[, 2] == 0]) #for LRT
          b.sub[i,j] <- -sum((design.matrix[,i]*design.xt[,j]*eps*(1-delta)*delta)[survt[, 2] == 0]) #for est
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
          d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta-eps^2*delta+eps^2*delta^2))[survt[, 2] == 0])
          #d.sub2[i,j] <- sum((design.xt[,i]*design.xt[,j]*(eps*delta^2))[survt[, 2] == 0])
          
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
    iterlim = iterlim, hessian=F);
  loglik <- -maximizer0$minimum  #in loglik function loglik was calculated as minus of actual loglik value
  
  
  #############################################################

  #create CI for profile likelihood, this option only outputs estimates and PL under specified model;

  
  ##loglik function for testing parameters of cure or surv part;
  loglik.mixture.part <- function(p, survt, design.matrix1, design.matrix0, part.cure,index.cure.var=index.cure.v,index.surv.var=index.surv.v, pl) {  
 
 ## Test function using univariate model ##
 ############################################
    # k=2; op.est=c(1.571611,0, -6.691835,  0.919785,  1.614850); design.matrix1=design.matrix; design.matrix0=design.matrix; 
    # p=op.est[-k]; index.cure.var=index.cure.v[-k];index.surv.var=index.surv.v; part.cure=T
    
      #design.matrix1-surv, design.matrix0-cure
    design.mtx.comb = cbind(design.matrix0,design.matrix1)
 
    #parameter and variable dep parameters;
    if (k > length(index.cure.v)) {
      theta = 1/(1+exp(-design.matrix[,index.cure.var]%*%as.matrix(p[index.cure.var])))
      eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]*p[-c(index.cure.var,index.gamma-1)])
                                  } else {
      theta = 1/(1+exp(-design.matrix[,index.cure.var]*p[-c(index.surv.var-1,index.gamma-1)]))
      eps = survt[,1]^(p[index.gamma-1])*exp(design.mtx.comb[,index.surv.var]%*%as.matrix(p[index.surv.var-1]))
                                         }

    eta = 1/((exp(eps)-1)*theta+1)
    delta = 1/(theta/(1-theta)*exp(eps)+1)
    kap = -theta*(1-theta)*(1-eta)+(1-theta)^2*eta*(1-eta) #for LRT, p-large
    #kap= (1-eta)*(1-theta)*(theta + eta) #for est, PLCI; LRT p-small
    pi = exp(eps)*eps*eta^2
    #lambda = (1-theta)^2*eta*(1-eta)*((2*eta-1)*(1-theta)+3)
    #phi = theta*(1-theta)*((2*eta-1)*(1-theta)+theta)*pi
    
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
    
    if (pl==F) {loglik.part = loglikelihood} else  {
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
      for (j in c(index.surv.var-max.len, max.len+1)) {
        #b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*theta*(1-theta)*pi)[survt[, 2] == 0])  #equivalent to expression below; *eps*(1-delta)*delta
        b.sub[i,j] <- -sum((as.matrix(design.matrix1)[,i]*design.xt0[,j]*eps*(1-delta)*delta)[survt[, 2] == 0])
        #b.sub[i,j] <- -sum((design.matrix1[,i]*design.xt0[,j]*eps*(1-eta)*eta*(1-theta))[survt[, 2] == 0])
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
        #d.sub2[i,j] <- sum((design.xt1[,i]*design.xt1[,j]*(eps*delta^2))[survt[, 2] == 0])
        
                                                       }
                                                     }
    d.sub = d.sub1 + d.sub2 + 
      matrix(c(rep(0, (n.elemd - 1)),sum(survt[, 2] == 1)/(p[index.gamma-1]^2)), 
             nrow = (max.len + 1))
    
    info.d = d.sub[c(index.surv.var-max.len,index.gamma-max.len),c(index.surv.var-max.len,index.gamma-max.len)]
    
    
    info.d.inv = mat.inv(info.d)
    
    #fisher.info = rbind(cbind(info.a,info.b),cbind(t(info.b),info.d))
    #hessian.mat = -fisher.info
    
    # #info.set0 is (A-BD^-1B^T), dif than used in modified score;
    if (part.cure == T) {info.set0 = info.a-t(info.b)%*%info.d.inv%*%info.b} else
      {info.set0 = info.a-info.b%*%info.d.inv%*%t(info.b)}
    
    #determinant of hessian matrix;
    det.info = matrix.det(info.set0)*matrix.det(info.d)
    
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
 # ll.est.cure <- array(0,c((2*dim.v+1),(2*dim.v+1)))
  #varmat.A.cure <-sample(0,3*dim.v,replace = TRUE)
  #varmat.A.cure <-array(0,c((2*dim.v+1),(2*dim.v+1),(2*dim.v+1))) #vcov matrix of each single parameter (cure.var1,cure.var2) likelihood ratio test;
  
 # dim(varmat.A.cure) = c(dim.v,dim.v,dim.v)
  #score.U.cure <- array(0,c((2*dim.v+1),(2*dim.v+1))) #gradient vector for all variables of each single parameter(cure.var1,cure.var2) likelihood ratio test;
  
  #init1=c(-0.1,rep(0,5),-5,rep(0,5),0.1)
  
  for (k in index.cure.v[-1]) {
   # mle under the reduced (null) model for cure parameter;
  #  if (k==5|k==6) {init1=c(-1,rep(0,5),-5,rep(-0.1,5),0.1)}
  #  if (k==6) {init1=c(-0.5,rep(0,5),-5,rep(0,5),0.1)}
     maximizer <- nlm(
      f = loglik.mixture.part, part.cure = T,
      p = init[-k], 
      survt = survt, design.matrix0 = design.matrix, 
      design.matrix1=design.matrix,
      index.cure.var=index.cure.v[-k], pl=pl,
      iterlim = iterlim, hessian=F
    );
  #  ll.est.cure[,k] <- maximizer$estimate
    loglik.part = -maximizer$minimum; 
    dif.ll = -2*(loglik.part-loglik);  #loglik is ll under Ha;
    pval = pchisq(abs(dif.ll),df=1,lower.tail=FALSE);
    ll.cure[k]<- loglik.part
    llr.cure[k]<- dif.ll
    pval.cure[k]<- pval
    # if (det(maximizer$hessian) < 1e-05) 
    #   diag(maximizer$hessian) <- diag(maximizer$hessian) + 1e-06
   # varmat.A.cure[,,k] <- solve(maximizer$hessian)  #if hessian matrix contains row/col of zeros, add small values;
   # score.U.cure [,k] <- maximizer$gradient
  
  }

 
  
  ### loglikelihood calculation for each surv part variable;
  
  ll.surv <- rep(0,ncol(design.matrix))
  llr.surv <- rep(0,ncol(design.matrix))
  pval.surv <- rep(0,ncol(design.matrix))
  
  #init=c(-0.1,rep(-0.1,1),-10,rep(-0.1,1),0.1)
  
  for (k in index.surv.v[-1]) {
    # mle under the reduced (null) model for surv parameter;
    is=k-length(index.cure.v)
    maximizer <- nlm(
      f = loglik.mixture.part, p =  init[-k], part.cure = F,  
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


#######################################
## Output tables from either method; ##
#######################################

# colnames(var.mat) <- c(
#   paste('cure.', colnames(design.matrix)), 
#   paste('surv.', colnames(design.matrix)),
#   'alpha'
# );
# rownames(var.mat) <- colnames(var.mat);

out <- list(
  coefficients = list(
    cure = coef.table.cure, 
    surv = coef.table.surv, 
    alpha = coef.table.alpha
 #   run.time
  )
  #cov = var.mat
);
class(out) <- c('mixcure', 'list');

return(out);

}


# #### print.mixcure #############################################################
# # DESCRIPTION
# #   To print a mixcure object.
# # INPUT
# #   object : a mixcure object, which is an outcome of function mixcure.
# #   digits : number of digits for printing, passed to print.default.
# #   ...    : other parameters passed to print.default.
# # OUTPUT
# #   NULL.   
# 
# print.mixcure <- function(object, digits = 3, ...) {
#   sep.line.cure   <- paste(c(rep('-', 37),  ' CURE ' , rep('-', 37)), collapse = '');
#   sep.line.surv   <- paste(c(rep('-', 37),  ' SURVIVAL ' , rep('-', 37)), collapse = '');
#   sep.line.alpha <- paste(c(rep('-', 36), ' ALPHA ', rep('-', 36)), collapse = '');
#   
#   message(sep.line.cure);
#   print.default(object$coefficients$cure,   digits = digits, ...);
#   
#   message(sep.line.surv);
#   print.default(object$coefficients$surv,   digits = digits, ...);
#   
#   message(sep.line.alpha);
#   print.default(object$coefficients$alpha, digits = digits, ...);
#   
#   return(NULL);
# };
# 
# 
# #### coef.mixcure ##############################################################
# coef.mixcure <- function(object) {
#   coefs <- c(
#     object$coefficients$cure[, 'coef'], 
#     object$coefficients$surv[, 'coef'], 
#     object$coefficients$alpha[, 'coef']
#   );
#   names(coefs) <- c( paste('cure.', rownames(object$coefficients$cure)), 
#                      paste('surv.', rownames(object$coefficients$surv)),
#                      rownames(object$coefficients$alpha) );
#   
#   return(coefs);
# }
# 
# 
# vcov.mixcure <- function(object) {
#   return(object$cov);
# }


