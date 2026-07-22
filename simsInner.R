##Inner functions 

library(dplyr)
library(mvtnorm)
library(MASS)
library(gamm4)
library(mgcv)
library(randomForest)

logit<-function(x){
  log(x/(1-x))
}

expit<-function(x){
  1/(1+exp(-x))
}

##Derivative of expit
d_expit<-function(x){
  expit(x)*(1-expit(x))
}

## Function to force positive semi-definiteness of square matrix x
make_pos_semi_def = function(x){
  x2 = svd(x)
  #remove negative vals
  x2$d=ifelse(x2$d>0, x2$d, 0)
  #remake matrix
  x2=x2$u%*%diag(x2$d)%*%t(x2$v)
  return(x2)
}

## Hall's GFPCA: Smooths mean and covariance using univariate and bivariate p-splines
## Input: Takes NxD matrix, curves, with each row a user's curve
## Output: Returns smooth mean and eigenfunctions, eigenvalues, and estimated covariance
gfpca<- function(curves){
  
  D = dim(curves)[2]
  N = dim(curves)[1]

  tt = seq(0,1, len=D)
  
  #estimate smooth mean function
  dta = data.frame(z=curves%>%t()%>%as.vector(), tt=tt%>%rep(N))
  gam1 <- gam(z~s(tt, bs = "cr", m=2), method="REML", family=binomial(link="logit"), data=dta)
  #because probability needs to ensure values are between 0 and 1
  alpha_t = gam1%>%predict(newdata=data.frame(tt=tt), type="response")%>%as.vector()
  alpha_t = ifelse(alpha_t<0.001, 0.001, alpha_t)
  alpha_t = ifelse(alpha_t>0.999, 0.999, alpha_t)
  
  #from Hall 2008 paper 
  v_t_hat = as.vector(logit(alpha_t))
  alpha_t_hat = as.vector(alpha_t)
  
  beta_ts  = crossprod(curves)/N
  
  #Format for smooth function 
  Bs = cbind(as.vector(beta_ts), rep(1:D, D), 
             as.vector(t(kronecker(matrix(1, ncol=D), 1:D))))
  #remove diagonal values
  Bs[which(Bs[,2]==Bs[,3]), 1] = NA
  Bs[,2:3] = Bs[,2:3]/D
  Bs = data.frame(val = Bs[,1], t1 = Bs[,2], t2 = Bs[,3])
  #smooth the covariance matrix, need large theta because we have lots of data points
  beta_ts2 = gam(val~te(t1, t2, 
                        bs = "cr", m=2), 
                 method = "REML", data = Bs)
  newdat = data.frame(t1 = (1:D)/D, t2 = (1:D)/D)
  diag.vals = predict(beta_ts2, newdata = newdat)
  return.mat = matrix(NA, nrow = D, ncol = D)
  diag(return.mat) = diag.vals
  return.mat[which(Bs[,2]!=Bs[,3])] = beta_ts2$fitted.value
  beta_ts_hat = return.mat  
  
  tau_ts = (beta_ts_hat - tcrossprod(alpha_t_hat))/
    (tcrossprod(d_expit(v_t_hat)))
  
  #make sym and semi-post def
  tau_ts_hat = round(make_pos_semi_def(tau_ts), 5)
  mu_t_hat = v_t_hat
  
  #decomposition of resulting covariance matrix
  #tau_eigens  =  estimate_eigenfunctions(tau_ts_hat, pve = 0.98)
  tau_eigen_funcs = eigen(tau_ts_hat)$vectors*sqrt(D)
  tau_eigen_vals = eigen(tau_ts_hat)$values/D
  
  ##Re-estimating mean function
  KX<-max(2, which(cumsum(tau_eigen_vals)/sum(tau_eigen_vals)>0.99)[1])
  
  dtaX<-data.frame(value = as.vector(t(curves)), id = rep(1:N, rep(D,N)),
                   index = rep(tt, N))
  for (k in 1:KX){
    dtaX<-cbind(dtaX, rep(tau_eigen_funcs[,k], N))
  }
  names(dtaX)[4:(4 + KX - 1)] <- c(paste0("psi", 1:KX))
  random.structure <- paste(paste0("psi", 1:KX), collapse = "+")
  random.formula <- formula(paste("~(0+", random.structure, "|| id)"))
  
  twoStep <- gamm4(value ~ s(index), family = "binomial", data = dtaX, random = random.formula)
  mu_t_hat <- as.vector(predict.gam(twoStep$gam, newdata = data.frame(index = tt)))
  
  return(list(mu=mu_t_hat, eigenfuncs=tau_eigen_funcs, eigenvals=tau_eigen_vals, sigma=tau_ts_hat))
}

## Pointwise weighted likelihood for censored, normally distributed, functional data;
## used by estMeanZInit when calling optim
## Input: Takes parameters beta, a 2-dim vector containing the mean (beta[1]) and 
##        pointwise standard dev. (exp(beta[2])), y an nx1 vector containing all user
##        data to be evaluated at the current time and prob, an nx1 vector containing
##        user specific probability of posting at the current time.
## Output: The negative log-likelihood
ffMargNew<-function(beta, y, prob){
  toRet<-0
  for (i in 1:length(y)){
    if (is.na(y[i])){
      toRet<-toRet
    }
    else if (y[i]==0){
      toRet<-toRet+(log(pnorm(0, mean=beta[1], sd=exp(beta[2])))/prob[i])
    }
    else if (y[i]==1){
      toRet<-toRet+(log(1-pnorm(1, mean=beta[1], sd=exp(beta[2])))/prob[i])
    }
    else{
      toRet<-toRet+(log(dnorm(y[i], mean=beta[1], sd=exp(beta[2])))/prob[i])
    }
  }
  if (toRet==-Inf){
    10000
  } else{
    -toRet
  }
}

## Estimates pointwise mean, variance of censored functional data
## Input: Takes a, nxJ matrix of user data, probEst, nxJ matrix of users' estimated
##        probability of posting, and susInd, giving indicators of user susceptibility
## Output: Returns pointwise estimates of mean, variance for susceptible and control users
estMeanZInit<-function(a, probEst, susInd){
  J<-ncol(a)
  muZInitSus<-rep(0,J)
  sdZInitSus<-rep(0,J)
  muZInitCont<-rep(0,J)
  sdZInitCont<-rep(0,J)
  for (j in 1:J){
    #Makes call to ffMargNew to compute likelihood
    betaOptJSus<-optim(c(0,0), ffMargNew, gr=NULL, a[which(susInd==1),j],  
                       probEst[which(susInd==1),j])$par
    muZInitSus[j]<-betaOptJSus[1]
    sdZInitSus[j]<-betaOptJSus[2]%>%exp()
    
    betaOptJCont<-optim(c(0,0), ffMargNew, gr=NULL, a[which(susInd==0),j],  
                        probEst[which(susInd==0),j])$par
    muZInitCont[j]<-betaOptJCont[1]
    sdZInitCont[j]<-betaOptJCont[2]%>%exp()
  }
  
  return(list(muZInitSus=muZInitSus, muZInitCont=muZInitCont, 
              sdZInitSus=sdZInitSus, sdZInitCont=sdZInitCont))
}

## Smooths pointwise mean of censored functional data for susceptible and control users
## Input: Takes muZInitSus, muZInitCont, the pointwise mean of susceptible and control users
## Output: Returns smooth mean function estimates for susceptible and control users
smthMeanZ<-function(muZInitSus, muZInitCont){
  J<-length(muZInitSus)
  estFrame<-data.frame(t=rep(1:J), mu=muZInitSus)
  muZEstSus <- as.vector(gam(mu~s(t, bs = "cr", m=2), data=estFrame, method="REML")$fitted.values)[1:J]
  
  estFrame<-data.frame(t=rep(1:J), mu=muZInitCont)
  muZEstCont <- as.vector(gam(mu~s(t, bs = "cr", m=2), data=estFrame, method="REML")$fitted.values)[1:J]
  
  return(list(muZEstSus=muZEstSus, muZEstCont=muZEstCont))
}

## Pairwise weighted likelihood for censored, normally distributed, functional data;
## used by estCovZ when calling optim
## Input: Takes parameters rho, the correlation at current pair of time points, 
##        y an nx2 vector containing all user data for current pair of time points,
##        mu and sigma, each two-dim vectors continaing and mean and variance of each
##        time point, and prob1, prob2, each nx1 vectors containing probability of posting
##        for time points 1 and 2 respectively.
## Output: The negative log-likelihood
ffPair<-function(rho, y, mu, sigma, prob1, prob2){
  toRet<-0
  multiVar<-matrix(c(sigma[1]^2, rho*sigma[1]*sigma[2], rho*sigma[1]*sigma[2], sigma[2]^2), nrow=2, ncol=2)
  for (i in 1:nrow(y)){
    if (is.na(y[i,1]) | is.na(y[i,2])){
      toRet<-toRet
    }
    else if (y[i,1]==0 & y[i,2]==0){
      toRet<-toRet+(log(pmvnorm(upper=c(0,0), mean = mu, sigma=multiVar)[1])/(prob1[i]*prob2[i]))
    }
    else if (y[i,1]>0 & y[i,1]<1 & y[i,2]==0){
      toRet<-toRet+((log(pnorm(0, mean=mu[2]+(rho*sigma[2]*(y[i,1]-mu[1])/sigma[1]), 
                               sd=sqrt(1-rho^2)*sigma[2]))+
                       log(dnorm(y[i,1], mean=mu[1], sd=sigma[1])))/(prob1[i]*prob2[i]))
    } 
    else if (y[i,1]==1 & y[i,2]==0){
      toRet<-toRet+(log(pmvnorm(lower=c(1, -Inf), upper=c(Inf,0), mean = mu, sigma=multiVar)[1])/(prob1[i]*prob2[i]))
    }
    else if (y[i,1]==0 & y[i,2]>0 & y[i,2]<1){
      toRet<-toRet+((log(pnorm(0, mean=mu[1]+(rho*sigma[1]*(y[i,2]-mu[2])/sigma[2]), 
                               sd=sqrt(1-rho^2)*sigma[1]))+
                       log(dnorm(y[i,2], mean=mu[2], sd=sigma[2])))/(prob1[i]*prob2[i]))
    }
    else if (y[i,1]>0 & y[i,1]<1 & y[i,2]>0 & y[i,2]<1){
      toRet<-toRet+(dmvnorm(y[i,], mean = mu, sigma=multiVar, log=T)/(prob1[i]*prob2[i]))
    }
    else if (y[i,1]==1 & y[i,2]>0 & y[i,2]<1){
      toRet<-toRet+((log(1-pnorm(1, mean=mu[1]+(rho*sigma[1]*(y[i,2]-mu[2])/sigma[2]), 
                                 sd=sqrt(1-rho^2)*sigma[1]))+
                       log(dnorm(y[i,2], mean=mu[2], sd=sigma[2])))/(prob1[i]*prob2[i]))
    }
    else if (y[i,1]==0 & y[i,2]==1){
      toRet<-toRet+(log(pmvnorm(lower=c(-Inf, 1), upper=c(0,Inf), mean = mu, sigma=multiVar)[1])/(prob1[i]*prob2[i]))
    }
    else if (y[i,2]==1 & y[i,1]>0 & y[i,1]<1){
      toRet<-toRet+((log(1-pnorm(1, mean=mu[2]+(rho*sigma[2]*(y[i,1]-mu[1])/sigma[1]), 
                                 sd=sqrt(1-rho^2)*sigma[2]))+
                       log(dnorm(y[i,1], mean=mu[1], sd=sqrt(sigma[1]))))/(prob1[i]*prob2[i]))
    }
    else{
      toRet<-toRet+(log(pmvnorm(lower=c(1, 1), mean = mu, sigma=multiVar)[1])/(prob1[i]*prob2[i]))
    }
  }
  -toRet
}

## Estimates each pointwise covariance function for censored functional data of susceptible and control users
## Input: Takes a, an nxJ matrix of user data, probEst, an nxJ matrix of user estiamted probability of 
##        posting, muZSus, sdZSus, the pointwise mean and variance estimates for susceptible users,
##        muZCont, sdZCont, the pointwise mean and variance estimates for control users, and susInd, 
##        which gives indicators of susceptibility.
## Output: Returns pointwise covariance function estimates for susceptible and control users
estCovZ<-function(a, probEst, muZSus, sdZSus, muZCont, sdZCont, susInd){
  J<-ncol(a)
  sigmaInitSus<-matrix(nrow = J, ncol=J)
  sigmaInitCont<-matrix(nrow = J, ncol=J)
  for (j1 in seq(1, J-1, by=2)){
    for (j2 in seq(j1, J, by=2)){
      if (j1==j2){
        #Should just be the variance if j1=j2
        sigmaInitSus[j1,j2]<-sdZSus[j1]^2
        sigmaInitCont[j1,j2]<-sdZCont[j1]^2
      }
      else{
        #Estimates correlation for j1,j2; makes call to ffPair to compute pairwise likelihood
        sigmaIter<-optim(0, ffPair, gr=NULL, cbind(a[which(susInd==1),j1], a[which(susInd==1),j2]), 
                         c(muZSus[j1], muZSus[j2]), c(sdZSus[j1], sdZSus[j2]), 
                         probEst[which(susInd==1),j1], probEst[which(susInd==1),j2], 
                         method="Brent", lower=-1, upper=1)$par
        #Multiply correlation by simga at j1, j2
        sigmaInitSus[j1,j2]<-sigmaIter*sdZSus[j1]*sdZSus[j2]
        sigmaInitSus[j2,j1]<-sigmaIter*sdZSus[j1]*sdZSus[j2]
        
        sigmaIter<-optim(0, ffPair, gr=NULL, cbind(a[which(susInd==0),j1], a[which(susInd==0),j2]), 
                         c(muZCont[j1], muZCont[j2]),  c(sdZCont[j1], sdZCont[j2]), 
                         probEst[which(susInd==0),j1], probEst[which(susInd==0),j2], 
                         method="Brent", lower=-1, upper=1)$par
        sigmaInitCont[j1,j2]<-sigmaIter*sdZCont[j1]*sdZCont[j2]
        sigmaInitCont[j2,j1]<-sigmaIter*sdZCont[j1]*sdZCont[j2]
      }
    }
  }
  return(list(sigmaInitSus=sigmaInitSus, sigmaInitCont=sigmaInitCont))
}

## Function to smooth sample covariance and extract eigen-components for censored functional data
## Input: Takes sigmaInit, a sample covariance to be smoothed
## Output: Returns smoothed covariance, eigenfunctions and eigenvalues
smoothAndEigen<-function(sigmaInit){
  J<-dim(sigmaInit)[1]
  
  Bs = cbind(as.vector(sigmaInit), rep(1:J, J)/J, as.vector(t(kronecker(matrix(1, ncol=J) , 1:J)))/J)
  #remove diagonal values
  Bs[which(Bs[,2]==Bs[,3]), 1] = NA
  
  estFrame<-data.frame(t1=Bs[,2], t2=Bs[,3], val=Bs[,1])
  
  sigmaSmthFit <- gam(val~te(t1, t2, 
                             bs = "cr", m=2), 
                      method = "REML", family="gaussian", data = estFrame)
  newdat = data.frame(t1 = (1:J)/J, t2 = (1:J)/J)
  diag.vals = predict(sigmaSmthFit, newdata = newdat)
  sigmaSmth = matrix(predict(sigmaSmthFit, data.frame(t1=Bs[,2], t2=Bs[,3]))%>%as.vector(), nrow=J, ncol=J)
  diag(sigmaSmth) = diag.vals
  
  sigmaSmth = round(make_pos_semi_def(sigmaSmth), 5)
  eigenfuncs<-eigen(sigmaSmth)$vectors*sqrt(J)
  eigenvals<-eigen(sigmaSmth)$values/J
  
  return(list(sigma=sigmaSmth, eigenfuncs=eigenfuncs, eigenvals=eigenvals))
}

## Function to calculate conditional binomial likelihood, given scores, using
## logit link; called by estScoresX, estScoresXTest
## Input: Takes mu, a vector of conditional mean values at each time point, and y
##        a vector of a users data at each time point
## Output: Returns log-likelihood
ffLogitBayes<-function(mu, y){
  toRet<-0
  for (j in 1:length(y)){
    if (is.na(y[j])){
      toRet<-toRet
    }
    else{
      toRet<-toRet+dbinom(y[j], size=1, prob=expit(mu[j]), log=T)
    }
  }
  toRet
}

## Function to calculate conditional normal censored likelihood, given scores;
## called by estScoresZ, estScoresZTest
## Input: Takes mu, a vector of conditional mean values at each time point, and y
##        a vector of a users data at each time point
## Output: Returns log-likelihood
ffMargBayes<-function(mu, y){
  toRet<-0
  for (j in 1:length(y)){
    if (is.na(y[j])){
      toRet<-toRet
    }
    else if (y[j]==0){
      toRet<-toRet+(log(pnorm(0, mean=mu[j], sd=0.25)))
    }
    else if (y[j]==1){
      toRet<-toRet+(log(1-pnorm(1, mean=mu[j], sd=0.25)))
    }
    else{
      toRet<-toRet+(log(dnorm(y[j], mean=mu[j], sd=0.25)))
    }
  }
  toRet
}

## Function to calculate conditional expectation of censored functional data scores, 
## given user data AND group value, using Monte Carlo integration and Bayes rule
## Input: Takes a, an nxJ matrix of user data, muZ, the smooth group-specific mean function, 
##        eigenfuncsZ, a JxK matrix of group-specific eigenfuncs, and eigenvalsZ, 
##        a K-dim vector of group-specific eigenvalues
## Output: Returns nxK matrix of scores, each row a user's specific scores
estScoresZ<-function(a, muZ, eigenfuncsZ, eigenvalsZ){
  J<-dim(eigenfuncsZ)[1]
  K<-dim(eigenfuncsZ)[2]
  n<-dim(a)[1]
  
  blupBayes<-matrix(nrow=n, ncol=K)
  xiBayes<-mvrnorm(10000, mu=rep(0, K), Sigma=diag(eigenvalsZ))
  for (s in 1:n){
    num<-matrix(nrow=K, ncol=10000)
    denom<-0
    for (k in 1:10000){
      #Compute conditional likelihood given kth score draw
      likeliK<-ffMargBayes(muZ+as.vector(eigenfuncsZ%*%xiBayes[k,]), a[s,])
      num[,k]<-xiBayes[k,]*exp(likeliK)
      denom<-denom+exp(likeliK)
    }
    if (denom==0){
      blupBayes[s,]<-colMeans(xiBayes)
    } else{
      blupBayes[s,]<-rowSums(num)/denom
    }
  }
  blupBayes
}

## Function to calculate conditional expectation of posting functional data scores, 
## given user data AND group value, using Monte Carlo integration and Bayes rule
## Input: Takes v, an nxJ matrix of user data, muX, the smooth group-specific mean function, 
##        eigenfuncsX, a JxK matrix of group-specific eigenfuncs, and eigenvalsX, 
##        a K-dim vector of group-specific eigenvalues
## Output: Returns nxK matrix of scores, each row a user's specific scores
estScoresX<-function(v, muX, eigenfuncsX, eigenvalsX){
  J<-dim(eigenfuncsX)[1]
  K<-dim(eigenfuncsX)[2]
  n<-dim(v)[1]
  
  blupBayes<-matrix(nrow=n, ncol=K)
  xiBayes<-mvrnorm(10000, mu=rep(0, K), Sigma=diag(eigenvalsX))
  for (s in 1:n){
    num<-matrix(nrow=K, ncol=10000)
    denom<-0
    for (k in 1:10000){
      #Compute conditional likelihood given kth score draw
      likeliK<-ffLogitBayes(muX+as.vector(eigenfuncsX%*%xiBayes[k,]), v[s,])
      num[,k]<-xiBayes[k,]*exp(likeliK)
      denom<-denom+exp(likeliK)
    }
    if (denom==0){
      blupBayes[s,]<-colMeans(xiBayes)
    } else{
      blupBayes[s,]<-rowSums(num)/denom
    }
  }
  blupBayes
}

## Function to calculate conditional expectation of censored functional data scores, 
## given user data using Monte Carlo integration and Bayes rule, by taking mixture 
## of susceptible and control.
## Input: Takes a, an nxJ matrix of user data, muZ, the marginal smooth mean function, eigenfuncsZ,
##        a JxK matrix of marginal eigenfuncs, scMeanSus, scVarSus, the apriori mean and covariance
##        of susceptible scores, scMeanCont, scVarCont, the apriori mean and covariance
##        of control scores, and pSus, giving p, the mixture probability
## Output: Returns nxK matrix of scores, each row a user's specific scores
estScoresZTest<-function(a, muZ, eigenfuncsZ, scMeanSus, scVarSus, scMeanCont, scVarCont, pSus){
  J<-dim(eigenfuncsZ)[1]
  K<-dim(eigenfuncsZ)[2]
  n<-dim(a)[1]
  
  #Force positive definite
  eigenScVarSus<-eigen(scVarSus)
  eigenScValsSus<-eigenScVarSus$values
  eigenScVecsSus<-eigenScVarSus$vectors
  eigenScValsSus[eigenScValsSus<0]<-0
  
  eigenScVarCont<-eigen(scVarCont)
  eigenScValsCont<-eigenScVarCont$values
  eigenScVecsCont<-eigenScVarCont$vectors
  eigenScValsCont[eigenScValsCont<0]<-0
  
  #Generate mixture
  blupBayes<-matrix(nrow=n, ncol=K)
  xiBayesSus<-mvrnorm(10000, mu=rep(0, length(scMeanSus)), Sigma=diag(rep(1, length(scMeanSus))))
  xiBayesSus<-t(matrix(scMeanSus, nrow=length(scMeanSus), ncol=10000))+
    (xiBayesSus%*%diag(sqrt(eigenScValsSus))%*%t(eigenScVecsSus))
  xiBayesCont<-mvrnorm(10000, mu=rep(0, length(scMeanCont)), Sigma=diag(rep(1, length(scMeanCont))))
  xiBayesCont<-t(matrix(scMeanCont, nrow=length(scMeanCont), ncol=10000))+
    (xiBayesCont%*%diag(sqrt(eigenScValsCont))%*%t(eigenScVecsCont))
  mixInd<-rbinom(10000, size=1, prob=pSus)
  xiBayes<-mixInd*xiBayesSus+(1-mixInd)*xiBayesCont
  for (s in 1:n){
    num<-matrix(nrow=K, ncol=10000)
    denom<-0
    for (k in 1:10000){
      #Compute conditional likelihood given kth score draw
      likeliK<-ffMargBayes(muZ+as.vector(eigenfuncsZ%*%xiBayes[k,]), a[s,])
      num[,k]<-xiBayes[k,]*exp(likeliK)
      denom<-denom+exp(likeliK)
    }
    if (denom==0){
      blupBayes[s,]<-colMeans(xiBayes)
    } else{
      blupBayes[s,]<-rowSums(num)/denom
    }
  }
  blupBayes
}

## Function to calculate conditional expectation of posting functional data scores, 
## given user data using Monte Carlo integration and Bayes rule, by taking mixture 
## of susceptible and control.
## Input: Takes v, an nxJ matrix of user data, muX, the marginal smooth mean function, eigenfuncsX,
##        a JxK matrix of marginal eigenfuncs, scMeanSus, scVarSus, the apriori mean and covariance
##        of susceptible scores, scMeanCont, scVarCont, the apriori mean and covariance
##        of control scores, and pSus, giving p, the mixture probability
## Output: Returns nxK matrix of scores, each row a user's specific scores
estScoresXTest<-function(v, muX, eigenfuncsX, scMeanSus, scVarSus, scMeanCont, scVarCont, pSus){
  J<-dim(eigenfuncsX)[1]
  K<-dim(eigenfuncsX)[2]
  n<-dim(v)[1]
  
  #Force positive definite
  eigenScVarSus<-eigen(scVarSus)
  eigenScValsSus<-eigenScVarSus$values
  eigenScVecsSus<-eigenScVarSus$vectors
  eigenScValsSus[eigenScValsSus<0]<-0
  
  eigenScVarCont<-eigen(scVarCont)
  eigenScValsCont<-eigenScVarCont$values
  eigenScVecsCont<-eigenScVarCont$vectors
  eigenScValsCont[eigenScValsCont<0]<-0
  
  #Generate mixture
  blupBayes<-matrix(nrow=n, ncol=K)
  xiBayesSus<-mvrnorm(10000, mu=rep(0, length(scMeanSus)), Sigma=diag(rep(1, length(scMeanSus))))
  xiBayesSus<-t(matrix(scMeanSus, nrow=length(scMeanSus), ncol=10000))+
    (xiBayesSus%*%diag(sqrt(eigenScValsSus))%*%t(eigenScVecsSus))
  xiBayesCont<-mvrnorm(10000, mu=rep(0, length(scMeanCont)), Sigma=diag(rep(1, length(scMeanCont))))
  xiBayesCont<-t(matrix(scMeanCont, nrow=length(scMeanCont), ncol=10000))+
    (xiBayesCont%*%diag(sqrt(eigenScValsCont))%*%t(eigenScVecsCont))
  mixInd<-rbinom(10000, size=1, prob=pSus)
  xiBayes<-mixInd*xiBayesSus+(1-mixInd)*xiBayesCont
  for (s in 1:n){
    num<-matrix(nrow=K, ncol=10000)
    denom<-0
    for (k in 1:10000){
      #Compute conditional likelihood given kth score draw
      likeliK<-ffLogitBayes(muX+as.vector(eigenfuncsX%*%xiBayes[k,]), v[s,])
      num[,k]<-xiBayes[k,]*exp(likeliK)
      denom<-denom+exp(likeliK)
    }
    if (denom==0){
      blupBayes[s,]<-colMeans(xiBayes)
    } else{
      blupBayes[s,]<-rowSums(num)/denom
    }
  }
  blupBayes
}

## Function to fit random forest model imputing missing data
## Input: Takes impDataNaive and impDataTestNaive, the raw training and test data,
##        as well as J, the number of time points
## Output: Returns prediction values for test set
fitRFNaive<-function(impDataNaive, impDataTestNaive, J){
  ##Impute binary with mode, and continuous with median, of training set
  #VMeds<-apply(impDataNaive$V, 2,
   #            function(x) unique(x[!is.na(x)])[which.max(tabulate(match(x[!is.na(x)], unique(x[!is.na(x)]))))])
  AMeds<-apply(impDataNaive$A, 2, median, na.rm=T)
  
  impDataTestNaive$A<-matrix(1:J)%>%apply(1, function(x) ifelse(is.na(impDataTestNaive$A[,x]),
                                                                AMeds[x], impDataTestNaive$A[,x]))
  impDataNaive$A<-matrix(1:J)%>%apply(1, function(x) ifelse(is.na(impDataNaive$A[,x]),
                                                            AMeds[x], impDataNaive$A[,x]))
  #impDataTestNaive$V<-matrix(1:J)%>%apply(1, function(x) ifelse(is.na(impDataTestNaive$V[,x]),
  #                                                              VMeds[x], impDataTestNaive$V[,x]))
  #impDataNaive$V<-matrix(1:J)%>%apply(1, function(x) ifelse(is.na(impDataNaive$V[,x]),
  #                                                          VMeds[x], impDataNaive$V[,x]))
  
  if ("twt"%in%names(impDataNaive)){
    wrdMeds<-apply(impDataNaive$wrd, 2, median, na.rm=T)
    twtMeds<-apply(impDataNaive$twt, 2, median, na.rm=T)
    impDataTestNaive$wrd<-matrix(1:J)%>%apply(1, function(x) ifelse(is.na(impDataTestNaive$wrd[,x]),
                                                                    wrdMeds[x], impDataTestNaive$wrd[,x]))
    impDataNaive$wrd<-matrix(1:J)%>%apply(1, function(x) ifelse(is.na(impDataNaive$wrd[,x]),
                                                                wrdMeds[x], impDataNaive$wrd[,x]))
    impDataTestNaive$twt<-matrix(1:J)%>%apply(1, function(x) ifelse(is.na(impDataTestNaive$twt[,x]),
                                                                    twtMeds[x], impDataTestNaive$twt[,x]))
    impDataNaive$twt<-matrix(1:J)%>%apply(1, function(x) ifelse(is.na(impDataNaive$twt[,x]),
                                                                twtMeds[x], impDataNaive$twt[,x]))
  }
  
  #Have to fit separately if matching, since state is a factor and the rest are numeric
  if ("state"%in%names(impDataNaive)){
    impDataTmp<-impDataNaive%>%subset(select=names(impDataNaive)[!names(impDataNaive)%in%c("treat", "state")])%>%as.matrix()%>%as.data.frame()
    impDataTmp$state<-impDataNaive$state
    impDataTmp$treat<-impDataNaive$treat
    impDataNaive<-impDataTmp
    
    impDataTestTmp<-impDataTestNaive%>%subset(select=names(impDataTestNaive)[!names(impDataTestNaive)%in%c("treat", "state")])%>%as.matrix()%>%as.data.frame()
    impDataTestTmp$state<-impDataTestNaive$state
    impDataTestNaive<-impDataTestTmp
    
    fitRFNaive<-randomForest(y=as.factor(impDataNaive$treat), x=impDataNaive%>%subset(select=names(impDataNaive)[!names(impDataNaive)%in%c("treat")]))
    pred<-as.vector(predict(fitRFNaive, newdata = impDataTestNaive, type="prob")[,2])
  } else{
    fitRFNaive<-randomForest(y=as.factor(impDataNaive$treat), x=impDataNaive%>%subset(select=names(impDataNaive)[!names(impDataNaive)%in%c("treat")])%>%as.matrix())
    pred<-as.vector(predict(fitRFNaive, newdata = impDataTestNaive%>%as.matrix(), type="prob")[,2])
  }
  pred
}