##Data application 

library(dplyr)
library(mgcv)
library(refund)
library(randomForest)
library(PRROC)
library(survival)

#Takes inputs 1...4800
# args=commandArgs(TRUE)
# input<-as.numeric(args[1])
input<-3506
source("~/TwitterProj/fileSub/simsInner.R")
totRawData<-readRDS("~/TwitterProj/fileSub/finalActData.rds")

##AInit is proportion of ats in tweet, vInit posting indicators
twtInitSus<-totRawData$twtInitSus
wrdInitSus<-totRawData$wrdInitSus
vInitSus<-totRawData$vInitSus
AInitSus<-totRawData$AInitSus

twtInitCont<-totRawData$twtInitCont
wrdInitCont<-totRawData$wrdInitCont
vInitCont<-totRawData$vInitCont
AInitCont<-totRawData$AInitCont

stateInit<-factor(totRawData$state)
joinDateInit<-as.numeric(as.Date(totRawData$joinDate))
joinDateInit<-joinDateInit-min(joinDateInit)+1

seed0<-235
nTest<-500

#ns give training set sample size
#rs determine JEnds, JStart; b gives MC replicate 1 to 100
ns<-apply(matrix(c(100, 500, 1000)), 1, rep, 1600)%>%as.vector()
rs<-apply(matrix(rep(c(1:16),3)), 1, rep, 100)%>%as.vector()
bs<-c(1:100)%>%rep(48)

n<-ns[input]
r<-rs[input]
b<-bs[input]

#JEnds are number of weeks in the future to predict susceptibility
#JStarts are number of weeks observed
JEnds<-c(rep(0, 4), rep(4, 4), rep(8, 4), rep(12, 4))
JStarts<-rep(c(28, 36, 44, 52), 4)

##Accuracy
praucRFNaive<-list()
praucRFNew<-list()
praucLinear<-list()
praucLinearCor<-list()
praucLinearWt<-list()
praucGAMUni<-list()
praucGAMMulti<-list()
praucCLogit<-list()
praucCLogitCor<-list()

#Store coefficients of LogiCenFD (M1)
coefEsts<-list()

set.seed(seed0*(b+1))
if (n==1000){
  #Sample from initial distribution of states when matching to include all states
  initStateFreq<-floor(table(stateInit)/length(stateInit)*(n/2))
  trainSamp<-c()
  for (s in 1:length(initStateFreq)){
    trainSamp<-trainSamp%>%append(sample(which(stateInit==names(initStateFreq)[s]), 
                                         size=as.vector(initStateFreq)[s]))
  }
  if (length(trainSamp)<n/2){
    trainSamp<-trainSamp%>%append(sample(c(1:nrow(AInitSus))[!(c(1:nrow(AInitSus))%in%trainSamp)], 
                                         size = n/2-length(trainSamp)))
  }
  
  aSus<-AInitSus[trainSamp, ]
  aCont<-AInitCont[trainSamp, ]
  vSus<-vInitSus[trainSamp,]
  vCont<-vInitCont[trainSamp,]
  twtSus<-twtInitSus[trainSamp,]
  twtCont<-twtInitCont[trainSamp,]
  wrdSus<-wrdInitSus[trainSamp,]
  wrdCont<-wrdInitCont[trainSamp,]
  
  a<-rbind(aSus, aCont)
  v<-rbind(vSus, vCont)
  twt<-rbind(twtSus, twtCont)
  wrd<-rbind(wrdSus, wrdCont)
  
  stateTrain<-stateInit[trainSamp]
  joinDateTrain<-joinDateInit[trainSamp]
} else{
  trainSamp<-sample(1:nrow(AInitSus), size = n/2)
  
  aSus<-AInitSus[trainSamp, ]
  aCont<-AInitCont[trainSamp, ]
  vSus<-vInitSus[trainSamp,]
  vCont<-vInitCont[trainSamp,]
  twtSus<-twtInitSus[trainSamp,]
  twtCont<-twtInitCont[trainSamp,]
  wrdSus<-wrdInitSus[trainSamp,]
  wrdCont<-wrdInitCont[trainSamp,]
  
  a<-rbind(aSus, aCont)
  v<-rbind(vSus, vCont)
  twt<-rbind(twtSus, twtCont)
  wrd<-rbind(wrdSus, wrdCont)
  
  stateTrain<-stateInit[trainSamp]
  joinDateTrain<-joinDateInit[trainSamp]
}

JStart<-JStarts[r]
JEnd<-JEnds[r]

##Sample training data to be JStart weeks so that in at most JEnd weeks, potential 
## interaction occurs
if (JEnd==0){
  trainEnds <- rep(52, n)
  trainSamples <- apply(matrix(trainEnds), 1, function(s) (s-JStart+1):s)%>%t()
} else if (JStart+JEnd>52){
  JEndTmp<-52-JStart
  if (JEndTmp==0){
    trainEnds <- rep(52, n)
    trainSamples <- apply(matrix(trainEnds), 1, function(s) (s-JStart+1):s)%>%t()
  } else{
    trainEnds <- sample((52-JEndTmp):52, size = n, replace = T)
    trainSamples <- apply(matrix(trainEnds), 1, function(s) (s-JStart+1):s)%>%t()
  }
} else{
  trainEnds <- sample((52-JEnd):52, size = n, replace = T)
  trainSamples <- apply(matrix(trainEnds), 1, function(s) (s-JStart+1):s)%>%t()
}

susInd<-c(rep(1, n/2), rep(0, n/2))

aTrain<-apply(matrix(1:n), 1,  function(i) a[i,trainSamples[i,]])%>%t()
vTrain<-apply(matrix(1:n), 1,  function(i) v[i,trainSamples[i,]])%>%t()
twtTrain<-apply(matrix(1:n), 1,  function(i) twt[i,trainSamples[i,]])%>%t()
wrdTrain<-apply(matrix(1:n), 1,  function(i) wrd[i,trainSamples[i,]])%>%t()

J<-JStart

##Step 1. Estimate X
vGFPCASus<-gfpca(vTrain[which(susInd==1),])
muXEstSus<-vGFPCASus$mu
eigenfuncsXSus<-vGFPCASus$eigenfuncs
eigenvalsXSus<-vGFPCASus$eigenvals
sigmaSmthSusX<-vGFPCASus$sigma

vGFPCACont<-gfpca(vTrain[which(susInd==0),])
muXEstCont<-vGFPCACont$mu
eigenfuncsXCont<-vGFPCACont$eigenfuncs
eigenvalsXCont<-vGFPCACont$eigenvals
sigmaSmthContX<-vGFPCACont$sigma

#Compute mixture mean and covariance
muXEst<-(muXEstSus*mean(susInd))+(muXEstCont*mean(1-susInd))
sigmaSmthX<-(sigmaSmthSusX*mean(susInd))+(sigmaSmthContX*mean(1-susInd))+
  (tcrossprod(muXEstSus-muXEst)*mean(susInd))+(tcrossprod(muXEstCont-muXEst)*mean(1-susInd))
eigenfuncsX<-eigen(sigmaSmthX)$vectors*sqrt(J)
eigenvalsX<-eigen(sigmaSmthX)$values/J
KX<-max(2, which(cumsum(eigenvalsX)/sum(eigenvalsX)>0.99)[1])

#Compute conditional expectation of scores
scMeanSusX<-as.vector((muXEstSus-muXEst)%*%eigenfuncsX[,1:KX])/J
scMeanContX<-as.vector((muXEstCont-muXEst)%*%eigenfuncsX[,1:KX])/J
scVarSusX<-t(eigenfuncsX[,1:KX])%*%sigmaSmthSusX%*%eigenfuncsX[,1:KX]/J^2
scVarContX<-t(eigenfuncsX[,1:KX])%*%sigmaSmthContX%*%eigenfuncsX[,1:KX]/J^2

scX<-estScoresXTest(vTrain, muXEst, eigenfuncsX[,1:KX], scMeanSusX, scVarSusX, scMeanContX, scVarContX, mean(susInd))

xEst<-t(matrix(muXEst, nrow=J, ncol=n))+scX%*%t(eigenfuncsX[,1:KX])
probEst<-expit(xEst) #Estimated probability of posting

#Reweighting scores for MFPCA to be on same scale.
WX<-diag(sigmaSmthX)
WX[WX<=0]<-10^(-8)
WX<-diag(1/sqrt(WX))
wtXCurves<-scX%*%t(eigenfuncsX[,1:KX])%*%WX
wtXEigenfuncs<-eigen(cov(wtXCurves))$vectors[,1:KX]*sqrt(J)
scX<-(wtXCurves%*%wtXEigenfuncs)/J

##Step 2. Estimate Z
#Estimate pointwise mean, variance
meanEstZInit<-estMeanZInit(aTrain, probEst, susInd)
muZInitSus<-meanEstZInit$muZInitSus
muZInitCont<-meanEstZInit$muZInitCont
sdZInitSus<-meanEstZInit$sdZInitSus
sdZInitCont<-meanEstZInit$sdZInitCont

#Smooth mean
meanEstZ<-smthMeanZ(muZInitSus, muZInitCont)
muZEstSus<-meanEstZ$muZEstSus
muZEstCont<-meanEstZ$muZEstCont

#Estimate pointwise covariance
covEstZInit<-estCovZ(aTrain, probEst, muZInitSus, sdZInitSus, muZInitCont, sdZInitCont, susInd)
sigmaInitSus<-covEstZInit$sigmaInitSus
sigmaInitCont<-covEstZInit$sigmaInitCont

#Smooth and eigendecomp covariance
eigenSus<-smoothAndEigen(sigmaInitSus)
eigenfuncsZSus<-eigenSus$eigenfuncs
eigenvalsZSus<-eigenSus$eigenvals
sigmaSmthSusZ<-eigenSus$sigma

eigenCont<-smoothAndEigen(sigmaInitCont)
eigenfuncsZCont<-eigenSus$eigenfuncs
eigenvalsZCont<-eigenSus$eigenvals
sigmaSmthContZ<-eigenCont$sigma


#Compute mixture mean and covariance
muZEst<-(muZEstSus*mean(susInd))+(muZEstCont*mean(1-susInd))
sigmaSmthZ<-(sigmaSmthSusZ*mean(susInd))+(sigmaSmthContZ*mean(1-susInd))+
  (tcrossprod(muZEstSus-muZEst)*mean(susInd))+(tcrossprod(muZEstCont-muZEst)*mean(1-susInd))
eigenfuncsZ<-eigen(sigmaSmthZ)$vectors*sqrt(J)
eigenvalsZ<-eigen(sigmaSmthZ)$values/J
KZ<-max(2,which(cumsum(eigenvalsZ)/sum(eigenvalsZ)>0.99)[1])

#Compute conditional expectation of scores
scMeanSusZ<-as.vector((muZEstSus-muZEst)%*%eigenfuncsZ[,1:KZ])/J
scMeanContZ<-as.vector((muZEstCont-muZEst)%*%eigenfuncsZ[,1:KZ])/J
scVarSusZ<-t(eigenfuncsZ[,1:KZ])%*%sigmaSmthSusZ%*%eigenfuncsZ[,1:KZ]/J^2
scVarContZ<-t(eigenfuncsZ[,1:KZ])%*%sigmaSmthContZ%*%eigenfuncsZ[,1:KZ]/J^2

scZ<-estScoresZTest(aTrain, muZEst, eigenfuncsZ[,1:KZ], scMeanSusZ, scVarSusZ, scMeanContZ, scVarContZ, mean(susInd))

#Reweighting scores for MFPCA to be on same scale.
WZ<-diag(sigmaSmthZ)
WZ[WZ<=0]<-10^(-8)
WZ<-diag(1/sqrt(WZ))
wtZCurves<-scZ%*%t(eigenfuncsZ[,1:KZ])%*%WZ
wtZEigenfuncs<-eigen(cov(wtZCurves))$vectors[,1:KZ]*sqrt(J)
scZ<-(wtZCurves%*%wtZEigenfuncs)/J

#Step 3. Estimate tweets
fpcaTwtSus<-fpca.sc(Y=twtTrain[which(susInd==1),], Y.pred=twtTrain, pve=0.99)
fpcaTwtCont<-fpca.sc(Y=twtTrain[which(susInd==0),], Y.pred=twtTrain, pve=0.99)

muTwtEst<-(fpcaTwtSus$mu*mean(susInd))+(fpcaTwtCont$mu*mean(1-susInd))
if (length(fpcaTwtSus$evalues)==1){
  sigmaSmthSusTwt<-fpcaTwtSus$evalues*tcrossprod(fpcaTwtSus$efunctions)
} else{
  sigmaSmthSusTwt<-fpcaTwtSus$efunctions%*%diag(fpcaTwtSus$evalues)%*%t(fpcaTwtSus$efunctions)
}
if (length(fpcaTwtCont$evalues)==1){
  sigmaSmthContTwt<-fpcaTwtCont$evalues*tcrossprod(fpcaTwtCont$efunctions)
} else{
  sigmaSmthContTwt<-fpcaTwtCont$efunctions%*%diag(fpcaTwtCont$evalues)%*%t(fpcaTwtCont$efunctions)
}
sigmaSmthTwt<-(sigmaSmthSusTwt*mean(susInd))+(sigmaSmthContTwt*mean(1-susInd))+
  (tcrossprod(fpcaTwtSus$mu-muTwtEst)*mean(susInd))+(tcrossprod(fpcaTwtCont$mu-muTwtEst)*mean(1-susInd))
eigenfuncsTwt<-eigen(sigmaSmthTwt)$vectors*sqrt(J)
eigenvalsTwt<-eigen(sigmaSmthTwt)$values/J
KTwt<-min(max(2,which(cumsum(eigenvalsTwt)/sum(eigenvalsTwt)>0.99)[1]), 
          ncol(fpcaTwtSus$efunctions), ncol(fpcaTwtCont$efunctions))

scTwtSus<-(fpcaTwtSus$Yhat-t(matrix(muTwtEst, nrow=J, ncol=n)))%*%eigenfuncsTwt[,1:KTwt]/J
scTwtCont<-(fpcaTwtCont$Yhat-t(matrix(muTwtEst, nrow=J, ncol=n)))%*%eigenfuncsTwt[,1:KTwt]/J
scTwt<-(scTwtSus*mean(susInd))+(scTwtCont*mean(1-susInd))

WTwt<-diag(sigmaSmthTwt)
WTwt[WTwt<=0]<=10^(-8)
WTwt<-diag(1/sqrt(WTwt))
wtTwtCurves<-scTwt%*%t(eigenfuncsTwt[,1:KTwt])%*%WTwt
wtTwtEigenfuncs<-eigen(cov(wtTwtCurves))$vectors[,1:KTwt]*sqrt(J)
scTwt<-(wtTwtCurves%*%wtTwtEigenfuncs)/J

#Step 4. Estimate words
fpcaWrdSus<-fpca.sc(Y=wrdTrain[which(susInd==1),], Y.pred=wrdTrain, pve=0.99)
fpcaWrdCont<-fpca.sc(Y=wrdTrain[which(susInd==0),], Y.pred=wrdTrain, pve=0.99)

muWrdEst<-(fpcaWrdSus$mu*mean(susInd))+(fpcaWrdCont$mu*mean(1-susInd))
if (length(fpcaWrdSus$evalues)==1){
  sigmaSmthSusWrd<-fpcaWrdSus$evalues*tcrossprod(fpcaWrdSus$efunctions)
} else{
  sigmaSmthSusWrd<-fpcaWrdSus$efunctions%*%diag(fpcaWrdSus$evalues)%*%t(fpcaWrdSus$efunctions)
}
if (length(fpcaWrdCont$evalues)==1){
  sigmaSmthContWrd<-fpcaWrdCont$evalues*tcrossprod(fpcaWrdCont$efunctions)
} else{
  sigmaSmthContWrd<-fpcaWrdCont$efunctions%*%diag(fpcaWrdCont$evalues)%*%t(fpcaWrdCont$efunctions)
}
sigmaSmthWrd<-(sigmaSmthSusWrd*mean(susInd))+(sigmaSmthContWrd*mean(1-susInd))+
  (tcrossprod(fpcaWrdSus$mu-muWrdEst)*mean(susInd))+(tcrossprod(fpcaWrdCont$mu-muWrdEst)*mean(1-susInd))
eigenfuncsWrd<-eigen(sigmaSmthWrd)$vectors*sqrt(J)
eigenvalsWrd<-eigen(sigmaSmthWrd)$values/J
KWrd<-min(max(2,which(cumsum(eigenvalsWrd)/sum(eigenvalsWrd)>0.99)[1]), 
          ncol(fpcaWrdSus$efunctions), ncol(fpcaWrdCont$efunctions))

scWrdSus<-(fpcaWrdSus$Yhat-t(matrix(muWrdEst, nrow=J, ncol=n)))%*%eigenfuncsWrd[,1:KWrd]/J
scWrdCont<-(fpcaWrdCont$Yhat-t(matrix(muWrdEst, nrow=J, ncol=n)))%*%eigenfuncsWrd[,1:KWrd]/J
scWrd<-(scWrdSus*mean(susInd))+(scWrdCont*mean(1-susInd))

WWrd<-diag(sigmaSmthWrd)
WWrd[WWrd<=0]<=10^(-8)
WWrd<-diag(1/sqrt(WWrd))
wtWrdCurves<-scWrd%*%t(eigenfuncsWrd[,1:KWrd])%*%WWrd
wtWrdEigenfuncs<-eigen(cov(wtWrdCurves))$vectors[,1:KWrd]*sqrt(J)
scWrd<-(wtWrdCurves%*%wtWrdEigenfuncs)/J

#Step 5: MFPCA
scEigen<-cbind(scZ, scX, scTwt, scWrd)%>%cov()%>%eigen()
scVals<-scEigen$values
scVecs<-scEigen$vectors
KSc<-max(2, which(cumsum(scVals)/sum(scVals)<0.99)%>%length())
KScAlt<-max(2, which(cumsum(scVals)/sum(scVals)<0.95)%>%length()) #For LogiCenFD(M3)

scTot<-cbind(scZ, scX, scTwt, scWrd)%*%scVecs
scTot<-scTot%>%apply(2, scale, center=F)

#Step 6: Fit models
trainFrame<-data.frame(treat=susInd, state=rep(stateTrain, 2), 
                       joinDate=rep(joinDateTrain, 2))
for (k in 1:KSc){
  trainFrame<-cbind(trainFrame, scTot[,k])
}
names(trainFrame)[4:(4 + KSc - 1)] <- c(paste0("sc", 1:KSc))
trainFrame$V<-vTrain
trainFrame$A<-aTrain
trainFrame$twt<-twtTrain
trainFrame$wrd<-wrdTrain
trainFrame$id<-rep(1:(n/2), 2)

##3. Linear
gamForm <- formula(paste0("treat~", paste(paste0("sc", 1:KSc), collapse = "+")))
fitLinear<-glm(gamForm, data=trainFrame, family = binomial(link="logit"))

coefEsts$unmatch<-as.vector(fitLinear$coefficients)

##3a: Linear Corrected
fitLinearCor<-fitLinear
fitLinearCor$coefficients[1]<-fitLinear$coefficients[1]-logit(mean(susInd))

coefEsts$cor<-as.vector(fitLinearCor$coefficients)

##3b. Linear Weighted
gamForm <- formula(paste0("treat~", paste(paste0("sc", 1:KSc), collapse = "+")))
fitLinearWt<-glm(gamForm, data=trainFrame, family = binomial(link="logit"), 
                 weights = c(rep(mean(susInd), sum(trainFrame$treat)), 
                             rep(mean(1-susInd), sum(1-trainFrame$treat))))

coefEsts$wt<-as.vector(fitLinearWt$coefficients)

##4. GAM Uni
gamForm <- formula(paste0("treat~", paste(paste0("s(sc", 1:KSc, ", bs='cr', m=2)"), collapse = "+")))
fitGAMUni<-gam(gamForm, data=trainFrame, family=binomial(link="logit"), method="REML")

##5. GAM Multi
gamForm <- formula(paste0("treat~s(", paste(paste0("sc", 1:KScAlt), collapse = ","), ", k=10, bs='tp')"))
fitGAMMulti<-gam(gamForm, data=trainFrame, family=binomial(link="logit"), method="REML")

##6. Conditional Logistic
gamForm <- formula(paste0("treat~strata(id)+", paste(paste0("sc", 1:KSc), collapse = "+")))
fitCLogit<-clogit(gamForm, data=trainFrame)

#Matched models
if (n==1000){
  ##3. Linear
  gamForm <- formula(paste0("treat~", paste(paste0("sc", 1:KSc), collapse = "+"), "+state+joinDate"))
  fitLinearMatch<-glm(gamForm, data=trainFrame, family = binomial(link="logit"))
  
  coefEsts$match<-as.vector(fitLinearMatch$coefficients)

  ##3a: Linear Corrected
  fitLinearCorMatch<-fitLinearMatch
  fitLinearCorMatch$coefficients[1]<-fitLinearMatch$coefficients[1]-logit(mean(susInd))
  
  coefEsts$matchCor<-as.vector(fitLinearCorMatch$coefficients)

  ##3b. Linear Weighted
  gamForm <- formula(paste0("treat~", paste(paste0("sc", 1:KSc), collapse = "+"), "+state+joinDate"))
  fitLinearWtMatch<-glm(gamForm, data=trainFrame, family = binomial(link="logit"), 
                        weights = c(rep(mean(susInd), sum(trainFrame$treat)), 
                                    rep(mean(1-susInd), sum(1-trainFrame$treat))))
  
  coefEsts$matchWt<-as.vector(fitLinearWtMatch$coefficients)

  ##4. GAM Uni
  gamForm <- formula(paste0("treat~", paste(paste0("s(sc", 1:KSc, ", bs='cr', m=2)"), collapse = "+"), "+state+joinDate"))
  fitGAMUniMatch<-gam(gamForm, data=trainFrame, family=binomial(link="logit"), method="REML")

  ##5. GAM Multi
  gamForm <- formula(paste0("treat~s(", paste(paste0("sc", 1:KScAlt), collapse = ","), ", k=10, bs='tp')+state+joinDate"))
  fitGAMMultiMatch<-gam(gamForm, data=trainFrame, family=binomial(link="logit"), method="REML")
}

##Testing
testSamp<-sample(c(1:nrow(AInitSus))[!(c(1:nrow(AInitSus))%in%trainSamp)], size = nTest/2)
aSusTest<-AInitSus[testSamp,]
aContTest<-AInitCont[testSamp,]
vSusTest<-vInitSus[testSamp,]
vContTest<-vInitCont[testSamp,]
twtSusTest<-twtInitSus[testSamp,]
twtContTest<-twtInitCont[testSamp,]
wrdSusTest<-wrdInitSus[testSamp,]
wrdContTest<-wrdInitCont[testSamp,]

aTest<-rbind(aSusTest, aContTest)
vTest<-rbind(vSusTest, vContTest)
twtTest<-rbind(twtSusTest, twtContTest)
wrdTest<-rbind(wrdSusTest, wrdContTest)

stateTest<-stateInit[testSamp]
joinDateTest<-joinDateInit[testSamp]

#Observe Jstart weeks, if JStart+JEnd>=52 and Yi=1, then susceptible; otherwise control
testStarts <- sample(1:(53-JStart), size = nTest, replace = T)
testSamples <- apply(matrix(testStarts), 1, function(s) s:(s + JStart - 1))%>%t()

susIndTest<-c(ifelse(testSamples[1:(nTest/2), JStart]+JEnd>=52, 1, 0), rep(0, nTest/2))

while (all(susIndTest==0)){
  testStarts <- sample(1:(53-JStart), size = nTest, replace = T)
  testSamples <- apply(matrix(testStarts), 1, function(s) s:(s + JStart - 1))%>%t()
  
  susIndTest<-c(ifelse(testSamples[1:(nTest/2), JStart]+JEnd>=52, 1, 0), rep(0, nTest/2))
}

aTest<-apply(matrix(1:nTest), 1,  function(i) aTest[i,testSamples[i,]])%>%t()
vTest<-apply(matrix(1:nTest), 1,  function(i) vTest[i,testSamples[i,]])%>%t()
twtTest<-apply(matrix(1:nTest), 1,  function(i) twtTest[i,testSamples[i,]])%>%t()
wrdTest<-apply(matrix(1:nTest), 1,  function(i) wrdTest[i,testSamples[i,]])%>%t()

#Estimate scores
scZTest<-estScoresZTest(aTest, muZEst, eigenfuncsZ[,1:KZ], scMeanSusZ, scVarSusZ, scMeanContZ, scVarContZ, mean(susInd))
wtZCurvesTest<-scZTest%*%t(eigenfuncsZ[,1:KZ])%*%WZ
scZTest<-(wtZCurvesTest%*%wtZEigenfuncs)/J

scXTest<-estScoresXTest(vTest, muXEst, eigenfuncsX[,1:KX], scMeanSusX, scVarSusX, scMeanContX, scVarContX, mean(susInd))
wtXCurvesTest<-scXTest%*%t(eigenfuncsX[,1:KX])%*%WX
scXTest<-(wtXCurvesTest%*%wtXEigenfuncs)/J

scTwtTestSus<-(fpca.sc(Y=twtTrain[which(susInd==1),], Y.pred=twtTest, pve=0.99)$Yhat-
                 t(matrix(muTwtEst, nrow=J, ncol=nTest)))%*%eigenfuncsTwt[,1:KTwt]/J
scTwtTestCont<-(fpca.sc(Y=twtTrain[which(susInd==0),], Y.pred=twtTest, pve=0.99)$Yhat-
                  t(matrix(muTwtEst, nrow=J, ncol=nTest)))%*%eigenfuncsTwt[,1:KTwt]/J
scTwtTest<-(scTwtTestSus*mean(susInd))+(scTwtTestCont*mean(1-susInd))
wtTwtCurvesTest<-scTwtTest%*%t(eigenfuncsTwt[,1:KTwt])%*%WTwt
scTwtTest<-(wtTwtCurvesTest%*%wtTwtEigenfuncs)/J

scWrdTestSus<-(fpca.sc(Y=wrdTrain[which(susInd==1),], Y.pred=wrdTest, pve=0.99)$Yhat-
                 t(matrix(muWrdEst, nrow=J, ncol=nTest)))%*%eigenfuncsWrd[,1:KWrd]/J
scWrdTestCont<-(fpca.sc(Y=wrdTrain[which(susInd==0),], Y.pred=wrdTest, pve=0.99)$Yhat-
                  t(matrix(muWrdEst, nrow=J, ncol=nTest)))%*%eigenfuncsWrd[,1:KWrd]/J
scWrdTest<-(scWrdTestSus*mean(susInd))+(scWrdTestCont*mean(1-susInd))
wtWrdCurvesTest<-scWrdTest%*%t(eigenfuncsWrd[,1:KWrd])%*%WWrd
scWrdTest<-(wtWrdCurvesTest%*%wtWrdEigenfuncs)/J

scTestTot<-cbind(scZTest, scXTest, scTwtTest, scWrdTest)%*%scVecs
scTestTot<-scTestTot%>%apply(2, scale, center=F)

#Prediction
testFrame<-data.frame(treat=susIndTest, state=rep(stateTest, 2),
                      joinDate=rep(joinDateTest, 2))
for (k in 1:KSc){
  testFrame<-cbind(testFrame, scTestTot[,k])
}
names(testFrame)[4:(4 + KSc - 1)] <- c(paste0("sc", 1:KSc))
testFrame$V<-vTest
testFrame$A<-aTest
testFrame$twt<-twtTest
testFrame$wrd<-wrdTest
testFrame$id<-rep(1:(nTest/2), 2)

##1. Naive Random Forest
impDataNaive<-trainFrame%>%subset(select=c("treat", "A", "twt", "wrd"))
impDataTestNaive<-testFrame%>%subset(select=c("A", "twt", "wrd"))
pred<-fitRFNaive(impDataNaive, impDataTestNaive, J)
praucRFNaive$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                               scores.class1 = as.vector(pred[testFrame$treat==0]),
                               curve=T)

##2. New Random Forest
fitRFNew<-randomForest(y=as.factor(trainFrame$treat),
                       x=trainFrame%>%subset(select=paste0("sc", 1:KSc))%>%as.matrix())
pred<-as.vector(predict(fitRFNew, newdata = testFrame%>%subset(select=paste0("sc", 1:KSc)), type="prob")[,2])
praucRFNew$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                     scores.class1 = as.vector(pred[testFrame$treat==0]),
                                     curve=T)

##3. Linear
pred<-as.vector(predict(fitLinear, newdata = testFrame, type = "response"))
praucLinear$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                      scores.class1 = as.vector(pred[testFrame$treat==0]),
                                      curve=T)

##3a. Linear Corrected
pred<-as.vector(predict(fitLinearCor, newdata = testFrame, type = "response"))
praucLinearCor$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                         scores.class1 = as.vector(pred[testFrame$treat==0]),
                                         curve=T)
##3b. Linear Weighted
pred<-as.vector(predict(fitLinearWt, newdata = testFrame, type = "response"))
praucLinearWt$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                        scores.class1 = as.vector(pred[testFrame$treat==0]),
                                        curve=T)

##4. GAM Uni
pred<-as.vector(predict(fitGAMUni, newdata = testFrame, type = "response"))
praucGAMUni$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                      scores.class1 = as.vector(pred[testFrame$treat==0]),
                                      curve=T)
##5. GAM Multi
pred<-as.vector(predict(fitGAMMulti, newdata = testFrame, type = "response"))
praucGAMMulti$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                           scores.class1 = as.vector(pred[testFrame$treat==0]),
                                           curve=T)

##6. Conditional Logistic
fitLinear$coefficients[2:length(fitLinear$coefficients)]<-fitCLogit$coefficients
pred<-as.vector(predict(fitLinear, newdata = testFrame))
praucCLogit$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                            scores.class1 = as.vector(pred[testFrame$treat==0]),
                            curve=T)

##6a. Conditional Logistic Corrected
fitLinearCor$coefficients[2:length(fitLinearCor$coefficients)]<-fitCLogit$coefficients
pred<-as.vector(predict(fitLinearCor, newdata = testFrame))
praucCLogitCor$unmatch<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                              scores.class1 = as.vector(pred[testFrame$treat==0]),
                              curve=T)



#Matching models
if (n==1000){
  ##1. Naive Random Forest
  impDataNaive<-trainFrame%>%subset(select=c("treat", "A", "twt", "wrd", "state", "joinDate"))
  impDataTestNaive<-testFrame%>%subset(select=c("A", "twt", "wrd", "state", "joinDate"))
  pred<-fitRFNaive(impDataNaive, impDataTestNaive, J)
  praucRFNaive$match<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                 scores.class1 = as.vector(pred[testFrame$treat==0]),
                                 curve=T)
  
  ##2. New Random Forest
  fitRFNewMatch<-randomForest(y=as.factor(trainFrame$treat),
                              x=trainFrame%>%subset(select=c(paste0("sc", 1:KSc), "state", "joinDate")))
  pred<-as.vector(predict(fitRFNewMatch, newdata = testFrame%>%subset(select=c(paste0("sc", 1:KSc), "state", "joinDate")), type="prob")[,2])
  praucRFNew$match<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                             scores.class1 = as.vector(pred[testFrame$treat==0]),
                             curve=T)
  
  ##3. Linear
  pred<-as.vector(predict(fitLinearMatch, newdata = testFrame, type = "response"))
  praucLinear$match<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                              scores.class1 = as.vector(pred[testFrame$treat==0]),
                              curve=T)
  
  ##3a. Linear Corrected
  pred<-as.vector(predict(fitLinearCorMatch, newdata = testFrame, type = "response"))
  praucLinearCor$match<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                 scores.class1 = as.vector(pred[testFrame$treat==0]),
                                 curve=T)
  
  ##3b. Linear Weighted
  pred<-as.vector(predict(fitLinearWtMatch, newdata = testFrame, type = "response"))
  praucLinearWt$match<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                scores.class1 = as.vector(pred[testFrame$treat==0]),
                                curve=T)
  
  ##4. GAM Uni
  pred<-as.vector(predict(fitGAMUniMatch, newdata = testFrame, type = "response"))
  praucGAMUni$match<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                              scores.class1 = as.vector(pred[testFrame$treat==0]),
                              curve=T)
  
  ##5. GAM Multi
  pred<-as.vector(predict(fitGAMMultiMatch, newdata = testFrame, type = "response"))
  praucGAMMulti$match<-pr.curve(scores.class0 = as.vector(pred[testFrame$treat==1]),
                                scores.class1 = as.vector(pred[testFrame$treat==0]),
                                curve=T)
}

saveRDS(list(praucRFNaive=praucRFNaive, praucRFNew=praucRFNew, praucLinear=praucLinear, 
             praucLinearCor=praucLinearCor, praucLinearWt=praucLinearWt,
             praucGAMUni=praucGAMUni, praucGAMMulti=praucGAMMulti,
             praucCLogit=praucCLogit, praucCLogitCor=praucCLogitCor, coefEsts=coefEsts),
        file=paste0("dataAppResults", input, ".rds"))
