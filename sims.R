##Simulations

library(dplyr)
library(mgcv)
library(randomForest)
library(PRROC)

#Takes inputs 1...6000
args=commandArgs(TRUE)
input<-as.numeric(args[1])
source("simsInner.R")

seed0<-235
B<-100
J<-52
nTest<-500

tt<-seq(0,1,length=J)
#Generate with fourier basis, scores with variance lambInit
fourier<-cbind(rep(1, J), sqrt(2)*cos(2*pi*seq(0,1,length=J)),
               sqrt(2)*sin(2*pi*seq(0,1,length=J)))/sqrt(2)
lambInit<-c(0.1, 0.05, 0.025)

#ns give training set sample size, deltas determine mean difference of susceptible and control
#mnews determine number of weeks observed, b gives MC replicate 1...100
ns<-matrix(c(100, 500, 1000))%>%apply(1, function(x) rep(x, 2000))%>%as.vector()
deltas<-matrix(seq(0,2,length=5))%>%apply(1, function(x) rep(x, 3))%>%as.vector()%>%rep(400)
mnews<-matrix(c(28, 36, 44, 52))%>%apply(1, function(x) rep(x, 15))%>%as.vector()%>%rep(100)
bs<-c(1:100)%>%rep(60)

n<-ns[input]
mnew<-mnews[input]
delta<-deltas[input]
b<-bs[input]

muZCont<-rep(expit(-0.4), J)
muZSus<-expit(0.2*(delta-2)+delta*tt^4)
muXCont<-logit(muZCont)
muXSus<-logit(muZSus)

muZ<-0.5*(muZSus+muZCont)
muX<-0.5*(muXSus+muXCont)

#Mean of susceptible, control marginal scores
meanSus<-(c(muZSus-muZ, muXSus-muX)%*%rbind(fourier, fourier))/J
meanCont<-(c(muZCont-muZ, muXCont-muX)%*%rbind(fourier, fourier))/J

#Actual marginal covariance
covInit<-(rbind(fourier, fourier)%*%diag(lambInit)%*%t(rbind(fourier, fourier)))+
  0.25*tcrossprod(c(muZSus, muXSus)-c(muZCont, muXCont))
eigenInit<-eigen(covInit)$vectors[,1:3]*sqrt(J)

#Oracle coefficients for LogiCenFD (M1)
orcCoefs<-(meanSus-meanCont)/lambInit

susInd<-c(rep(1, n/2), rep(0, n/2))

#Generate data
set.seed(seed0*(b+1))
xiSus<-matrix(1:3)%>%apply(1, function(k) rnorm(n/2, mean=meanSus[k], sd=sqrt(lambInit[k])))
xiCont<-matrix(1:3)%>%apply(1, function(k) rnorm(n/2, mean=meanCont[k], sd=sqrt(lambInit[k])))

xiTot<-rbind(xiSus, xiCont)

sc1Gen<-xiTot[,1]
sc2Gen<-xiTot[,2]
sc3Gen<-xiTot[,3]

zSus<-t(matrix(muZ, nrow=J, ncol=n/2))+xiSus%*%t(fourier)
zCont<-t(matrix(muZ, nrow=J, ncol=n/2))+xiCont%*%t(fourier)

xSus<-t(matrix(muX, nrow=J, ncol=n/2))+xiSus%*%t(fourier)
xCont<-t(matrix(muX, nrow=J, ncol=n/2))+xiCont%*%t(fourier)

epsilonSus<-rnorm(n/2*J, sd=0.25)
epsilonCont<-rnorm(n/2*J, sd=0.25)

aSus<-zSus+matrix(epsilonSus, nrow=n/2, ncol=J)
aCont<-zCont+matrix(epsilonCont, nrow=n/2, ncol=J)
a<-rbind(aSus, aCont)
a[a<0]<-0
a[a>1]<-1

vSus<-matrix(rbinom(n/2*J, size=1, prob=expit(as.vector(xSus))), nrow=n/2, ncol=J)
vCont<-matrix(rbinom(n/2*J, size=1, prob=expit(as.vector(xCont))), nrow=n/2, ncol=J)
v<-rbind(vSus, vCont)
a[!v]<-NA

##Step 1. Estimate X
vGFPCASus<-gfpca(vSus)
muXEstSus<-vGFPCASus$mu
eigenfuncsXSus<-vGFPCASus$eigenfuncs
eigenvalsXSus<-vGFPCASus$eigenvals
sigmaSmthSusX<-vGFPCASus$sigma

vGFPCACont<-gfpca(vCont)
muXEstCont<-vGFPCACont$mu
eigenfuncsXCont<-vGFPCACont$eigenfuncs
eigenvalsXCont<-vGFPCACont$eigenvals
sigmaSmthContX<-vGFPCACont$sigma

#Compute mixture mean and covariance
muXEst<-(muXEstSus+muXEstCont)/2
sigmaSmthX<-((sigmaSmthSusX+sigmaSmthContX)/2)+(0.25*(tcrossprod(muXEstSus-muXEstCont)))
eigenfuncsX<-eigen(sigmaSmthX)$vectors*sqrt(J)
KX<-3

#Compute group-specific conditional expectation of scores for recovery
scXSus<-estScoresX(vSus, muXEstSus, eigenfuncsXSus[,1:KX], eigenvalsXSus[1:KX])
xEstSus<-t(matrix(muXEstSus, nrow=J, ncol=n/2))+scXSus%*%t(eigenfuncsXSus[,1:KX])

scXCont<-estScoresX(vCont, muXEstCont, eigenfuncsXCont[,1:KX], eigenvalsXCont[1:KX])
xEstCont<-t(matrix(muXEstCont, nrow=J, ncol=n/2))+scXCont%*%t(eigenfuncsXCont[,1:KX])

muRecX<-muXEst

probEst<-expit(rbind(xEstSus, xEstCont)) #Estimated probability of posting

##Step 2. Estimate Z
#Estimate pointwise mean, variance
meanEstZInit<-estMeanZInit(a, probEst, susInd)
muZInitSus<-meanEstZInit$muZInitSus
muZInitCont<-meanEstZInit$muZInitCont
sdZInitSus<-meanEstZInit$sdZInitSus
sdZInitCont<-meanEstZInit$sdZInitCont

#Smooth mean
meanEstZ<-smthMeanZ(muZInitSus, muZInitCont)
muZEstSus<-meanEstZ$muZEstSus
muZEstCont<-meanEstZ$muZEstCont

##Estimate pointwise covariance
covEstZInit<-estCovZ(a, probEst, muZInitSus, sdZInitSus, muZInitCont, sdZInitCont, susInd)
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

#Compute mixture mean
muZEst<-(muZEstSus+muZEstCont)/2
muRecZ<-muZEst

#Compute mixture covariance
sigmaSmthZ<-((sigmaSmthSusZ+sigmaSmthContZ)/2)+(0.25*(tcrossprod(muZEstSus-muZEstCont)))
eigenfuncsZ<-eigen(sigmaSmthZ)$vectors*sqrt(J)
KZ<-3

#Compute group-specific conditional expectation of scores for recovery
scZSus<-estScoresZ(aSus, muZEstSus, eigenfuncsZSus[,1:KZ], eigenvalsZSus[1:KZ])
scZCont<-estScoresZ(aCont, muZEstCont, eigenfuncsZCont[,1:KZ], eigenvalsZCont[1:KZ])

#Step 3: MFPCA
scEigen<-cbind(rbind(scZSus, scZCont), rbind(scXSus, scXCont))%>%cov()%>%eigen()
scVals<-scEigen$values
scVecs<-scEigen$vectors
KSc<-3

#Recover eigenfunctions with correct sign
eigenTot<-rbind(eigenfuncsZ[,1:KZ]%*%scVecs[1:KSc,], eigenfuncsX[,1:KX]%*%scVecs[1:KSc+KZ,])

zMin1<-which.min(c(mean((eigenTot[1:J,1]-eigenInit[1:J,1])^2), mean((eigenTot[1:J,1]+eigenInit[1:J,1])^2)))
xMin1<-which.min(c(mean((eigenTot[1:J+J,1]-eigenInit[1:J+J,1])^2), mean((eigenTot[1:J+J,1]+eigenInit[1:J+J,1])^2)))
zMin2<-which.min(c(mean((eigenTot[1:J,2]-eigenInit[1:J,2])^2), mean((eigenTot[1:J,2]+eigenInit[1:J,2])^2)))
xMin2<-which.min(c(mean((eigenTot[1:J+J,2]-eigenInit[1:J+J,2])^2), mean((eigenTot[1:J+J,2]+eigenInit[1:J+J,2])^2)))
zMin3<-which.min(c(mean((eigenTot[1:J,3]-eigenInit[1:J,3])^2), mean((eigenTot[1:J,3]+eigenInit[1:J,3])^2)))
xMin3<-which.min(c(mean((eigenTot[1:J+J,3]-eigenInit[1:J+J,3])^2), mean((eigenTot[1:J+J,3]+eigenInit[1:J+J,3])^2)))

psi1Rec<-c(eigenTot[1:J,1]*(-1)*((2*(zMin1-1))-1),eigenTot[1:J+J,1]*(-1)*((2*(xMin1-1))-1))
psi2Rec<-c(eigenTot[1:J,2]*(-1)*((2*(zMin2-1))-1),eigenTot[1:J+J,2]*(-1)*((2*(xMin2-1))-1))
psi3Rec<-c(eigenTot[1:J,3]*(-1)*((2*(zMin3-1))-1),eigenTot[1:J+J,3]*(-1)*((2*(xMin3-1))-1))

#Recover scores with correct sign
scTot<-cbind(rbind(scZSus, scZCont), rbind(scXSus, scXCont))%*%scVecs

sc1Rec<-scTot[,1]*(-1)*((2*(which.min(c(mean((scTot[,1]-xiTot[,1])^2), 
                                            mean((scTot[,1]+xiTot[,1])^2)))-1))-1)
sc2Rec<-scTot[,2]*(-1)*((2*(which.min(c(mean((scTot[,2]-xiTot[,2])^2), 
                                            mean((scTot[,2]+xiTot[,2])^2)))-1))-1)
sc3Rec<-scTot[,3]*(-1)*((2*(which.min(c(mean((scTot[,3]-xiTot[,3])^2), 
                                            mean((scTot[,3]+xiTot[,3])^2)))-1))-1)

#Recover curves
xEst<-t(matrix(muXEst, nrow=J, ncol=n))+scTot%*%t(eigenTot[1:J+J,])
zEst<-t(matrix(muZEst, nrow=J, ncol=n))+scTot%*%t(eigenTot[1:J,])

xRec<-xEst
zRec<-zEst

#Recovery estimated oracle coefficients for LogiCenFD (M1)
meanSusEst<-(c(muZEstSus-muZEst, muXEstSus-muXEst)%*%eigenTot[,1:KSc])/J
meanContEst<-(c(muZEstCont-muZEst, muXEstCont-muXEst)%*%eigenTot[,1:KSc])/J

orcCoefEsts<-(meanSusEst-meanContEst)/scVals[1:KSc]

#Now re-estimating marginal conditional expectation of scores by taking mixture
scMeanSusZ<-as.vector((muZEstSus-muZEst)%*%eigenfuncsZ[,1:KZ])/J
scMeanContZ<-as.vector((muZEstCont-muZEst)%*%eigenfuncsZ[,1:KZ])/J
scVarSusZ<-t(eigenfuncsZ[,1:KZ])%*%sigmaSmthSusZ%*%eigenfuncsZ[,1:KZ]/J^2
scVarContZ<-t(eigenfuncsZ[,1:KZ])%*%sigmaSmthContZ%*%eigenfuncsZ[,1:KZ]/J^2

scZ<-estScoresZTest(a[,1:mnew], muZEst, eigenfuncsZ[,1:KZ], 
                    scMeanSusZ, scVarSusZ, scMeanContZ, scVarContZ, 1/2)

scMeanSusX<-as.vector((muXEstSus-muXEst)%*%eigenfuncsX[,1:KX])/J
scMeanContX<-as.vector((muXEstCont-muXEst)%*%eigenfuncsX[,1:KX])/J
scVarSusX<-t(eigenfuncsX[,1:KX])%*%sigmaSmthSusX%*%eigenfuncsX[,1:KX]/J^2
scVarContX<-t(eigenfuncsX[,1:KX])%*%sigmaSmthContX%*%eigenfuncsX[,1:KX]/J^2

scX<-estScoresXTest(v[,1:mnew], muXEst, eigenfuncsX[,1:KX], 
                    scMeanSusX, scVarSusX, scMeanContX, scVarContX, 1/2)

scTot<-cbind(scZ, scX)%*%scVecs
scTot<-scTot%>%apply(2, scale, center=F)

#Step 4: Fit models
trainFrame<-data.frame(treat=c(rep(1, n/2), rep(0, n/2)))
for (k in 1:KSc){
  trainFrame<-cbind(trainFrame, scTot[,k])
}
names(trainFrame)[2:(2 + KSc - 1)] <- c(paste0("sc", 1:KSc))
trainFrame$V<-v[,1:mnew]
trainFrame$A<-a[,1:mnew]

##3. Linear
gamForm <- formula(paste0("treat~", paste(paste0("sc", 1:KSc), collapse = "+")))
fitLinear<-glm(gamForm, data=trainFrame, family = binomial(link="logit"))

coefEsts<-as.vector(fitLinear$coefficients[1:KSc+1])

##4. GAM Uni
gamForm <- formula(paste0("treat~", paste(paste0("s(sc", 1:KSc, ", bs='cr', m=2)"), collapse = "+")))
fitGAMUni<-gam(gamForm, data=trainFrame, family=binomial(link="logit"))

##5. GAM Multi
gamForm <- formula(paste0("treat~s(", paste(paste0("sc", 1:KSc), collapse = ","), ", k=10, bs='tp')"))
fitGAMMulti<-gam(gamForm, data=trainFrame, family=binomial(link="logit"))

##Testing
xiSusTest<-matrix(1:3)%>%apply(1, function(k) rnorm(nTest/2, mean=meanSus[k], sd=sqrt(lambInit[k])))
xiContTest<-matrix(1:3)%>%apply(1, function(k) rnorm(nTest/2, mean=meanCont[k], sd=sqrt(lambInit[k])))

zSusTest<-t(matrix(muZ, nrow=J, ncol=nTest/2))+xiSusTest%*%t(fourier)
zContTest<-t(matrix(muZ, nrow=J, ncol=nTest/2))+xiContTest%*%t(fourier)

xSusTest<-t(matrix(muX, nrow=J, ncol=nTest/2))+xiSusTest%*%t(fourier)
xContTest<-t(matrix(muX, nrow=J, ncol=nTest/2))+xiContTest%*%t(fourier)

epsilonSusTest<-rnorm(nTest/2*J, sd=0.25)
epsilonContTest<-rnorm(nTest/2*J, sd=0.25)

aSusTest<-zSusTest+matrix(epsilonSusTest, nrow=nTest/2, ncol=J)
aContTest<-zContTest+matrix(epsilonContTest, nrow=nTest/2, ncol=J)
aTest<-rbind(aSusTest, aContTest)[,1:mnew]
aTest[aTest<0]<-0
aTest[aTest>1]<-1

vSusTest<-matrix(rbinom(nTest/2*J, size=1, prob=expit(as.vector(xSusTest))), nrow=nTest/2, ncol=J)
vContTest<-matrix(rbinom(nTest/2*J, size=1, prob=expit(as.vector(xContTest))), nrow=nTest/2, ncol=J)
vTest<-rbind(vSusTest, vContTest)[,1:mnew]
aTest[!vTest]<-NA

#Estimate scores
scZTest<-estScoresZTest(aTest, muZEst, eigenfuncsZ[,1:KZ], 
                        scMeanSusZ, scVarSusZ, scMeanContZ, scVarContZ, 1/2)
scXTest<-estScoresXTest(vTest, muXEst, eigenfuncsX[,1:KX], 
                        scMeanSusX, scVarSusX, scMeanContX, scVarContX, 1/2)

scTestTot<-cbind(scZTest, scXTest)%*%scVecs
scTestTot<-scTestTot%>%apply(2, scale, center=F)

##Prediction
susIndTest<-c(rep(1, nTest/2), rep(0, nTest/2))
testFrame<-data.frame(treat=c(rep(1, nTest/2), rep(0, nTest/2)))
for (k in 1:KSc){
  testFrame<-cbind(testFrame, scTestTot[,k])
}
names(testFrame)[2:(2 + KSc - 1)] <- c(paste0("sc", 1:KSc))
testFrame$V<-vTest
testFrame$A<-aTest

##1. Naive Random Forest
impDataNaive<-trainFrame%>%subset(select=c("treat", "A"))
impDataTestNaive<-testFrame%>%subset(select=c("A"))
pred<-fitRFNaive(impDataNaive, impDataTestNaive, mnew)
aucRFNaive<-roc.curve(scores.class0 = pred[testFrame$treat==1],
                           scores.class1 = pred[testFrame$treat==0],
                           curve=T)

##2. New Random Forest
fitRFNew<-randomForest(y=as.factor(trainFrame$treat), 
                       x=trainFrame%>%subset(select=paste0("sc", 1:KSc))%>%as.matrix())
pred<-as.vector(predict(fitRFNew, newdata = testFrame%>%subset(select=paste0("sc", 1:KSc)), type="prob")[,2])
aucRFNew<-roc.curve(scores.class0 = pred[testFrame$treat==1],
                         scores.class1 = pred[testFrame$treat==0],
                         curve=T)

##3. Linear
pred<-as.vector(predict(fitLinear, newdata = testFrame, type = "response"))
aucLinear<-roc.curve(scores.class0 = pred[testFrame$treat==1],
                          scores.class1 = pred[testFrame$treat==0],
                          curve=T)

##4. GAM Uni
pred<-as.vector(predict(fitGAMUni, newdata = testFrame, type = "response"))
aucGAMUni<-roc.curve(scores.class0 = pred[testFrame$treat==1],
                          scores.class1 = pred[testFrame$treat==0],
                          curve=T)

##5. GAM Multi
pred<-as.vector(predict(fitGAMMulti, newdata = testFrame, type = "response"))
aucGAMMulti<-roc.curve(scores.class0 = pred[testFrame$treat==1],
                            scores.class1 = pred[testFrame$treat==0],
                            curve=T)

saveRDS(list(muRecZ=muRecZ, muRecX=muRecX, psi1Rec=psi1Rec, psi2Rec=psi2Rec, 
             psi3Rec=psi3Rec, sc1Gen=sc1Gen, sc2Gen=sc2Gen, sc3Gen=sc3Gen,
             sc1Rec=sc1Rec, sc2Rec=sc2Rec, sc3Rec=sc3Rec, xRec=xRec, zRec=zRec, 
             aucRFNaive=aucRFNaive, aucRFNew=aucRFNew, aucLinear=aucLinear,
             aucGAMUni=aucGAMUni, aucGAMMulti=aucGAMMulti,
             orcCoefs=orcCoefs, orcCoefEsts=orcCoefEsts, coefEsts=coefEsts), 
        file=paste0("results", input, ".rds"))
