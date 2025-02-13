
# Load Libraries & Options
rm(list=ls())
library(OpenMx)
library(dplyr)
library(tidyverse)
source("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/02_scripts/miFunctions.R")

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data
vcu25w <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/01_tidy_data/vcu25w.csv")

# Select Variables for Analysis
vars <- c('rssrs','rsds') # list of variables names
nv <- 2 # number of variables
ntv <- nv*2 # number of total variables
selVars <- paste(vars,c(rep("_T1",nv),rep("_T2",nv)),sep="")

mzfData <- subset(vcu25w, (zygo == "MZ" & sex_T1 == 2 & sex_T2 == 2), select = selVars)
dzfData <- subset(vcu25w, (zygo == "DZ" & sex_T1 == 2 & sex_T2 == 2), select = selVars)
mzmData <- subset(vcu25w, (zygo == "MZ" & sex_T1 == 1 & sex_T2 == 1), select = selVars)
dzmData <- subset(vcu25w, (zygo == "DZ" & sex_T1 == 1 & sex_T2 == 1), select = selVars)
dzoData <- subset(vcu25w, (zygo == "DZ" & sex_T1 != sex_T2), select = selVars)

colMeans(mzfData,na.rm=TRUE)
colMeans(dzfData,na.rm=TRUE)
colMeans(mzmData,na.rm=TRUE)
colMeans(dzmData,na.rm=TRUE)
colMeans(dzoData,na.rm=TRUE)

cov(mzfData,use="complete")
cov(dzfData,use="complete")
cov(mzmData,use="complete")
cov(dzmData,use="complete")
cov(dzoData,use="complete")

cor(mzfData,use="complete")
cor(dzfData,use="complete")
cor(mzmData,use="complete")
cor(dzmData,use="complete")
cor(dzoData,use="complete")

# Set Starting Values
svMe <- c(0.01,0.01) # start value for means
svVa <- c(.5,.8) # start value for variances
lbVa <- .0001 # lower bound for variances

# Create Labels
labMeMZf <- labVars("meanMZf",selVars)
labMeDZf <- labVars("meanDZf",selVars)
labMeMZm <- labVars("meanMZm",selVars)
labMeDZm <- labVars("meanDZm",selVars)
labMeDZo <- labVars("meanDZo",selVars)
labMeZf <- labVars("meanZf",selVars) # Added this 1/29 to try to run EMVP
labMeZm <- labVars("meanZm",selVars) # Added this 1/29 to try to run EMVP
labMeZ <- labVars("meanZ",selVars)
labCvMZf <- labLower("covMZf",ntv)
labCvDZf <- labLower("covDZf",ntv)
labCvMZm <- labLower("covMZm",ntv)
labCvDZm <- labLower("covDZm",ntv)
labCvDZo <- labLower("covDZo",ntv)
labCvZf <- labLower("covZf",ntv) # Added this 1/29 to try to run EMVP
labCvZm <- labLower("covZm",ntv) # Added this 1/29 to try to run EMVP
labCvZ <- labLower("covZ",ntv)
labVaMZf <- labDiag("covMZf",ntv)
labVaDZf <- labDiag("covDZf",ntv)
labVaMZm <- labDiag("covMZm",ntv)
labVaDZm <- labDiag("covDZm",ntv)
labVaDZo <- labDiag("covDZo",ntv)
labVaZf <- labDiag("covZf",ntv) # Added this 1/29 to try to run EMVP
labVaZm <- labDiag("covZm",ntv) # Added this 1/29 to try to run EMVP
labVaZ <- labDiag("covZ",ntv)

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE MODEL

# Create Algebra for expected Mean Matrices
meanMZf <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeMZf, name="meanMZf" )
meanDZf <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeDZf, name="meanDZf" )
meanMZm <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeMZm, name="meanMZm" )
meanDZm <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeDZm, name="meanDZm" )
meanDZo <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labMeDZo, name="meanDZo" )

# Create Algebra for expected Variance/Covariance Matrices
covMZf <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=labCvMZf, name="covMZf" )
covDZf <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=labCvDZf, name="covDZf" )
covMZm <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=labCvMZm, name="covMZm" )
covDZm <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=labCvDZm, name="covDZm" )
covDZo <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), lbound=valDiag(lbVa,ntv),
                    labels=labCvDZo, name="covDZo" )

# Create Data Objects for Multiple Groups
dataMZf <- mxData( observed=mzfData, type="raw" )
dataDZf <- mxData( observed=dzfData, type="raw" )
dataMZm <- mxData( observed=mzmData, type="raw" )
dataDZm <- mxData( observed=dzmData, type="raw" )
dataDZo <- mxData( observed=dzoData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZf <- mxExpectationNormal( covariance="covMZf", means="meanMZf", dimnames=selVars )
expDZf <- mxExpectationNormal( covariance="covDZf", means="meanDZf", dimnames=selVars )
expMZm <- mxExpectationNormal( covariance="covMZm", means="meanMZm", dimnames=selVars )
expDZm <- mxExpectationNormal( covariance="covDZm", means="meanDZm", dimnames=selVars )
expDZo <- mxExpectationNormal( covariance="covDZo", means="meanDZo", dimnames=selVars )
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
modelMZf <- mxModel( meanMZf, covMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZf <- mxModel( meanDZf, covDZf, dataDZf, expDZf, funML, name="DZf" )
modelMZm <- mxModel( meanMZm, covMZm, dataMZm, expMZm, funML, name="MZm" )
modelDZm <- mxModel( meanDZm, covDZm, dataDZm, expDZm, funML, name="DZm" )
modelDZo <- mxModel( meanDZo, covDZo, dataDZo, expDZo, funML, name="DZo" )
multi <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

# Create Confidence Interval Objects
ciCov <- mxCI( c('MZf.covMZf','DZf.covDZf','MZm.covMZm','DZm.covDZm','DZo.covDZo') )
ciMean <- mxCI( c('MZf.meanMZf','DZf.meanDZf','MZm.meanMZm','DZm.meanDZm','DZo.meanDZo') )

#to run on Neda's laptop
#mxOption(NULL,"Default optimizer","SLSQP")

# Build Saturated Model with Confidence Intervals
modelSAT <- mxModel( "twoSAT", modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, ciCov, ciMean )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run Saturated Model
fitSAT <- mxRun( modelSAT, intervals=F )
sumSAT <- summary( fitSAT )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitSAT)
fitEsts(fitSAT)
mxGetExpected(fitSAT, "means")
mxGetExpected(fitSAT, "covariance")

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Constrain expected Means to be equal across Twin Order
modelEMO <- mxModel( fitSAT, name="twoEMO" )
for (i in 1:nv) {
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeMZf[nv+i],labMeMZf[i]), free=TRUE, values=svMe, newlabels=labMeMZf[i] )
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeDZf[nv+i],labMeDZf[i]), free=TRUE, values=svMe, newlabels=labMeDZf[i] )
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeMZm[nv+i],labMeMZm[i]), free=TRUE, values=svMe, newlabels=labMeMZm[i] )
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeDZm[nv+i],labMeDZm[i]), free=TRUE, values=svMe, newlabels=labMeDZm[i] )
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeDZo[nv+i],labMeDZo[i]), free=TRUE, values=svMe, newlabels=labMeDZo[i] )}
fitEMO <- mxRun( modelEMO, intervals=F )
fitGofs(fitEMO); fitEsts(fitEMO)

# Constrain expected Means and Variances to be equal across Twin Order
modelEMVO <- mxModel( fitEMO, name="twoEMVO" )
for (i in 1:nv) {
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaMZf[nv+i],labVaMZf[i]), free=TRUE, values=c(.75,.75), newlabels=labVaMZf[i] )
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaDZf[nv+i],labVaDZf[i]), free=TRUE, values=c(.72,.72), newlabels=labVaDZf[i] )
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaMZm[nv+i],labVaMZm[i]), free=TRUE, values=c(.87,.87), newlabels=labVaMZm[i] )
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaDZm[nv+i],labVaDZm[i]), free=TRUE, values=c(0.88,0.88), newlabels=labVaDZm[i] )
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaDZo[nv+i],labVaDZo[i]), free=TRUE, values=c(0.83,0.83), newlabels=labVaDZo[i] )}
fitEMVO <- mxRun( modelEMVO, intervals=F )
fitGofs(fitEMVO); fitEsts(fitEMVO)

# Constrain expected Means and Variances to be equal across Twin Order and Zygosity
modelEMVZ <- mxModel( fitEMO, name="twoEMVZ" )
for (i in 1:nv) {
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labMeMZf[i],labMeDZf[i]), free=TRUE, values=svMe, newlabels=labMeZf[i] ) # Added f to this to try to run EMVP
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labVaMZf[i],labVaDZf[i]), free=TRUE, values=c(.81,.81), newlabels=labVaZf[i] )
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labMeMZm[i],labMeDZm[i]), free=TRUE, values=svMe, newlabels=labMeZm[i] ) # Added m to this to try to run EMVP
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labVaMZm[i],labVaDZm[i]), free=TRUE, values=c(.81,.81), newlabels=labVaZm[i] )}
fitEMVZ <- mxRun( modelEMVZ, intervals=F )
fitGofs(fitEMVZ); fitEsts(fitEMVZ)

# Constrain expected Means and Variances to be equal across twin order and zygosity and SS/OS
modelEMVP <- mxModel( fitEMVZ, name="twoEMVP" ) # 1/29 made strict=FALSE, not sure what this means
for (i in 1:nv) {
  modelEMVP <- omxSetParameters( modelEMVP, label=c(labMeZf[i],labMeDZo[i]), free=TRUE, strict=FALSE, values=svMe, newlabels=labMeZf[i] )
  modelEMVP <- omxSetParameters( modelEMVP, label=c(labVaZf[i],labVaDZo[i]), free=TRUE, strict=FALSE, values=c(.77,.77), newlabels=labVaZf[i] )
  modelEMVP <- omxSetParameters( modelEMVP, label=c(labMeZm[i],labMeDZo[i]), free=TRUE, strict=FALSE, values=svMe, newlabels=labMeZm[i] )
  modelEMVP <- omxSetParameters( modelEMVP, label=c(labVaZm[i],labVaDZo[i]), free=TRUE, strict=FALSE, values=c(.77,.77), newlabels=labVaZm[i] )}
fitEMVP <- mxRun( modelEMVP, intervals=F )
fitGofs(fitEMVP); fitEsts(fitEMVP)

# Constrain expected Means and Variances to be equal across twin order and zygosity and SS/OS and sex
modelEMVS <- mxModel( fitEMVP, name="twoEMVS" )
for (i in 1:nv) {
  modelEMVS <- omxSetParameters( modelEMVS, label=c(labMeZf[i],labMeZm[i]), free=TRUE, values=svMe, newlabels=labMeZ[i] )
  modelEMVS <- omxSetParameters( modelEMVS, label=c(labVaZf[i],labVaZm[i]), free=TRUE, values=c(.77,.77), newlabels=labVaZ[i] )}
fitEMVS <- mxRun( modelEMVS, intervals=F )
fitGofs(fitEMVS); fitEsts(fitEMVS)

# Print Comparative Fit Statistics
mxCompare( fitSAT, subs <- list(fitEMO, fitEMVO, fitEMVZ, fitEMVP, fitEMVS) )

# ----------------------------------------------------------------------------------------------------------------------
# ACEra MODEL STARTS HERE

# Set Starting Values
svMe      <- c(0.01, 0.01)                     # start value for means
svPa      <- .4                    # start value for variance component A
svPe      <- .5                   # start value for variance component E

# Create Algebra for expected Mean Matrices
meanGf      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSf","meanSDSf", "meanSRSf","meanSDSf"), name="meanGf" )
meanGm      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSm","meanSDSm", "meanSRSm","meanSDSm"), name="meanGm" )
meanGo      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSo","meanSDSo", "meanSRSo","meanSDSo"), name="meanGo" )

# Create Matrices for Variance Components
covAf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAf",nv), name="VAf" ) # lbound=valDiag(1e-4,ntv), 
covCf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0.1,nv),
                        label=labLower("VCf",nv), name="VCf" )
covEf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv), 
                        label=labLower("VEf",nv), name="VEf" )
covAm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAm",nv), name="VAm" )
covCm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0.1,nv), 
                        label=labLower("VCm",nv), name="VCm" )
covEm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv), 
                        label=labLower("VEm",nv), name="VEm" )
covAms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0,nv),
                        label=labLower("VAms",nv), name="VAms" )
covCms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
                         label=labLower("VCms",nv), name="VCms" )

# Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
zeroMat <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values = 0, name = "zeroMat")
signA <- mxAlgebra( (-1)^omxLessThan(VAf, zeroMat) * (-1)^omxLessThan(VAm, zeroMat), name="signA") # 2/4: had to change this code to make a zero matrix. added zeroMat to parsZo
signC <- mxAlgebra( (-1)^omxLessThan(VCf, zeroMat) * (-1)^omxLessThan(VCm, zeroMat), name="signC")
# Calculate absolute covariation between Males and Females and then un-absoluting the product
covAos <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name="VAos")
covCos <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm)))), name="VCos")
# Calculate rg/rc from reparameterized model
pathRg <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
pathRc <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm))))/sqrt(VCf*(VCm+VCms)), name="rc")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covPf <- mxAlgebra( expression= VAf+VCf+VEf, name="Vf" )
covPm <- mxAlgebra( expression= VAm+VCm+VEm+VAms+VCms, name="Vm" )
covMZf <- mxAlgebra( expression= VAf+VCf, name="cMZf" )
covDZf <- mxAlgebra( expression= 0.5%x%VAf+ VCf, name="cDZf" )
covMZm <- mxAlgebra( expression= VAm+VCm+VAms+VCms, name="cMZm" )
covDZm <- mxAlgebra( expression= 0.5%x%VAm+ VCm+ 0.5%x%VAms+VCms, name="cDZm" )
covDZo <- mxAlgebra( expression= 0.5%x%VAos+VCos, name="cDZo" )
expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
expCovDZo <- mxAlgebra( expression= rbind( cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name="expCovDZo" )

# Create Data Objects for Multiple Groups
dataMZf <- mxData( observed=mzfData, type="raw" )
dataDZf <- mxData( observed=dzfData, type="raw" )
dataMZm <- mxData( observed=mzmData, type="raw" )
dataDZm <- mxData( observed=dzmData, type="raw" )
dataDZo <- mxData( observed=dzoData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZf <- mxExpectationNormal( covariance="expCovMZf", means="meanGf", dimnames=selVars )
expDZf <- mxExpectationNormal( covariance="expCovDZf", means="meanGf", dimnames=selVars )
expMZm <- mxExpectationNormal( covariance="expCovMZm", means="meanGm", dimnames=selVars )
expDZm <- mxExpectationNormal( covariance="expCovDZm", means="meanGm", dimnames=selVars )
expDZo <- mxExpectationNormal( covariance="expCovDZo", means="meanGo", dimnames=selVars )
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
parsZf <- list( covAf, covCf, covEf, covPf )
parsZm <- list( covAm, covCm, covEm, covPm, covAms, covCms )
parsZo <- list( parsZm, parsZf, zeroMat, signA, signC, covAos, covCos, pathRg, pathRc )
modelMZf <- mxModel( parsZf, meanGf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZf <- mxModel( parsZf, meanGf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf" )
modelMZm <- mxModel( parsZm, meanGm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm" )
modelDZm <- mxModel( parsZm, meanGm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm" )
modelDZo <- mxModel( parsZo, meanGo, covDZo, expCovDZo, dataDZo, expDZo, funML, name="DZo" )
multi <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

## Create Algebra for Standardization
matIf      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="If" )
matIm      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="Im" )

## Calculate genetic,  environmental, and phenotypic correlations
corAf      <- mxAlgebra( expression=solve(sqrt(If*VAf))%&%VAf, name ="rAf" )
corCf      <- mxAlgebra( expression=solve(sqrt(If*VCf))%&%VCf, name ="rCf" )
corEf      <- mxAlgebra( expression=solve(sqrt(If*VEf))%&%VEf, name ="rEf" )
corPf      <- mxAlgebra( expression=solve(sqrt(If*Vf))%&%Vf, name ="rPf")

corAm      <- mxAlgebra( expression=solve(sqrt(Im*VAm))%&%VAm, name ="rAm" )
corCm      <- mxAlgebra( expression=solve(sqrt(Im*VCm))%&%VCm, name ="rCm" )
corEm      <- mxAlgebra( expression=solve(sqrt(Im*VEm))%&%VEm, name ="rEm" )
corPm      <- mxAlgebra( expression=solve(sqrt(Im*Vm))%&%Vm, name ="rPm")

## Create Algebra for Variance Components ( formatting output)
#Unstandardized variance components
rowUV     <- rep('UV',nv)
colUV     <- rep(c('VA','VC','VE'),each=nv)
estUVf     <- mxAlgebra( expression=cbind(VAf,VCf,VEf), name="UVf", dimnames=list(rowUV,colUV) )
estUVm     <- mxAlgebra( expression=cbind(VAm,VCm,VEm), name="UVm", dimnames=list(rowUV,colUV) )

#Standardized variance components
rowSV     <- rep('SV',nv)
colSV     <- rep(c('SA','SC','SE'),each=nv)
estSVf     <- mxAlgebra( expression=cbind(VAf/Vf,VCf/Vf,VEf/Vf), name="SVf", dimnames=list(rowSV,colSV) )
estSVm     <- mxAlgebra( expression=cbind(VAm/Vm,VCm/Vm,VEm/Vm), name="SVm", dimnames=list(rowSV,colSV) )

#Correlations
rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rC','rE','rP'),each=nv)
estcorrf     <- mxAlgebra( expression=cbind(rAf,rCf,rEf,rPf), name="corrf", dimnames=list(rowcorr,colcorr) )
estcorrm     <- mxAlgebra( expression=cbind(rAm,rCm,rEm,rPm), name="corrm", dimnames=list(rowcorr,colcorr) )

## Create Confidence Interval Objects
ciUV      <- mxCI( c("UV") )
ciSV      <- mxCI( c("SVf", "SVm") )
cicorr      <- mxCI( c("corrf", "corrm") )

## Build Model with Confidence Intervals
calc      <- list( matIf,  corAf, corCf, corEf,corPf, estUVf, estSVf,estcorrf,
                   matIm,  corAm, corCm, corEm,corPm, estUVm, estSVm,estcorrm, ciSV, cicorr)
modelACEra  <- mxModel( "twoACEra5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACEra Model - Qualitative (Ra) & Quantative Sex Differences ACE model
fitACEra <- mxRun( modelACEra, intervals=T )
sumACEra <- summary( fitACEra )

# Compare with Saturated Model
mxCompare( fitSAT, fitACEra )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitACEra)
fitEsts(fitACEra)
fitEstCis(fitACEra)

# Print Covariance & Correlation Matrices
fitACEra$algebras$UVf
fitACEra$algebras$UVm
fitACEra$algebras$SVf
fitACEra$algebras$SVm
fitACEra$algebras$corrf
fitACEra$algebras$corrm

# fitACEra2 <- mxRun( fitACEra, intervals=T )
# fitACEra3 <- mxRun( fitACEra2, intervals=T )
# fitACEra4 <- mxRun( fitACEra3, intervals=T )
# fitACEra5 <- mxRun( fitACEra4, intervals=T )
# fitACEra6 <- mxRun( fitACEra5, intervals=T )
# fitACEra7 <- mxRun( fitACEra6, intervals=T )
# fitACEra8 <- mxRun( fitACEra7, intervals=T )
# fitACEra9 <- mxRun( fitACEra8, intervals=T )
# fitACEra10 <- mxRun( fitACEra9, intervals=T )
# fitACEra11 <- mxRun( fitACEra10, intervals=T )
# fitACEra12 <- mxRun( fitACEra11, intervals=T )
# fitACEra13 <- mxRun( fitACEra12, intervals=T )
# fitACEra14 <- mxRun( fitACEra13, intervals=T )
# fitACEra15 <- mxRun( fitACEra14, intervals=T )
# 
# 
# fitGofs(fitACEra10)
# fitEsts(fitACEra10)
# fitEstCis(fitACEra10)
# 
# fitACEra10$algebras$UVf
# fitACEra10$algebras$UVm
# fitACEra10$algebras$SVf
# fitACEra10$algebras$SVm
# fitACEra10$algebras$corrf
# fitACEra10$algebras$corrm

# ----------------------------------------------------------------------------------------------------------------------
# ACErc MODEL HERE

# Set Starting Values 
svMe      <- c(0.01, 0.01)                     # start value for means
svPa      <- .4                    # start value for variance component A
svPe      <- .5                   # start value for variance component E

# Create Algebra for expected Mean Matrices
meanGf      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, 
                         labels=c("meanSRSf","meanSDSf", "meanSRSf","meanSDSf"), name="meanGf" )
meanGm      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, 
                         labels=c("meanSRSm","meanSDSm", "meanSRSm","meanSDSm"), name="meanGm" )
meanGo      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, 
                         labels=c("meanSRSo","meanSDSo", "meanSRSo","meanSDSo"), name="meanGo" )

# Create Matrices for Variance Components
covAf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAf",nv), name="VAf" ) # lbound=valDiag(1e-4,ntv),  
covCf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0.1,nv),
                        label=labLower("VCf",nv), name="VCf" )
covEf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                        label=labLower("VEf",nv), name="VEf" )
covAm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAm",nv), name="VAm" ) 
covCm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0.1,nv),
                        label=labLower("VCm",nv), name="VCm" )
covEm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                        label=labLower("VEm",nv), name="VEm" )
covAms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
                         label=labLower("VAms",nv), name="VAms" )
covCms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0.1,nv),
                         label=labLower("VCms",nv), name="VCms" )

# Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
# signA <- mxAlgebra( ((-1)^omxLessThan(VAf,0))*((-1)^omxLessThan(VAm,0)), name="signA") # original code, not working
# signC <- mxAlgebra( ((-1)^omxLessThan(VCf,0))*((-1)^omxLessThan(VCm,0)), name="signC")
zeroMat <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values = 0, name = "zeroMat")
signA <- mxAlgebra( (-1)^omxLessThan(VAf, zeroMat) * (-1)^omxLessThan(VAm, zeroMat), name="signA") # 2/4: had to change this code to make a zero matrix. added zeroMat to parsZo
signC <- mxAlgebra( (-1)^omxLessThan(VCf, zeroMat) * (-1)^omxLessThan(VCm, zeroMat), name="signC")
# Calculate absolute covariation between Males and Females and then un-absoluting the product
covAos <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name="VAos")
covCos <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm)))), name="VCos")
# Calculate rg/rc from reparameterized model
pathRg <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
pathRc <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm))))/sqrt(VCf*(VCm+VCms)), name="rc")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covPf <- mxAlgebra( expression= VAf+VCf+VEf, name="Vf" )
covPm <- mxAlgebra( expression= VAm+VCm+VEm+VAms+VCms, name="Vm" )
covMZf <- mxAlgebra( expression= VAf+VCf, name="cMZf" )
covDZf <- mxAlgebra( expression= 0.5%x%VAf+ VCf, name="cDZf" )
covMZm <- mxAlgebra( expression= VAm+VCm+VAms+VCms, name="cMZm" )
covDZm <- mxAlgebra( expression= 0.5%x%VAm+ VCm+ 0.5%x%VAms+VCms, name="cDZm" )
covDZo <- mxAlgebra( expression= 0.5%x%VAos+VCos, name="cDZo" )
expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
expCovDZo <- mxAlgebra( expression= rbind( cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name="expCovDZo" )

# Create Data Objects for Multiple Groups
dataMZf <- mxData( observed=mzfData, type="raw" )
dataDZf <- mxData( observed=dzfData, type="raw" )
dataMZm <- mxData( observed=mzmData, type="raw" )
dataDZm <- mxData( observed=dzmData, type="raw" )
dataDZo <- mxData( observed=dzoData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZf <- mxExpectationNormal( covariance="expCovMZf", means="meanGf", dimnames=selVars )
expDZf <- mxExpectationNormal( covariance="expCovDZf", means="meanGf", dimnames=selVars )
expMZm <- mxExpectationNormal( covariance="expCovMZm", means="meanGm", dimnames=selVars )
expDZm <- mxExpectationNormal( covariance="expCovDZm", means="meanGm", dimnames=selVars )
expDZo <- mxExpectationNormal( covariance="expCovDZo", means="meanGo", dimnames=selVars )
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
parsZf <- list( covAf, covCf, covEf, covPf )
parsZm <- list( covAm, covCm, covEm, covPm, covAms, covCms )
parsZo <- list( parsZm, parsZf, zeroMat, signA, signC, covAos, covCos, pathRg, pathRc )
modelMZf <- mxModel( parsZf, meanGf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZf <- mxModel( parsZf, meanGf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf" )
modelMZm <- mxModel( parsZm, meanGm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm" )
modelDZm <- mxModel( parsZm, meanGm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm" )
modelDZo <- mxModel( parsZo, meanGo, covDZo, expCovDZo, dataDZo, expDZo, funML, name="DZo" )
multi <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

## Create Algebra for Standardization 
matIf      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="If" )
matIm      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="Im" )

## Calculate genetic,  environmental, and phenotypic correlations
corAf      <- mxAlgebra( expression=solve(sqrt(If*VAf))%&%VAf, name ="rAf" ) 
corCf      <- mxAlgebra( expression=solve(sqrt(If*VCf))%&%VCf, name ="rCf" )
corEf      <- mxAlgebra( expression=solve(sqrt(If*VEf))%&%VEf, name ="rEf" )
corPf      <- mxAlgebra( expression=solve(sqrt(If*Vf))%&%Vf, name ="rPf") 

corAm      <- mxAlgebra( expression=solve(sqrt(Im*VAm))%&%VAm, name ="rAm" ) 
corCm      <- mxAlgebra( expression=solve(sqrt(Im*VCm))%&%VCm, name ="rCm" )
corEm      <- mxAlgebra( expression=solve(sqrt(Im*VEm))%&%VEm, name ="rEm" )
corPm      <- mxAlgebra( expression=solve(sqrt(Im*Vm))%&%Vm, name ="rPm") 

## Create Algebra for Variance Components ( formatting output)
#Unstandardized variance components
rowUV     <- rep('UV',nv)
colUV     <- rep(c('VA','VC','VE'),each=nv)
estUVf     <- mxAlgebra( expression=cbind(VAf,VCf,VEf), name="UVf", dimnames=list(rowUV,colUV) )
estUVm     <- mxAlgebra( expression=cbind(VAm,VCm,VEm), name="UVm", dimnames=list(rowUV,colUV) )

#Standardized variance components
rowSV     <- rep('SV',nv)
colSV     <- rep(c('SA','SC','SE'),each=nv)
estSVf     <- mxAlgebra( expression=cbind(VAf/Vf,VCf/Vf,VEf/Vf), name="SVf", dimnames=list(rowSV,colSV) )
estSVm     <- mxAlgebra( expression=cbind(VAm/Vm,VCm/Vm,VEm/Vm), name="SVm", dimnames=list(rowSV,colSV) )

#Correlations
rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rC','rE','rP'),each=nv)
estcorrf     <- mxAlgebra( expression=cbind(rAf,rCf,rEf,rPf), name="corrf", dimnames=list(rowcorr,colcorr) )
estcorrm     <- mxAlgebra( expression=cbind(rAm,rCm,rEm,rPm), name="corrm", dimnames=list(rowcorr,colcorr) )

## Create Confidence Interval Objects
ciUV      <- mxCI( c("UV") )
ciSV      <- mxCI( c("SVf", "SVm") )
cicorr      <- mxCI( c("corrf", "corrm") )

## Build Model with Confidence Intervals
calc      <- list( matIf,  corAf, corCf, corEf,corPf, estUVf, estSVf,estcorrf,
                   matIm,  corAm, corCm, corEm,corPm, estUVm, estSVm,estcorrm, ciSV, cicorr)
modelACErc  <- mxModel( "twoACErc5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACErc Model - Qualitative (Rc) & Quantative Sex Differences ACE model
fitACErc <- mxRun( modelACErc, intervals=T )
sumACErc <- summary( fitACErc )

# Compare with Saturated Model
mxCompare( fitSAT, fitACErc )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitACErc)
fitEsts(fitACErc)
fitEstCis(fitACErc)

# Print Covariance & Correlation Matrices
fitACErc$algebras$UVf
fitACErc$algebras$UVm
fitACErc$algebras$SVf
fitACErc$algebras$SVm
fitACErc$algebras$corrf
fitACErc$algebras$corrm

# ----------------------------------------------------------------------------------------------------------------------
# ACEq MODEL STARTS HERE

# Set Starting Values
svMe      <- c(0.01, 0.01)                     # start value for means
svPa      <- .4                    # start value for variance component A
svPe      <- .5                   # start value for variance component E

# Create Algebra for expected Mean Matrices
meanGf      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSf","meanSDSf", "meanSRSf","meanSDSf"), name="meanGf" )
meanGm      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSm","meanSDSm", "meanSRSm","meanSDSm"), name="meanGm" )
meanGo      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSo","meanSDSo", "meanSRSo","meanSDSo"), name="meanGo" )

# Create Matrices for Variance Components
covAf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAf",nv), name="VAf" ) # lbound=valDiag(1e-4,ntv), 
covCf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0.1,nv),
                        label=labLower("VCf",nv), name="VCf" )
covEf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                        label=labLower("VEf",nv), name="VEf" )
covAm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAm",nv), name="VAm" )
covCm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(0.1,nv),
                        label=labLower("VCm",nv), name="VCm" )
covEm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                        label=labLower("VEm",nv), name="VEm" )
covAms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
                        label=labLower("VAms",nv), name="VAms" ) 
covCms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
                         label=labLower("VCms",nv), name="VCms" )

# Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
# signA <- mxAlgebra( ((-1)^omxLessThan(VAf,0))*((-1)^omxLessThan(VAm,0)), name="signA") # original code, not working
# signC <- mxAlgebra( ((-1)^omxLessThan(VCf,0))*((-1)^omxLessThan(VCm,0)), name="signC")
zeroMat <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values = 0, name = "zeroMat")
signA <- mxAlgebra( (-1)^omxLessThan(VAf, zeroMat) * (-1)^omxLessThan(VAm, zeroMat), name="signA") # 2/4: had to change this code to make a zero matrix. added zeroMat to parsZo
signC <- mxAlgebra( (-1)^omxLessThan(VCf, zeroMat) * (-1)^omxLessThan(VCm, zeroMat), name="signC")
# Calculate absolute covariation between Males and Females and then un-absoluting the product
covAos <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name="VAos")
covCos <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm)))), name="VCos")
# Calculate rg/rc from reparameterized model
pathRg <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
pathRc <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm))))/sqrt(VCf*(VCm+VCms)), name="rc")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covPf <- mxAlgebra( expression= VAf+VCf+VEf, name="Vf" )
covPm <- mxAlgebra( expression= VAm+VCm+VEm+VAms+VCms, name="Vm" )
covMZf <- mxAlgebra( expression= VAf+VCf, name="cMZf" )
covDZf <- mxAlgebra( expression= 0.5%x%VAf+ VCf, name="cDZf" )
covMZm <- mxAlgebra( expression= VAm+VCm+VAms+VCms, name="cMZm" )
covDZm <- mxAlgebra( expression= 0.5%x%VAm+ VCm+ 0.5%x%VAms+VCms, name="cDZm" )
covDZo <- mxAlgebra( expression= 0.5%x%VAos+VCos, name="cDZo" )
expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
expCovDZo <- mxAlgebra( expression= rbind( cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name="expCovDZo" )

# Create Data Objects for Multiple Groups
dataMZf <- mxData( observed=mzfData, type="raw" )
dataDZf <- mxData( observed=dzfData, type="raw" )
dataMZm <- mxData( observed=mzmData, type="raw" )
dataDZm <- mxData( observed=dzmData, type="raw" )
dataDZo <- mxData( observed=dzoData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZf <- mxExpectationNormal( covariance="expCovMZf", means="meanGf", dimnames=selVars )
expDZf <- mxExpectationNormal( covariance="expCovDZf", means="meanGf", dimnames=selVars )
expMZm <- mxExpectationNormal( covariance="expCovMZm", means="meanGm", dimnames=selVars )
expDZm <- mxExpectationNormal( covariance="expCovDZm", means="meanGm", dimnames=selVars )
expDZo <- mxExpectationNormal( covariance="expCovDZo", means="meanGo", dimnames=selVars )
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
parsZf <- list( covAf, covCf, covEf, covPf )
parsZm <- list( covAm, covCm, covEm, covPm, covAms, covCms )
parsZo <- list( parsZm, parsZf, zeroMat, signA, signC, covAos, covCos, pathRg, pathRc )
modelMZf <- mxModel( parsZf, meanGf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZf <- mxModel( parsZf, meanGf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf" )
modelMZm <- mxModel( parsZm, meanGm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm" )
modelDZm <- mxModel( parsZm, meanGm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm" )
modelDZo <- mxModel( parsZo, meanGo, covDZo, expCovDZo, dataDZo, expDZo, funML, name="DZo" )
multi <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

## Create Algebra for Standardization
matIf      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="If" )
matIm      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="Im" )

## Calculate genetic,  environmental, and phenotypic correlations
corAf      <- mxAlgebra( expression=solve(sqrt(If*VAf))%&%VAf, name ="rAf" )
corCf      <- mxAlgebra( expression=solve(sqrt(If*VCf))%&%VCf, name ="rCf" )
corEf      <- mxAlgebra( expression=solve(sqrt(If*VEf))%&%VEf, name ="rEf" )
corPf      <- mxAlgebra( expression=solve(sqrt(If*Vf))%&%Vf, name ="rPf")

corAm      <- mxAlgebra( expression=solve(sqrt(Im*VAm))%&%VAm, name ="rAm" )
corCm      <- mxAlgebra( expression=solve(sqrt(Im*VCm))%&%VCm, name ="rCm" )
corEm      <- mxAlgebra( expression=solve(sqrt(Im*VEm))%&%VEm, name ="rEm" )
corPm      <- mxAlgebra( expression=solve(sqrt(Im*Vm))%&%Vm, name ="rPm")

## Create Algebra for Variance Components ( formatting output)
#Unstandardized variance components
rowUV     <- rep('UV',nv)
colUV     <- rep(c('VA','VC','VE'),each=nv)
estUVf     <- mxAlgebra( expression=cbind(VAf,VCf,VEf), name="UVf", dimnames=list(rowUV,colUV) )
estUVm     <- mxAlgebra( expression=cbind(VAm,VCm,VEm), name="UVm", dimnames=list(rowUV,colUV) )

#Standardized variance components
rowSV     <- rep('SV',nv)
colSV     <- rep(c('SA','SC','SE'),each=nv)
estSVf     <- mxAlgebra( expression=cbind(VAf/Vf,VCf/Vf,VEf/Vf), name="SVf", dimnames=list(rowSV,colSV) )
estSVm     <- mxAlgebra( expression=cbind(VAm/Vm,VCm/Vm,VEm/Vm), name="SVm", dimnames=list(rowSV,colSV) )

#Correlations
rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rC','rE','rP'),each=nv)
estcorrf     <- mxAlgebra( expression=cbind(rAf,rCf,rEf,rPf), name="corrf", dimnames=list(rowcorr,colcorr) )
estcorrm     <- mxAlgebra( expression=cbind(rAm,rCm,rEm,rPm), name="corrm", dimnames=list(rowcorr,colcorr) )

## Create Confidence Interval Objects
ciUV      <- mxCI( c("UV") )
ciSV      <- mxCI( c("SVf", "SVm") )
cicorr      <- mxCI( c("corrf", "corrm") )

## Build Model with Confidence Intervals
calc      <- list( matIf,  corAf, corCf, corEf,corPf, estUVf, estSVf,estcorrf,
                   matIm,  corAm, corCm, corEm,corPm, estUVm, estSVm,estcorrm, ciSV, cicorr)
modelACEq  <- mxModel( "twoACEq5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACEq model - Quantitative non-scalar Sex Differences ACE model
fitACEq <- mxRun( modelACEq, intervals=T )
sumACEq <- summary( fitACEq )

# Compare with Saturated Model
mxCompare( fitSAT, fitACEq )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitACEq)
fitEsts(fitACEq)
fitEstCis(fitACEq)

# Print Covariance & Correlation Matrices
fitACEq$algebras$UVf
fitACEq$algebras$UVm
fitACEq$algebras$SVf
fitACEq$algebras$SVm
fitACEq$algebras$corrf
fitACEq$algebras$corrm

# ----------------------------------------------------------------------------------------------------------------------
# ACE (no sex diff) MODEL STARTS HERE

# Set Starting Values
svMe      <- c(0.01, 0.01)                     # start value for means
svPa      <- .4                    # start value for variance component A
svPe      <- .5                   # start value for variance component E

# Create Algebra for expected Mean Matrices
meanGf      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSf","meanSDSf", "meanSRSf","meanSDSf"), name="meanGf" )
meanGm      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSm","meanSDSm", "meanSRSm","meanSDSm"), name="meanGm" )
meanGo      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSo","meanSDSo", "meanSRSo","meanSDSo"), name="meanGo" )

# Create Matrices for Variance Components
covAf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAf",nv), name="VAf" ) # lbound=valDiag(1e-4,ntv), 
covCf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VCf",nv), name="VCf" )
covEf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                        label=labLower("VEf",nv), name="VEf" )
covAm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAm",nv), name="VAm" )
covCm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VCm",nv), name="VCm" )
covEm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                        label=labLower("VEm",nv), name="VEm" )
covAms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
                         label=labLower("VAms",nv), name="VAms" ) 
covCms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
                         label=labLower("VCms",nv), name="VCms" )

# Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
# signA <- mxAlgebra( ((-1)^omxLessThan(VAf,0))*((-1)^omxLessThan(VAm,0)), name="signA") # original code, not working
# signC <- mxAlgebra( ((-1)^omxLessThan(VCf,0))*((-1)^omxLessThan(VCm,0)), name="signC")
zeroMat <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values = 0, name = "zeroMat")
signA <- mxAlgebra( (-1)^omxLessThan(VAf, zeroMat) * (-1)^omxLessThan(VAm, zeroMat), name="signA") # 2/4: had to change this code to make a zero matrix. added zeroMat to parsZo
signC <- mxAlgebra( (-1)^omxLessThan(VCf, zeroMat) * (-1)^omxLessThan(VCm, zeroMat), name="signC")
# Calculate absolute covariation between Males and Females and then un-absoluting the product
covAos <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name="VAos")
covCos <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm)))), name="VCos")
# Calculate rg/rc from reparameterized model
pathRg <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
pathRc <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm))))/sqrt(VCf*(VCm+VCms)), name="rc")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covPf <- mxAlgebra( expression= VAf+VCf+VEf, name="Vf" )
covPm <- mxAlgebra( expression= VAm+VCm+VEm+VAms+VCms, name="Vm" )
covMZf <- mxAlgebra( expression= VAf+VCf, name="cMZf" )
covDZf <- mxAlgebra( expression= 0.5%x%VAf+ VCf, name="cDZf" )
covMZm <- mxAlgebra( expression= VAm+VCm+VAms+VCms, name="cMZm" )
covDZm <- mxAlgebra( expression= 0.5%x%VAm+ VCm+ 0.5%x%VAms+VCms, name="cDZm" )
covDZo <- mxAlgebra( expression= 0.5%x%VAos+VCos, name="cDZo" )
expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
expCovDZo <- mxAlgebra( expression= rbind( cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name="expCovDZo" )

# Create Data Objects for Multiple Groups
dataMZf <- mxData( observed=mzfData, type="raw" )
dataDZf <- mxData( observed=dzfData, type="raw" )
dataMZm <- mxData( observed=mzmData, type="raw" )
dataDZm <- mxData( observed=dzmData, type="raw" )
dataDZo <- mxData( observed=dzoData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZf <- mxExpectationNormal( covariance="expCovMZf", means="meanGf", dimnames=selVars )
expDZf <- mxExpectationNormal( covariance="expCovDZf", means="meanGf", dimnames=selVars )
expMZm <- mxExpectationNormal( covariance="expCovMZm", means="meanGm", dimnames=selVars )
expDZm <- mxExpectationNormal( covariance="expCovDZm", means="meanGm", dimnames=selVars )
expDZo <- mxExpectationNormal( covariance="expCovDZo", means="meanGo", dimnames=selVars )
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
parsZf <- list( covAf, covCf, covEf, covPf )
parsZm <- list( covAm, covCm, covEm, covPm, covAms, covCms )
parsZo <- list( parsZm, parsZf, zeroMat, signA, signC, covAos, covCos) #, pathRg, pathRc )
modelMZf <- mxModel( parsZf, meanGf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZf <- mxModel( parsZf, meanGf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf" )
modelMZm <- mxModel( parsZm, meanGm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm" )
modelDZm <- mxModel( parsZm, meanGm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm" )
modelDZo <- mxModel( parsZo, meanGo, covDZo, expCovDZo, dataDZo, expDZo, funML, name="DZo" )
multi <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

## Create Algebra for Standardization
matIf      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="If" )
matIm      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="Im" )

## Calculate genetic,  environmental, and phenotypic correlations
corAf      <- mxAlgebra( expression=solve(sqrt(If*VAf))%&%VAf, name ="rAf" )
corCf      <- mxAlgebra( expression=solve(sqrt(If*VCf))%&%VCf, name ="rCf" )
corEf      <- mxAlgebra( expression=solve(sqrt(If*VEf))%&%VEf, name ="rEf" )
corPf      <- mxAlgebra( expression=solve(sqrt(If*Vf))%&%Vf, name ="rPf")

corAm      <- mxAlgebra( expression=solve(sqrt(Im*VAm))%&%VAm, name ="rAm" )
corCm      <- mxAlgebra( expression=solve(sqrt(Im*VCm))%&%VCm, name ="rCm" )
corEm      <- mxAlgebra( expression=solve(sqrt(Im*VEm))%&%VEm, name ="rEm" )
corPm      <- mxAlgebra( expression=solve(sqrt(Im*Vm))%&%Vm, name ="rPm")

## Create Algebra for Variance Components ( formatting output)
#Unstandardized variance components
rowUV     <- rep('UV',nv)
colUV     <- rep(c('VA','VC','VE'),each=nv)
estUVf     <- mxAlgebra( expression=cbind(VAf,VCf,VEf), name="UVf", dimnames=list(rowUV,colUV) )
estUVm     <- mxAlgebra( expression=cbind(VAm,VCm,VEm), name="UVm", dimnames=list(rowUV,colUV) )

#Standardized variance components
rowSV     <- rep('SV',nv)
colSV     <- rep(c('SA','SC','SE'),each=nv)
estSVf     <- mxAlgebra( expression=cbind(VAf/Vf,VCf/Vf,VEf/Vf), name="SVf", dimnames=list(rowSV,colSV) )
estSVm     <- mxAlgebra( expression=cbind(VAm/Vm,VCm/Vm,VEm/Vm), name="SVm", dimnames=list(rowSV,colSV) )

#Correlations
rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rC','rE','rP'),each=nv)
estcorrf     <- mxAlgebra( expression=cbind(rAf,rCf,rEf,rPf), name="corrf", dimnames=list(rowcorr,colcorr) )
estcorrm     <- mxAlgebra( expression=cbind(rAm,rCm,rEm,rPm), name="corrm", dimnames=list(rowcorr,colcorr) )

## Create Confidence Interval Objects
ciUV      <- mxCI( c("UV") )
ciSV      <- mxCI( c("SVf", "SVm") )
cicorr      <- mxCI( c("corrf", "corrm") )

## Build Model with Confidence Intervals
calc      <- list( matIf,  corAf, corCf, corEf,corPf, estUVf, estSVf,estcorrf,
                   matIm,  corAm, corCm, corEm,corPm, estUVm, estSVm,estcorrm, ciSV, cicorr)
modelACE  <- mxModel( "twoACE5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACE model - No Sex differences ACE model
fitACE <- mxRun( modelACE, intervals=T )
sumACE <- summary( fitACE )

# Compare with Saturated Model
mxCompare( fitSAT, fitACE )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitACE)
fitEsts(fitACE)
fitEstCis(fitACE)

# Print Covariance & Correlation Matrices
fitACE$algebras$UVf
fitACE$algebras$UVm
fitACE$algebras$SVf
fitACE$algebras$SVm
fitACE$algebras$corrf
fitACE$algebras$corrm

# Print Comparative Fit Statistics
mxCompare( fitACEra, nested <- list( fitACErc, fitACEq, fitACE) )
#colnames <- c('name','VAf','VCf','VEf','SAf','SCf','SEf','VAm','VCm','VEm','SAm','SCm','SEm','rg','rc')
#print(cbind(rbind(fitACEra$name,fitACErc$name,fitACEq$name,fitACE$name),
            #round(rbind(fitACEra$US$result,fitACErc$US$result,fitACEq$US$result,fitACE$US$result),4)),quote=F)

# Based on this, it seems like ACE5 is the best fit --> males and females can be equated?

# ----------------------------------------------------------
# Test Significance of Sources of Variance of ACE model without Sex differences

# Run AE model
modelAE <- mxModel( fitACE, name="twoAE5" )
modelAE <- omxSetParameters( modelAE, labels=labLower("VCf",nv), free=FALSE, values=0 )
modelAE <- omxSetParameters( modelAE, labels=labLower("VCm",nv), free=FALSE, values=0 )
fitAE <- mxRun( modelAE, intervals=T )
fitGofs(fitAE); fitEsts(fitAE)

# Run CE model
modelCE <- mxModel( fitACE, name="twoCE" )
modelCE <- omxSetParameters( modelCE, labels=labLower("VAf",nv), free=FALSE, values=0 )
modelCE <- omxSetParameters( modelCE, labels=labLower("VCf",nv), free=TRUE, values=.6 )
modelCE <- omxSetParameters( modelCE, labels=labLower("VAm",nv), free=FALSE, values=0 )
modelCE <- omxSetParameters( modelCE, labels=labLower("VCm",nv), free=TRUE, values=.6 )
fitCE <- mxRun( modelCE, intervals=T )
fitGofs(fitCE); fitEstCis(fitCE)

# Run E model
modelE <- mxModel( fitCE, name="twoE" )
modelE <- modelCE <- omxSetParameters( modelE, labels=labLower("VCf",nv), free=FALSE, values=0 )
modelE <- modelCE <- omxSetParameters( modelE, labels=labLower("VCm",nv), free=FALSE, values=0 )
fitE <- mxRun( modelE, intervals=T )
fitGofs(fitE); fitEstCis(fitE)

# Print Comparative Fit Statistics
mxCompare( fitACE, nested <- list( fitAE, fitCE, fitE) )

# ----------------------------------------------------------------------------------------------------------------------
# ADE (no sex diff) MODEL STARTS HERE

# Set Starting Values
svMe      <- c(0.01, 0.01)                     # start value for means
svPa      <- .2                    # start value for variance component A
svPe      <- .5                   # start value for variance component E

# Create Algebra for expected Mean Matrices
meanGf      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSf","meanSDSf", "meanSRSf","meanSDSf"), name="meanGf" )
meanGm      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSm","meanSDSm", "meanSRSm","meanSDSm"), name="meanGm" )
meanGo      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                         labels=c("meanSRSo","meanSDSo", "meanSRSo","meanSDSo"), name="meanGo" )

# Create Matrices for Variance Components
covAf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAf",nv), name="VAf" ) #lbound=1e-4 lbound=valDiag(-1,ntv),
covDf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VDf",nv), name="VDf" )
covEf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                        label=labLower("VEf",nv), name="VEf" )
covAm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VAm",nv), name="VAm" )
covDm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                        label=labLower("VDm",nv), name="VDm" )
covEm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                        label=labLower("VEm",nv), name="VEm" )
covAms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
                         label=labLower("VAms",nv), name="VAms" ) #lbound=valDiag(-0.5,ntv)
covDms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
                         label=labLower("VDms",nv), name="VDms" )

# Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
# signA <- mxAlgebra( ((-1)^omxLessThan(VAf,0))*((-1)^omxLessThan(VAm,0)), name="signA") # original code, not working
# signC <- mxAlgebra( ((-1)^omxLessThan(VCf,0))*((-1)^omxLessThan(VCm,0)), name="signC")
zeroMat <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values = 0, name = "zeroMat")
signA <- mxAlgebra( (-1)^omxLessThan(VAf, zeroMat) * (-1)^omxLessThan(VAm, zeroMat), name="signA") # 2/4: had to change this code to make a zero matrix. added zeroMat to parsZo
signD <- mxAlgebra( (-1)^omxLessThan(VDf, zeroMat) * (-1)^omxLessThan(VDm, zeroMat), name="signD")
# Calculate absolute covariation between Males and Females and then un-absoluting the product
covAos <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name="VAos")
covDos <- mxAlgebra( signD*(sqrt(abs(VDf))*t(sqrt(abs(VDm)))), name="VDos")
# Calculate rg/rc from reparameterized model
pathRg <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
pathRd <- mxAlgebra( signD*(sqrt(abs(VDf))*t(sqrt(abs(VDm))))/sqrt(VDf*(VDm+VDms)), name="rd")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covPf <- mxAlgebra( expression= VAf+VDf+VEf, name="Vf" )
covPm <- mxAlgebra( expression= VAm+VDm+VEm+VAms+VDms, name="Vm" )
covMZf <- mxAlgebra( expression= VAf+VDf, name="cMZf" )
covDZf <- mxAlgebra( expression= 0.5%x%VAf+ 0.25%x%VDf, name="cDZf" )
covMZm <- mxAlgebra( expression= VAm+VDm+VAms+VDms, name="cMZm" )
covDZm <- mxAlgebra( expression= 0.5%x%VAm+ 0.25%x%VDm+ 0.5%x%VAms+ 0.25%x%VDms, name="cDZm" )
covDZo <- mxAlgebra( expression= 0.5%x%VAos+VDos, name="cDZo" )
expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
expCovDZo <- mxAlgebra( expression= rbind( cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name="expCovDZo" )

# Create Data Objects for Multiple Groups
dataMZf <- mxData( observed=mzfData, type="raw" )
dataDZf <- mxData( observed=dzfData, type="raw" )
dataMZm <- mxData( observed=mzmData, type="raw" )
dataDZm <- mxData( observed=dzmData, type="raw" )
dataDZo <- mxData( observed=dzoData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZf <- mxExpectationNormal( covariance="expCovMZf", means="meanGf", dimnames=selVars )
expDZf <- mxExpectationNormal( covariance="expCovDZf", means="meanGf", dimnames=selVars )
expMZm <- mxExpectationNormal( covariance="expCovMZm", means="meanGm", dimnames=selVars )
expDZm <- mxExpectationNormal( covariance="expCovDZm", means="meanGm", dimnames=selVars )
expDZo <- mxExpectationNormal( covariance="expCovDZo", means="meanGo", dimnames=selVars )
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
parsZf <- list( covAf, covDf, covEf, covPf )
parsZm <- list( covAm, covDm, covEm, covPm, covAms, covDms )
parsZo <- list( parsZm, parsZf, zeroMat, signA, signD, covAos, covDos, pathRg, pathRd )
modelMZf <- mxModel( parsZf, meanGf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf" )
modelDZf <- mxModel( parsZf, meanGf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf" )
modelMZm <- mxModel( parsZm, meanGm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm" )
modelDZm <- mxModel( parsZm, meanGm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm" )
modelDZo <- mxModel( parsZo, meanGo, covDZo, expCovDZo, dataDZo, expDZo, funML, name="DZo" )
multi <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )

## Create Algebra for Standardization
matIf      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="If" )
matIm      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="Im" )

## Calculate genetic,  environmental, and phenotypic correlations
corAf      <- mxAlgebra( expression=solve(sqrt(If*VAf))%&%VAf, name ="rAf" )
corDf      <- mxAlgebra( expression=solve(sqrt(If*VDf))%&%VDf, name ="rDf" )
corEf      <- mxAlgebra( expression=solve(sqrt(If*VEf))%&%VEf, name ="rEf" )
corPf      <- mxAlgebra( expression=solve(sqrt(If*Vf))%&%Vf, name ="rPf")

corAm      <- mxAlgebra( expression=solve(sqrt(Im*VAm))%&%VAm, name ="rAm" )
corDm      <- mxAlgebra( expression=solve(sqrt(Im*VDm))%&%VDm, name ="rDm" )
corEm      <- mxAlgebra( expression=solve(sqrt(Im*VEm))%&%VEm, name ="rEm" )
corPm      <- mxAlgebra( expression=solve(sqrt(Im*Vm))%&%Vm, name ="rPm")

## Create Algebra for Variance Components ( formatting output)
#Unstandardized variance components
rowUV     <- rep('UV',nv)
colUV     <- rep(c('VA','VD','VE'),each=nv)
estUVf     <- mxAlgebra( expression=cbind(VAf,VDf,VEf), name="UVf", dimnames=list(rowUV,colUV) )
estUVm     <- mxAlgebra( expression=cbind(VAm,VDm,VEm), name="UVm", dimnames=list(rowUV,colUV) )

#Standardized variance components
rowSV     <- rep('SV',nv)
colSV     <- rep(c('SA','SD','SE'),each=nv)
estSVf     <- mxAlgebra( expression=cbind(VAf/Vf,VDf/Vf,VEf/Vf), name="SVf", dimnames=list(rowSV,colSV) )
estSVm     <- mxAlgebra( expression=cbind(VAm/Vm,VDm/Vm,VEm/Vm), name="SVm", dimnames=list(rowSV,colSV) )

#Correlations
rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rD','rE','rP'),each=nv)
estcorrf     <- mxAlgebra( expression=cbind(rAf,rDf,rEf,rPf), name="corrf", dimnames=list(rowcorr,colcorr) )
estcorrm     <- mxAlgebra( expression=cbind(rAm,rDm,rEm,rPm), name="corrm", dimnames=list(rowcorr,colcorr) )

## Create Confidence Interval Objects
ciUV      <- mxCI( c("UV") )
ciSV      <- mxCI( c("SVf", "SVm") )
cicorr      <- mxCI( c("corrf", "corrm") )

## Build Model with Confidence Intervals
calc      <- list( matIf,  corAf, corDf, corEf,corPf, estUVf, estSVf,estcorrf,
                   matIm,  corAm, corDm, corEm,corPm, estUVm, estSVm,estcorrm, ciSV, cicorr)
modelADE  <- mxModel( "twoADE5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ADE model - No Sex differences ADE model
fitADE <- mxRun( modelADE, intervals=T )
sumADE <- summary( fitADE )

# Compare with Saturated Model
mxCompare( fitSAT, nested <- list( fitACE, fitADE) )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(fitADE)
fitEsts(fitADE)
fitEstCis(fitADE)

# Print Covariance & Correlation Matrices
fitADE$algebras$UVf
fitADE$algebras$UVm
fitADE$algebras$SVf
fitADE$algebras$SVm
fitADE$algebras$corrf
fitADE$algebras$corrm

