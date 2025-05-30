# Load Libraries & Options
rm(list=ls())
library(OpenMx)
library(dplyr)
library(tidyverse)
source("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/02_scripts/miFunctions.R")

# ----------------------------------------------------------------------------------------------------------------------
# PREPARE DATA

# Load Data
vcu25w_imp_07 <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/01_tidy_data/vcu25w_imp_07.csv")
vcu25w_imp_log_07 <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/01_tidy_data/vcu25w_imp_log_07.csv")
vcu25w_imp_log_inc_07 <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/01_tidy_data/vcu25w_imp_log_inc_07.csv")
vcu25w_med_imp_log_07 <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/01_tidy_data/vcu25w_med_imp_log_07.csv")
vcu25w_imp_06 <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/01_tidy_data/vcu25w_imp_06.csv")

# Select Variables for Analysis
vars <- c('rsrs_imp_new','rsds_imp_new') # list of variables names
nv <- 2 # number of variables
ntv <- nv*2 # number of total variables
selVars <- paste(vars,c(rep("_T1",nv),rep("_T2",nv)),sep="")

# mzfData <- subset(vcu25w_imp_log_inc_07, (zygo == "MZ" & sex_T1 == 2 & sex_T2 == 2), select = selVars)
# dzfData <- subset(vcu25w_imp_log_inc_07, (zygo == "DZ" & sex_T1 == 2 & sex_T2 == 2), select = selVars)
# mzmData <- subset(vcu25w_imp_log_inc_07, (zygo == "MZ" & sex_T1 == 1 & sex_T2 == 1), select = selVars)
# dzmData <- subset(vcu25w_imp_log_inc_07, (zygo == "DZ" & sex_T1 == 1 & sex_T2 == 1), select = selVars)
# dzoData <- subset(vcu25w_imp_log_inc_07, (zygo == "DZ" & sex_T1 != sex_T2), select = selVars)

mzfData <- subset(vcu25w_imp_new_2, (zygo == "MZ" & sex_T1 == 2 & sex_T2 == 2), select = selVars)
dzfData <- subset(vcu25w_imp_new_2, (zygo == "DZ" & sex_T1 == 2 & sex_T2 == 2), select = selVars)
mzmData <- subset(vcu25w_imp_new_2, (zygo == "MZ" & sex_T1 == 1 & sex_T2 == 1), select = selVars)
dzmData <- subset(vcu25w_imp_new_2, (zygo == "DZ" & sex_T1 == 1 & sex_T2 == 1), select = selVars)
dzoData <- subset(vcu25w_imp_new_2, (zygo == "DZ" & sex_T1 != sex_T2), select = selVars)

# mzfData <- subset(vcu25w_med_imp_log_07, (zygo == "MZ" & sex_T1 == 2 & sex_T2 == 2), select = selVars)
# dzfData <- subset(vcu25w_med_imp_log_07, (zygo == "DZ" & sex_T1 == 2 & sex_T2 == 2), select = selVars)
# mzmData <- subset(vcu25w_med_imp_log_07, (zygo == "MZ" & sex_T1 == 1 & sex_T2 == 1), select = selVars)
# dzmData <- subset(vcu25w_med_imp_log_07, (zygo == "DZ" & sex_T1 == 1 & sex_T2 == 1), select = selVars)
# dzoData <- subset(vcu25w_med_imp_log_07, (zygo == "DZ" & sex_T1 != sex_T2), select = selVars)

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
# lbVa <- .0001 # lower bound for variances

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
covMZf <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvMZf, name="covMZf" ) # lbound=valDiag(lbVa,ntv),
covDZf <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvDZf, name="covDZf" )
covMZm <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvMZm, name="covMZm" )
covDZm <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvDZm, name="covDZm" )
covDZo <- mxMatrix( type="Symm", nrow=ntv, ncol=ntv, free=TRUE, values=valDiag(svVa,ntv), labels=labCvDZo, name="covDZo" )

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
modelSAT <- mxModel( "vcu_twoSAT", modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, ciCov, ciMean )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run Saturated Model
vcu_fitSAT <- mxRun( modelSAT, intervals=F )
sumSAT <- summary( vcu_fitSAT )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(vcu_fitSAT)
fitEsts(vcu_fitSAT)
mxGetExpected(vcu_fitSAT, "means")
mxGetExpected(vcu_fitSAT, "covariance")

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# Constrain expected Means to be equal across Twin Order
modelEMO <- mxModel( vcu_fitSAT, name="vcu_twoEMO" )
for (i in 1:nv) {
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeMZf[nv+i],labMeMZf[i]), free=TRUE, values=svMe, newlabels=labMeMZf[i] )
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeDZf[nv+i],labMeDZf[i]), free=TRUE, values=svMe, newlabels=labMeDZf[i] )
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeMZm[nv+i],labMeMZm[i]), free=TRUE, values=svMe, newlabels=labMeMZm[i] )
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeDZm[nv+i],labMeDZm[i]), free=TRUE, values=svMe, newlabels=labMeDZm[i] )
  modelEMO <- omxSetParameters( modelEMO, label=c(labMeDZo[nv+i],labMeDZo[i]), free=TRUE, values=svMe, newlabels=labMeDZo[i] )}
vcu_fitEMO <- mxRun( modelEMO, intervals=F )
fitGofs(vcu_fitEMO); fitEsts(vcu_fitEMO)

# Constrain expected Means and Variances to be equal across Twin Order
modelEMVO <- mxModel( vcu_fitEMO, name="vcu_twoEMVO" )
for (i in 1:nv) {
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaMZf[nv+i],labVaMZf[i]), free=TRUE, values=c(.75,.75), newlabels=labVaMZf[i] )
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaDZf[nv+i],labVaDZf[i]), free=TRUE, values=c(.72,.72), newlabels=labVaDZf[i] )
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaMZm[nv+i],labVaMZm[i]), free=TRUE, values=c(.87,.87), newlabels=labVaMZm[i] )
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaDZm[nv+i],labVaDZm[i]), free=TRUE, values=c(0.88,0.88), newlabels=labVaDZm[i] )
  modelEMVO <- omxSetParameters( modelEMVO, label=c(labVaDZo[nv+i],labVaDZo[i]), free=TRUE, values=c(0.83,0.83), newlabels=labVaDZo[i] )}
vcu_fitEMVO <- mxRun( modelEMVO, intervals=F )
fitGofs(vcu_fitEMVO); fitEsts(vcu_fitEMVO)

# Constrain expected Means and Variances to be equal across Twin Order and Zygosity
modelEMVZ <- mxModel( vcu_fitEMO, name="vcu_twoEMVZ" )
for (i in 1:nv) {
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labMeMZf[i],labMeDZf[i]), free=TRUE, values=svMe, newlabels=labMeZf[i] ) # Added f to this to try to run EMVP
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labVaMZf[i],labVaDZf[i]), free=TRUE, values=c(.81,.81), newlabels=labVaZf[i] )
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labMeMZm[i],labMeDZm[i]), free=TRUE, values=svMe, newlabels=labMeZm[i] ) # Added m to this to try to run EMVP
  modelEMVZ <- omxSetParameters( modelEMVZ, label=c(labVaMZm[i],labVaDZm[i]), free=TRUE, values=c(.81,.81), newlabels=labVaZm[i] )}
vcu_fitEMVZ <- mxRun( modelEMVZ, intervals=F )
fitGofs(vcu_fitEMVZ); fitEsts(vcu_fitEMVZ)

# Constrain expected Means and Variances to be equal across twin order and zygosity and SS/OS
modelEMVP <- mxModel( vcu_fitEMVZ, name="vcu_twoEMVP" ) # 1/29 made strict=FALSE, not sure what this means
for (i in 1:nv) {
  modelEMVP <- omxSetParameters( modelEMVP, label=c(labMeZf[i],labMeDZo[i]), free=TRUE, strict=FALSE, values=svMe, newlabels=labMeZf[i] )
  modelEMVP <- omxSetParameters( modelEMVP, label=c(labVaZf[i],labVaDZo[i]), free=TRUE, strict=FALSE, values=c(.77,.77), newlabels=labVaZf[i] )
  modelEMVP <- omxSetParameters( modelEMVP, label=c(labMeZm[i],labMeDZo[i]), free=TRUE, strict=FALSE, values=svMe, newlabels=labMeZm[i] )
  modelEMVP <- omxSetParameters( modelEMVP, label=c(labVaZm[i],labVaDZo[i]), free=TRUE, strict=FALSE, values=c(.77,.77), newlabels=labVaZm[i] )}
vcu_fitEMVP <- mxRun( modelEMVP, intervals=F )
fitGofs(vcu_fitEMVP); fitEsts(vcu_fitEMVP)

# Constrain expected Means and Variances to be equal across twin order and zygosity and SS/OS and sex
modelEMVS <- mxModel( vcu_fitEMVP, name="vcu_twoEMVS" )
for (i in 1:nv) {
  modelEMVS <- omxSetParameters( modelEMVS, label=c(labMeZf[i],labMeZm[i]), free=TRUE, values=svMe, newlabels=labMeZ[i] )
  modelEMVS <- omxSetParameters( modelEMVS, label=c(labVaZf[i],labVaZm[i]), free=TRUE, values=c(.77,.77), newlabels=labVaZ[i] )}
vcu_fitEMVS <- mxRun( modelEMVS, intervals=F )
fitGofs(vcu_fitEMVS); fitEsts(vcu_fitEMVS)

# Print Comparative Fit Statistics
mxCompare( vcu_fitSAT, subs <- list(vcu_fitEMO, vcu_fitEMVZ, vcu_fitEMVP, vcu_fitEMVS) ) #  vcu_fitEMVO,

# ----------------------------------------------------------------------------------------------------------------------
# ACEra MODEL STARTS HERE

# Set Starting Values
svMe      <- c(0.01, 0.01)                     # start value for means
svPa      <- .4                    # start value for variance component A
svPe      <- .8                   # start value for variance component E

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
modelACEra  <- mxModel( "vcu_twoACEra5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACEra Model - Qualitative (Ra) & Quantative Sex Differences ACE model
vcu_fitACEra <- mxTryHard( modelACEra, intervals=T )
sumACEra <- summary( vcu_fitACEra )

# Compare with Saturated Model
mxCompare( vcu_fitSAT, vcu_fitACEra )

# Print Goodness-of-fit Statistics & Parameter Estimates
fitGofs(vcu_fitACEra)
fitEsts(vcu_fitACEra)
fitEstCis(vcu_fitACEra)

# Print Covariance & Correlation Matrices
vcu_fitACEra$algebras$UVf
vcu_fitACEra$algebras$UVm
vcu_fitACEra$algebras$SVf
vcu_fitACEra$algebras$SVm
vcu_fitACEra$algebras$corrf
vcu_fitACEra$algebras$corrm

# ----------------------------------------------------------------------------------------------------------------------
# RUN SUBMODELS

# # Run ACErc Model - Qualitative (Rc) & Quantative Sex Differences ACE model
# modelACErc <- mxModel( vcu_fitACEra, name="ACErc5" )
# modelACErc <- omxSetParameters(modelACErc, labels=c("VAms11", "VAms21", "VAms22"), free=FALSE, values=valDiag(0,nv))
# modelACErc <- omxSetParameters(modelACErc, labels=c("VCms11", "VCms21", "VCms22"), free=TRUE, values=valDiag(0,nv))
# modelACErc <- omxSetParameters(modelACErc, labels=c("VEf11", "VEf21", "VEf22", "VEm11", "VEm21", "VEm22"), free=TRUE, values=valDiag(svPe,nv))
# vcu_fitACErc <- mxTryHard( modelACErc, intervals=T )
# fitGofs(vcu_fitACErc); fitEsts(vcu_fitACErc)

# ----------------------------------------------------------------------------------------------------------------------
# ACErc MODEL HERE

# Set Starting Values 
svMe      <- c(0.01, 0.01)                     # start value for means
svPa      <- .4                    # start value for variance component A
svPe      <- .8                   # start value for variance component E

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
modelACErc  <- mxModel( "vcu_twoACErc5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACErc Model - Qualitative (Rc) & Quantative Sex Differences ACE model
vcu_fitACErc <- mxTryHard( modelACErc, intervals=T )
sumACErc <- summary( vcu_fitACErc )


# ----------------------------------------------------------------------------------------------------------------------
# ACEq MODEL STARTS HERE

# Set Starting Values
svMe      <- c(0.01, 0.01)                     # start value for means
svPa      <- .4                    # start value for variance component A
svPe      <- .8                   # start value for variance component E

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
modelACEq  <- mxModel( "vcu_twoACEq5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACEq model - Quantitative non-scalar Sex Differences ACE model
vcu_fitACEq <- mxTryHard( modelACEq, intervals=T )
sumACEq <- summary( vcu_fitACEq )

# ----------------------------------------------------------------------------------------------------------------------
# ACE (no sex diff) MODEL STARTS HERE

# vcu25w_imp_log_inc_sex_07 <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/01_tidy_data/vcu25w_imp_log_inc_sex_07.csv")
# vcu25w_med_imp_log_sex_07 <- read_csv("/Volumes/SoCoDeN/SSCDN_Private/data/VCU_MATR/Amanda_analysis/01_tidy_data/vcu25w_med_imp_log_07.csv")

# mzData <- subset(vcu25w_imp_log_inc_sex_07, (zygo == "MZ"), select = selVars)
# dzData <- subset(vcu25w_imp_log_inc_sex_07, (zygo == "DZ"), select = selVars)

mzData <- subset(vcu25w_imp_new_sex_2, (zygo == "MZ"), select = selVars)
dzData <- subset(vcu25w_imp_new_sex_2, (zygo == "DZ"), select = selVars)

# mzData <- subset(vcu25w_med_imp_log_sex_07, (zygo == "MZ"), select = selVars)
# dzData <- subset(vcu25w_med_imp_log_sex_07, (zygo == "DZ"), select = selVars)

colMeans(mzData,na.rm=TRUE)
colMeans(dzData,na.rm=TRUE)

# Set Starting Values
svMe      <- c(-0.03, 0.06)                     # start value for means (changed from 0.01, 0.01)
svPa      <- .4                    # start value for variance component A
svPe      <- .8                   # start value for variance component E

# Create Algebra for expected Mean Matrices
meanG      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
                        labels=c("meanSRS","meanSDS", "meanSRS","meanSDS"), name="meanG" )
# meanGf      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
#                          labels=c("meanSRSf","meanSDSf", "meanSRSf","meanSDSf"), name="meanGf" )
# meanGm      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
#                          labels=c("meanSRSm","meanSDSm", "meanSRSm","meanSDSm"), name="meanGm" )
# meanGo      <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe,
#                          labels=c("meanSRSo","meanSDSo", "meanSRSo","meanSDSo"), name="meanGo" )

# Create Matrices for Variance Components
covA      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                       label=labLower("VA",nv), name="VA" ) # lbound=valDiag(1e-4,ntv), 
covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
                       label=labLower("VC",nv), name="VC" )
# covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
#                        label=labLower("VC",nv), lbound=valDiag(1e-4,ntv), name="VC" ) # added lower bound
# covC      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
#                        label=labLower("VC",nv), name="VC" ) # fix to zero --> estimates made more sense with lower bound
covE      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
                       label=labLower("VE",nv), name="VE" )
# covAf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
#                         label=labLower("VAf",nv), name="VAf" ) # lbound=valDiag(1e-4,ntv), 
# covCf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
#                         label=labLower("VCf",nv), name="VCf" )
# covEf      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
#                         label=labLower("VEf",nv), name="VEf" )
# covAm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
#                         label=labLower("VAm",nv), name="VAm" )
# covCm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv),
#                         label=labLower("VCm",nv), name="VCm" )
# covEm      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPe,nv),
#                         label=labLower("VEm",nv), name="VEm" )
# covAms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
#                          label=labLower("VAms",nv), name="VAms" ) 
# covCms      <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=FALSE, values=valDiag(0,nv),
#                          label=labLower("VCms",nv), name="VCms" )

# # Produce a vector which is = to 1 if variance is positive and -1 if variance is negative
# zeroMat <- mxMatrix(type = "Full", nrow = 2, ncol = 2, values = 0, name = "zeroMat")
# signA <- mxAlgebra( (-1)^omxLessThan(VAf, zeroMat) * (-1)^omxLessThan(VAm, zeroMat), name="signA") # 2/4: had to change this code to make a zero matrix. added zeroMat to parsZo
# signC <- mxAlgebra( (-1)^omxLessThan(VCf, zeroMat) * (-1)^omxLessThan(VCm, zeroMat), name="signC")
# # Calculate absolute covariation between Males and Females and then un-absoluting the product
# covAos <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm)))), name="VAos")
# covCos <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm)))), name="VCos")
# # Calculate rg/rc from reparameterized model
# pathRg <- mxAlgebra( signA*(sqrt(abs(VAf))*t(sqrt(abs(VAm))))/sqrt(VAf*(VAm+VAms)), name="rg")
# pathRc <- mxAlgebra( signC*(sqrt(abs(VCf))*t(sqrt(abs(VCm))))/sqrt(VCf*(VCm+VCms)), name="rc")

# Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
covP <- mxAlgebra( expression= VA+VC+VE, name="V" )
# covPf <- mxAlgebra( expression= VAf+VCf+VEf, name="Vf" )
# covPm <- mxAlgebra( expression= VAm+VCm+VEm+VAms+VCms, name="Vm" )
covMZ <- mxAlgebra( expression= VA+VC, name="cMZ" )
covDZ <- mxAlgebra( expression= 0.5%x%VA+ VC, name="cDZ" )
# covMZf <- mxAlgebra( expression= VAf+VCf, name="cMZf" )
# covDZf <- mxAlgebra( expression= 0.5%x%VAf+ VCf, name="cDZf" )
# covMZm <- mxAlgebra( expression= VAm+VCm+VAms+VCms, name="cMZm" )
# covDZm <- mxAlgebra( expression= 0.5%x%VAm+ VCm+ 0.5%x%VAms+VCms, name="cDZm" )
# covDZo <- mxAlgebra( expression= 0.5%x%VAos+VCos, name="cDZo" )
expCovMZ <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
expCovDZ <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
# expCovMZf <- mxAlgebra( expression= rbind( cbind(Vf, cMZf), cbind(t(cMZf), Vf)), name="expCovMZf" )
# expCovDZf <- mxAlgebra( expression= rbind( cbind(Vf, cDZf), cbind(t(cDZf), Vf)), name="expCovDZf" )
# expCovMZm <- mxAlgebra( expression= rbind( cbind(Vm, cMZm), cbind(t(cMZm), Vm)), name="expCovMZm" )
# expCovDZm <- mxAlgebra( expression= rbind( cbind(Vm, cDZm), cbind(t(cDZm), Vm)), name="expCovDZm" )
# expCovDZo <- mxAlgebra( expression= rbind( cbind(Vf, cDZo), cbind(t(cDZo), Vm)), name="expCovDZo" )

# Create Data Objects for Multiple Groups
dataMZ <- mxData( observed=mzData, type="raw" )
dataDZ <- mxData( observed=dzData, type="raw" )
# dataMZf <- mxData( observed=mzfData, type="raw" )
# dataDZf <- mxData( observed=dzfData, type="raw" )
# dataMZm <- mxData( observed=mzmData, type="raw" )
# dataDZm <- mxData( observed=dzmData, type="raw" )
# dataDZo <- mxData( observed=dzoData, type="raw" )

# Create Expectation Objects for Multiple Groups
expMZ <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
expDZ <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
# expMZf <- mxExpectationNormal( covariance="expCovMZf", means="meanGf", dimnames=selVars )
# expDZf <- mxExpectationNormal( covariance="expCovDZf", means="meanGf", dimnames=selVars )
# expMZm <- mxExpectationNormal( covariance="expCovMZm", means="meanGm", dimnames=selVars )
# expDZm <- mxExpectationNormal( covariance="expCovDZm", means="meanGm", dimnames=selVars )
# expDZo <- mxExpectationNormal( covariance="expCovDZo", means="meanGo", dimnames=selVars )
funML <- mxFitFunctionML()

# Create Model Objects for Multiple Groups
parsZ <- list( covA, covC, covE, covP )
# parsZf <- list( covAf, covCf, covEf, covPf )
# parsZm <- list( covAm, covCm, covEm, covPm, covAms, covCms )
# parsZo <- list( parsZm, parsZf, zeroMat, signA, signC, covAos, covCos) #, pathRg, pathRc )
modelMZ <- mxModel( parsZ, meanG, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
modelDZ <- mxModel( parsZ, meanG, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
# modelMZf <- mxModel( parsZf, meanGf, covMZf, expCovMZf, dataMZf, expMZf, funML, name="MZf" )
# modelDZf <- mxModel( parsZf, meanGf, covDZf, expCovDZf, dataDZf, expDZf, funML, name="DZf" )
# modelMZm <- mxModel( parsZm, meanGm, covMZm, expCovMZm, dataMZm, expMZm, funML, name="MZm" )
# modelDZm <- mxModel( parsZm, meanGm, covDZm, expCovDZm, dataDZm, expDZm, funML, name="DZm" )
# modelDZo <- mxModel( parsZo, meanGo, covDZo, expCovDZo, dataDZo, expDZo, funML, name="DZo" )
# multi <- mxFitFunctionMultigroup( c("MZf","DZf","MZm","DZm","DZo") )
multi <- mxFitFunctionMultigroup( c("MZ","DZ") )

## Create Algebra for Standardization
matI      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I" )
# matIf      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="If" )
# matIm      <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="Im" )

## Calculate genetic,  environmental, and phenotypic correlations
corA      <- mxAlgebra( expression=solve(sqrt(I*VA))%&%VA, name ="rA" )
corC      <- mxAlgebra( expression=solve(sqrt(I*VC))%&%VC, name ="rC" )
corE      <- mxAlgebra( expression=solve(sqrt(I*VE))%&%VE, name ="rE" )
corP      <- mxAlgebra( expression=solve(sqrt(I*V))%&%V, name ="rP")

# corAf      <- mxAlgebra( expression=solve(sqrt(If*VAf))%&%VAf, name ="rAf" )
# corCf      <- mxAlgebra( expression=solve(sqrt(If*VCf))%&%VCf, name ="rCf" )
# corEf      <- mxAlgebra( expression=solve(sqrt(If*VEf))%&%VEf, name ="rEf" )
# corPf      <- mxAlgebra( expression=solve(sqrt(If*Vf))%&%Vf, name ="rPf")
# 
# corAm      <- mxAlgebra( expression=solve(sqrt(Im*VAm))%&%VAm, name ="rAm" )
# corCm      <- mxAlgebra( expression=solve(sqrt(Im*VCm))%&%VCm, name ="rCm" )
# corEm      <- mxAlgebra( expression=solve(sqrt(Im*VEm))%&%VEm, name ="rEm" )
# corPm      <- mxAlgebra( expression=solve(sqrt(Im*Vm))%&%Vm, name ="rPm")

## Create Algebra for Variance Components ( formatting output)
#Unstandardized variance components
rowUV     <- rep('UV',nv)
colUV     <- rep(c('VA','VC','VE'),each=nv)
estUV     <- mxAlgebra( expression=cbind(VA,VC,VE), name="UV", dimnames=list(rowUV,colUV) )
# estUVf     <- mxAlgebra( expression=cbind(VAf,VCf,VEf), name="UVf", dimnames=list(rowUV,colUV) )
# estUVm     <- mxAlgebra( expression=cbind(VAm,VCm,VEm), name="UVm", dimnames=list(rowUV,colUV) )

#Standardized variance components
rowSV     <- rep('SV',nv)
colSV     <- rep(c('SA','SC','SE'),each=nv)
estSV     <- mxAlgebra( expression=cbind(VA/V,VC/V,VE/V), name="SV", dimnames=list(rowSV,colSV) )
# estSVf     <- mxAlgebra( expression=cbind(VAf/Vf,VCf/Vf,VEf/Vf), name="SVf", dimnames=list(rowSV,colSV) )
# estSVm     <- mxAlgebra( expression=cbind(VAm/Vm,VCm/Vm,VEm/Vm), name="SVm", dimnames=list(rowSV,colSV) )

#Correlations
rowcorr     <- rep('corr',nv)
colcorr     <- rep(c('rA','rC','rE','rP'),each=nv)
estcorr     <- mxAlgebra( expression=cbind(rA,rC,rE,rP), name="corr", dimnames=list(rowcorr,colcorr) )
# estcorrf     <- mxAlgebra( expression=cbind(rAf,rCf,rEf,rPf), name="corrf", dimnames=list(rowcorr,colcorr) )
# estcorrm     <- mxAlgebra( expression=cbind(rAm,rCm,rEm,rPm), name="corrm", dimnames=list(rowcorr,colcorr) )

## Create Confidence Interval Objects
ciUV      <- mxCI( c("UV") )
ciSV      <- mxCI( c("SV") )
cicorr      <- mxCI( c("corr") )
# ciSV      <- mxCI( c("SVf", "SVm") )
# cicorr      <- mxCI( c("corrf", "corrm") )

## Build Model with Confidence Intervals
calc      <- list( matI,  corA, corC, corE,corP, estUV, estSV, estcorr, ciSV, cicorr)
modelACE  <- mxModel( "vcu_twoACE5", parsZ, modelMZ, modelDZ, multi, calc )
# calc      <- list( matIf,  corAf, corCf, corEf,corPf, estUVf, estSVf,estcorrf,
#                    matIm,  corAm, corCm, corEm,corPm, estUVm, estSVm,estcorrm, ciSV, cicorr)
# modelACE  <- mxModel( "vcu_twoACE5", parsZf, parsZm, parsZo, modelMZf, modelDZf, modelMZm, modelDZm, modelDZo, multi, calc )

# ----------------------------------------------------------------------------------------------------------------------
# RUN MODEL

# Run ACE model - No Sex differences ACE model
vcu_fitACE <- mxTryHard( modelACE, intervals=T )
sumACE <- summary( vcu_fitACE )


# Print Comparative Fit Statistics
mxCompare( vcu_fitSAT, nested <- list( vcu_fitACEra, vcu_fitACErc, vcu_fitACEq, vcu_fitACE) )
#colnames <- c('name','VAf','VCf','VEf','SAf','SCf','SEf','VAm','VCm','VEm','SAm','SCm','SEm','rg','rc')



# Run AE model
modelAE <- mxModel( vcu_fitACE, name="vcu_twoAE5" )
modelAE <- omxSetParameters( modelAE, labels=labLower("VC",nv), free=FALSE, values=0 )
vcu_fitAE <- mxRun( modelAE, intervals=T )
fitGofs(vcu_fitAE); fitEsts(vcu_fitAE)

# Run CE model
modelCE <- mxModel( vcu_fitACE, name="vcu_twoCE" )
modelCE <- omxSetParameters( modelCE, labels=labLower("VA",nv), free=FALSE, values=0 )
modelCE <- omxSetParameters( modelCE, labels=labLower("VC",nv), free=TRUE, values=.6 )
vcu_fitCE <- mxRun( modelCE, intervals=T )
fitGofs(vcu_fitCE); fitEstCis(vcu_fitCE)

# Run E model
modelE <- mxModel( vcu_fitAE, name="vcu_twoE" )
modelE <- modelCE <- omxSetParameters( modelE, labels=labLower("VC",nv), free=FALSE, values=0 )
vcu_fitE <- mxRun( modelE, intervals=T )
fitGofs(vcu_fitE); fitEstCis(vcu_fitE)

# Print Comparative Fit Statistics
mxCompare( vcu_fitSAT, nested <- list( vcu_fitACE, vcu_fitAE, vcu_fitCE, vcu_fitE) )


# # ----------------------------------------------------------------------------------------------------------------------
# # 2-GROUP (NO SEX-LIM) ADE MODEL STARTS HERE
# 
# # Set Starting Values
# svMe <- c(0.01, 0.01) # start value for means
# svPa <- .4 # start value for path coefficient
# svPe <- .8 # start value for path coefficient for e
# 
# # ----------------------------------------------------------------------------------------------------------------------
# # PREPARE MODEL
# 
# # Create Algebra for expected Mean Matrices
# meanG <- mxMatrix( type="Full", nrow=1, ncol=ntv, free=TRUE, values=svMe, labels=labVars("mean",vars), name="meanG" )
# 
# # Create Matrices for Variance Components
# covA <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labLower("VA",nv), name="VA" )
# covD <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labLower("VD",nv), name="VD" )
# covE <- mxMatrix( type="Symm", nrow=nv, ncol=nv, free=TRUE, values=valDiag(svPa,nv), label=labLower("VE",nv), name="VE" )
# 
# # Create Algebra for expected Variance/Covariance Matrices in MZ & DZ twins
# covP <- mxAlgebra( expression= VA+VD+VE, name="V" )
# covMZ <- mxAlgebra( expression= VA+VD, name="cMZ" )
# covDZ <- mxAlgebra( expression= 0.5%x%VA+ 0.25%x%VD, name="cDZ" )
# expCovMZ <- mxAlgebra( expression= rbind( cbind(V, cMZ), cbind(t(cMZ), V)), name="expCovMZ" )
# expCovDZ <- mxAlgebra( expression= rbind( cbind(V, cDZ), cbind(t(cDZ), V)), name="expCovDZ" )
# 
# # Create Data Objects for Multiple Groups
# dataMZ <- mxData( observed=mzData, type="raw" )
# dataDZ <- mxData( observed=dzData, type="raw" )
# 
# # Create Expectation Objects for Multiple Groups
# expMZ <- mxExpectationNormal( covariance="expCovMZ", means="meanG", dimnames=selVars )
# expDZ <- mxExpectationNormal( covariance="expCovDZ", means="meanG", dimnames=selVars )
# funML <- mxFitFunctionML()
# 
# # Create Model Objects for Multiple Groups
# pars <- list( meanG, covA, covD, covE, covP )
# modelMZ <- mxModel( pars, covMZ, expCovMZ, dataMZ, expMZ, funML, name="MZ" )
# modelDZ <- mxModel( pars, covDZ, expCovDZ, dataDZ, expDZ, funML, name="DZ" )
# multi <- mxFitFunctionMultigroup( c("MZ","DZ") )
# 
# # Create Algebra for Standardization
# matI <- mxMatrix( type="Iden", nrow=nv, ncol=nv, name="I" )
# invSD <- mxAlgebra( expression=solve(sqrt(I*V)), name="iSD" )
# 
# # Calculate genetic and environmental correlations
# corA <- mxAlgebra( expression=solve(sqrt(I*VA))%&%VA, name ="rA" ) #cov2cor()
# corD <- mxAlgebra( expression=solve(sqrt(I*VD))%&%VD, name ="rD" )
# corE <- mxAlgebra( expression=solve(sqrt(I*VE))%&%VE, name ="rE" )
# 
# # Create Algebra for Unstandardized and Standardized Variance Components
# rowUS <- rep('US',nv)
# colUS <- rep(c('VA','VD','VE','SA','SD','SE'),each=nv)
# estUS <- mxAlgebra( expression=cbind(VA,VD,VE,VA/V,VD/V,VE/V), name="US", dimnames=list(rowUS,colUS) )
# 
# # Create Confidence Interval Objects
# odd <- seq(1+3*nv,2*3*nv,nv)
# even <- seq(2+3*nv,2*3*nv,nv)
# ciADE <- mxCI( c("US[1,odd]","US[2,odd]","US[2,even]") )
# 
# # Build Model with Confidence Intervals
# calc <- list( matI, invSD, corA, corD, corE, estUS, ciADE )
# modelADE <- mxModel( "vcu_twoADE2", pars, modelMZ, modelDZ, multi, calc )
# 
# # ----------------------------------------------------------------------------------------------------------------------
# # RUN MODEL
# 
# # Run ADE Model
# vcu_fitADE <- mxRun( modelADE, intervals=T )
# sumADE <- summary( vcu_fitADE )
# 
# # Compare with Saturated Model
# mxCompare( vcu_fitSAT, vcu_fitADE )
# 
# # Print Goodness-of-fit Statistics & Parameter Estimates
# fitGofs(vcu_fitADE)
# fitEstCis(vcu_fitADE)
# 
# mxCompare( vcu_fitSAT, nested <- list( vcu_fitACE, vcu_fitADE) )