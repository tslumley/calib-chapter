#
# R Commands used with NWTS data to simulate a single case-cohort study and perform two-phase
# stratified analyses using the three methods described in the text
#
# Set working directory
#
setwd("C:\\IEA08")
#
# Load previously installed packages 
#
library(Hmisc)
library(survival)
library(survey)
library(Epi)
#
# Load NWTS data
#
nwt <- read.table("nwts-share.txt",header=T)
#
# Create additional variables used in various analyses
#
nwt$idx<-1:length(nwt$age)     # ID is sequence number
nwt$age1<-pmin(nwt$age,1)      # First age variable in linear spline
nwt$age2<-(nwt$age-1)*as.numeric(nwt$age>1)     # Second age variable in linear spline
nwt$UH<-nwt$histol             # Unfavorable histology (binary indicator)
nwt$stg34<-as.numeric(nwt$stage>2)    # Stage III or IV (binary indicator)
nwt$UHage1<-nwt$UH*nwt$age1    # Interaction between UH and first age variable
nwt$UHage2<-nwt$UH*nwt$age2    # Interaction between UH and second age variable
nwt$stgdiam<-nwt$stg34*nwt$tumdiam     # Interaction beteween tumor diameter and Stage III/IV
#
# Calculate stratum (1-8) based on age, stage and institutional histology
#
nwt$strt<-1+4*nwt$instit+2*nwt$stg34+as.numeric(nwt$age>1)
#
# Add a ninth stratum containing all (relapsed) cases
#
nwt$strt[nwt$relaps==1]<-9      # Strata 1-8 of non-relapsed controls 
#
# Create vector with true Phase II sampling fractions for 9 strata
#
p.true<-c(120/452,160/1620,1,120/914,1,1,1,1,1)
nwt$p.true<-p.true[nwt$strt]
#
# Fit model to main cohort (full) data with robust standard error
#
full<-coxph(Surv(trel,relaps)~UH*(age1+age2)+stg34+tumdiam+stgdiam,data=nwt,robust=T)
summary(full)
#
# Randomly select subjects for phase two sample
#
idx.case<-nwt$idx[nwt$relaps==1]   # Vector of ID's for relapsed cases
idx.cont0<-nwt$idx[nwt$relaps==0&(nwt$strt==3|nwt$instit==1)]   # Vector of ID's for controls sampled at 100% 
idx.cont1<-nwt$idx[nwt$relaps==0&nwt$strt==1]    # Vector of ID's for controls in stratum 1
idx.cont2<-nwt$idx[nwt$relaps==0&nwt$strt==2]    # Vector of ID's for controls in stratum 2
idx.cont4<-nwt$idx[nwt$relaps==0&nwt$strt==4]    # Vector of ID's for controls in stratum 4
idx.scont1<-sample(idx.cont1,size=120)           # Sample of ID's of size 120 from controls in stratum 1
idx.scont2<-sample(idx.cont2,size=160)           # Sample of ID's of size 160 from controls in stratum 2
idx.scont4<-sample(idx.cont4,size=120)           # Sample of ID's of size 120 from controls in stratum 4
idx.sample<-sort(cbind(idx.case,idx.cont0,idx.scont1,idx.scont2,idx.scont4)) # ID's for subjects sampled at Phase II
#
# Create indicator of whether or not sampled at phase two
#
nwt$ins<-F
nwt$ins[idx.sample]<-T
#
# Construct dataset containing Phase II observations only
#
nwtII<-nwt[nwt$ins==T,]
#
# Create inverse probability weights based on known sampling fractions
#
nwtII$wt<-1/nwtII$p.true
#
# Naive analysis of Phase II data using a priori weights
#
naive<-coxph(Surv(trel,relaps)~UH*(age1+age2)+stg34+tumdiam+stgdiam,
weights=wt,data=nwtII,robust=T)
sum.naive<-summary(naive)
round(sum.naive$coef,4)

# Define standard two phase stratified sampling design based on simple
# random sampling at phase one and stratified sampling at phase two
# (See documentation for survey package)
#
dstrat<-twophase(id=list(~1,~1),strata=list(NULL,~strt),subset=~ins,data=nwt)
#
# Prediction model (weighted logistic regression) for unfavorable histology according to central lab
# based on unfavorable histology according to institution, age, stage and study as described in text
#
Hmodel<-svyglm(histol~instit*I(stage>3)+I(age>10)+factor(study),family=quasibinomial,design=dstrat)
#
# Impute values of unfavorable histology (central) for all subjects in main cohort using prediction model
# just fit (Hmodel). Variable estH contains the imputed (predicted) values.
#
nwt$estH<-predict(Hmodel,type="response",newdata=nwt,se=F)
#
# Fit Cox regression model of interest to main cohort using imputed values for unfavorable histology
# (central) and actual values for other covariates. This is the "calibration model".
#
calmodel<-coxph(Surv(trel,relaps)~estH*(age1+age2)+stg34+tumdiam+stgdiam,data=nwt)
sumcal<-summary(calmodel)
round(sumcal$coef,4)
#
# Extract influence function contributions (delta-betas) (plus one for numerical reasons) from the calibration
# model for use as auxiliary variables in adjusting the weights in Horwitz-Thompson analysis of phase II data
#
db<-resid(calmodel,"dfbeta")+1
colnames(db)<-paste("db",1:ncol(db),sep="")
#
# Add the influence function contributions to main cohort dataset
#
nwtDB<-cbind(nwt,db)
#
# Redefine standard two phase stratified sampling design so that it contains the newly constructed variablesw
#
dstrt<-twophase(id=list(~1,~1),strata=list(NULL,~strt),subset=~ins,data=nwtDB)
#
# Calibrate the design so that the grand total of subjects (3915) and the main cohort totals of the
# augmented influence function contributions (which should all equal 1) are exactly estimated
#
dcal<-calibrate(dstrt,formula=make.formula(colnames(db)),pop=c(`(Intercept)`=3915,colSums(db)),calfun="raking",eps=0.0001)
#
# Construct design that uses estimated weights where estimates are based on stratum levels only.
# This should yield results identical to those based on the standard Horwitz-Thompson analysis
# without adustment of sampling weights, i.e. to estimates based on dstrt.
#
drrz0<-estWeights(nwt,strata=~strt,subset=~ins,formula=~factor(strt))
#
# Construct design that uses estimated weights where the auxiliary variables include both the
# sampling stratum indicators and also the influence function contributions
#
drrz<-estWeights(nwtDB,strata=~strt,subset=~ins,formula=update(make.formula(colnames(db)),~.+factor(strt)))
#
# Fit model to two-phase stratified data using the standard weights.
# This corresponds to Estimator II of Borgan et al. (LDA, 2000)
#
borgan<-svycoxph(Surv(trel,relaps)~UH*(age1+age2)+factor(stg34)*tumdiam,design=dstrt)
sum.borgan<-summary(borgan)
round(sum.borgan$coef,4)
round(sqrt(diag(attr(borgan$var,"phases")$phase1)),3)
round(sqrt(diag(attr(borgan$var,"phases")$phase2)),3)
# 
# Fit model to two-phase stratified data using estimated weights where estimates are
# determined solely from sampling strata. This should yield results identical to those
# from the standard Horwitz-Thompson (Borgan) approach. 
#
rrz0<-svycoxph(Surv(trel,relaps)~UH*(age1+age2)+factor(stg34)*tumdiam,design=drrz0)
# 
# Fit model using weights calibrated to the estimated influence function contributions
#
cal<-svycoxph(Surv(trel,relaps)~UH*(age1+age2)+factor(stg34)*tumdiam,design=dcal)
sum.cal<-summary(cal)
round(sum.cal$coef,4)
round(sqrt(diag(attr(cal$var,"phases")$phase1)),3)
round(sqrt(diag(attr(cal$var,"phases")$phase2)),3)
#
# Fit model using weights estimated from sampling strata and estimated influence function contributions
#
rrz<-svycoxph(Surv(trel,relaps)~UH*(age1+age2)+factor(stg34)*tumdiam,design=drrz)
sum.rrz<-summary(rrz)
round(sum.rrz$coef,4)
round(sqrt(diag(attr(rrz$var,"phases")$phase1)),3)
round(sqrt(diag(attr(rrz$var,"phases")$phase2)),3)
#
#



