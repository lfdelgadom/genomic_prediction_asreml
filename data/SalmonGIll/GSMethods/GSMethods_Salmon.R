###############################
## Genomic Selection: BGLR
###############################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/SalmonGIll/GSRR")
library(ASRgenomics)
library(BGLR)

# Reading Phenotypic data
data(pheno.salmon)
datag <- pheno.salmon
head(datag)

# Reading full marker data
markerFULL <- read.table('../Mclean.txt', header=TRUE)
markerFULL[1:5,1:5]
dim(markerFULL)  # 1481 x 11146
# NOTE: order does not match 
markerFULL <- markerFULL[order(rownames(markerFULL)),]
markerFULL[1:5,1:5]
sum(rownames(markerFULL) == datag$ID)

##############################
# Fitting a Model - Bayes B

#nIter <- 15000; burnIn <- 5000  # Minimum recommended
nIter <- 6000; burnIn <- 2000    # For Testing
#X2 <- model.matrix(datag$Block)  # Incidence matrix for factor Block!
X <- as.matrix(markerFULL)
y <- datag$mean_gill_score

ETA <- list(list(X=X, model="BayesB"))
gsBayesB <- BGLR(y=y, ETA=ETA, nIter=nIter, burnIn=burnIn, saveAt='TTest_')
summary(gsBayesB)
ls(gsBayesB)

# Total Predictions
predGS <- gsBayesB$yHat
head(predGS)
plot(y, predGS)
(corr_pearson <- cor(y, predGS, method='pearson', use="complete.obs"))

# Getting bHAT Parameters 
bHAT <- gsBayesB$ETA[[1]]$b
head(bHAT)
boxplot(bHAT)
hist(bHAT)

############################################################
## Several Bayesian Methods

## Bayesian Ridge Regression (Gaussian prior), "equivalent" to G-BLUP
ETA <- list(list(X=X, model='BRR'))
gsBRR <- BGLR(y=y, ETA=ETA, nIter=nIter, burnIn=burnIn)
summary(gsBRR)

## Bayes A (Scaled-t prior)
ETA <- list(list(X=X, model='BayesA'))
gsBayesA <- BGLR(y=y, ETA=ETA, nIter=nIter)
summary(gsBayesA)

## Bayes B
ETA <- list(list(X=X, model="BayesB", df0=5))
gsBayesB <- BGLR(y=y, ETA=ETA, nIter=nIter, burnIn=burnIn)
summary(gsBayesB)

## Bayesian LASSO (Laplace prior)
ETA <- list(list(X=X, model='BL'))
gsBL <- BGLR(y=y, ETA=ETA, nIter=nIter, burnIn=burnIn)
summary(gsBL)

## Bayes C (point of mass at zero + Gaussian slab)
ETA <- list(list(X=X, model='BayesC'))
gsBayesC <- BGLR(y=y, ETA=ETA, nIter=nIter, burnIn=burnIn)
summary(gsBayesC)

#######################################################
# Calculating goodness-of-fit statistics

# Heritability (GS)
(Vary <- var(y,na.rm=TRUE))
(VarE <- gsBayesB$varE)
(h2_GS <- 1 - VarE/Vary)

# Predictive Ability  corr(yadj,ghat)
predGS <- gsBayesB$yHat
(PA <- cor(y, predGS, method='pearson', use="complete.obs"))

# Predictive Accuracy corr(greal,ghat)
# Need a vector with the greal
(ACC <- PA/sqrt(h2_GS))  # Approximation
