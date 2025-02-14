###############################
## Genomic Selection: BGLR
###############################

rm(list=ls()) 
setwd("D:/OneDrive - CGIAR/Data Analysis/genomic_seelction_training_vsni/data/Apple/GSRR")
library(ASRgenomics)
library(BGLR)

# Reading Phenotypic data
data(pheno.apple)
datag <- pheno.apple
head(datag)

# Reading full marker data
markerFULL <- read.table('../Mclean.txt', header=TRUE)
markerFULL[1:5,1:10]
dim(markerFULL)  # 247 x 2732
# NOTE: order needs to match 
sum(rownames(markerFULL) == datag$INDIV)

##############################

# Fitting a Model - Bayes B
#nIter <- 15000; burnIn <- 5000  # Minimum recommended
nIter <- 6000; burnIn <- 700    # For Testing
#X2 <- model.matrix(datagFULL$Block)  # Incidence matrix for factor Block!
X <- as.matrix(markerFULL)
y <- datag$WT_HB

ETA <- list(list(X=X, model="BayesB", df0=5, S0=500))
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

# other way to obtain yHat
# yHAT <- mean(datag$WT_HB, na.rm=TRUE) + X %*% bHAT
# head(yHAT)
# cor(yHAT, y, use='complete.obs')

# Predicted Individuals (not measured)
dataP <- cbind(datag[,1:8], predGS)
(dataP[gsBayesB$whichNa,])

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
(Vary <- var(y, na.rm=TRUE))
(VarE <- gsBayesB$varE)
(h2_GS <- 1 - VarE/Vary)

# Predictive Ability  corr(yadj,ghat)
predGS <- gsBayesB$yHat
(PA <- cor(y, predGS, method='pearson', use="complete.obs"))

# Predictive Accuracy corr(greal,ghat)
# Need a vector with the greal
(ACC <- PA/sqrt(h2_GS))  # Approximation
