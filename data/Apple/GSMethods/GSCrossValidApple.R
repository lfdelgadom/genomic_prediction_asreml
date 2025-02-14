###############################
## Genomic Selection
## k-fold cross-validation
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
markerFULL[1:5,1:5]
dim(markerFULL)  # 247 x 2732
# NOTE: order needs to match 
sum(rownames(markerFULL) == datag$INDIV)

##############################

# Fitting a Model - Bayes B
#nIter <- 15000; burnIn <- 5000  # Minimum recommended
nIter <- 15000; burnIn <- 5000    # For Testing
X <- as.matrix(markerFULL)
y <- datag$WT_HB

#ETA <- list(list(X=X, model="BayesB", df0=5, S0=500))
ETA <- list(list(X=X, model="BayesB"))

#############################
# k-fold cross-validation

k <- 5
n <- length(datag$INDIV)
sel <- rep(1:k, length.out = n) 
group <- sample(sel, n)
table(group)  

ypred_cv <- matrix(data=NA, nrow=n, ncol=1)
h2_cv <- matrix(data=NA, nrow=k, ncol=1)

for (g in 1:k)
{ 
  # Reading Phenotypic data
  ycv <- y
  
  for (j in 1:n) {
    if(group[j] == g) { ycv[j] <- NA } 
  }
  
  gsBayesBcv <- BGLR(y=ycv, ETA=ETA, nIter=nIter, burnIn=burnIn)
  predGScv <- gsBayesBcv$yHat

  (Vary <- var(y, na.rm=TRUE))
  (VarE <- gsBayesBcv$varE)
  (h2_cv[g] <- 1 - VarE/Vary)
  
  for (j in 1:n) {
    if(group[j] == g) { ypred_cv[j] <- predGScv[j] } 
  }
}

(corr_cv <- cor(y, ypred_cv, method='pearson', use="complete.obs"))
plot(y, ypred_cv)
h2_cv
mean(h2_cv)

# Predictive Ability  corr(yadj,ghat)
(PA <- corr_cv)

# Predictive Accuracy corr(greal,ghat)
(ACC <- PA/sqrt(mean(h2_cv)))  # Approximation
