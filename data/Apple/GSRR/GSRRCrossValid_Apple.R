####################################
##  Genomic Selection             ##
##    Library: rrBLUP             ##
##    K-fold cross-validation     ##
####################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/Apple/GSRR")
library(ASRgenomics)
library(rrBLUP)

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

# Defining Marker Matrix for RRBLUP
Z <- as.matrix(markerFULL)-1   # -1, 0, 1
Z[1:10,1:10]

#############################
# k-fold cross-validation

k <- 5
n <- length(datag$ID)
sel <- rep(1:5, length.out = n) 
group <- sample(sel, n)
table(group) 

ypred_cv <- matrix(data=NA, nrow=n, ncol=1)
h2_cv <- matrix(data=NA, nrow=k, ncol=1)

y <- datag$WT_HB

for (g in 1:k)
{ 
  ycv <- y
  
  for (j in 1:n) {
    if(group[j] == g) { ycv[j] <- NA } 
  }
  
  gsRRcv <- mixed.solve(ycv, Z=Z, K=NULL, X=NULL, SE=FALSE, method='REML', return.Hinv=FALSE)
  beta0cv <- gsRRcv$beta
  predGSRRcv <- matrix(data=beta0cv, nrow=length(y), ncol=1) + Z %*% gsRRcv$u
  
  (Vary <- var(ycv, na.rm=TRUE))
  (VarE <- gsRRcv$Ve)
  (h2_cv[g] <- 1-VarE/Vary)
  
  for (j in 1:n) {
    if(group[j] == g) { ypred_cv[j] <- predGSRRcv[j] } 
  }
}

head(cbind(group, y, predGSRRcv))
(PA <- cor(y, ypred_cv, method='pearson', use="complete.obs"))
plot(y, ypred_cv)

# Heritability (GS cross-validation)
h2_cv
mean(h2_cv)

# Predictive Ability  corr(yadj,ghat)
PA

# Predictive Accuracy corr(greal,ghat)
(ACC <- PA/sqrt(mean(h2_cv)))  # Approximation
