####################################
##  Genomic Selection             ##
##    K-fold cross-validation     ##
####################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/Apple/GBLUP")
library(asreml)
library(ASRgenomics)

# Reading Phenotypic data
data(pheno.apple)
datag <- pheno.apple
head(datag)

# Reading AHAT inverse (sparse format) 
load(file = '../Gmat.RData')

# Checking Ginverse
head(Ginv.bend)
head(attr(Ginv.bend, "rowNames"))
head(attr(Ginv.bend, "colNames"))
attr(Ginv.bend, "INVERSE")

# Matching data
check <- match.kinship2pheno(K=G_bend, pheno.data=datag,
                             indiv='INDIV', clean=FALSE, mism=TRUE)

########################
# Fitting GBLUP Model in ASReml-R

# Defining factors
head(datag)
datag$INDIV <-as.factor(datag$INDIV)
datag$Family <- as.factor(datag$Family) 
str(datag)

#############################
# k-fold cross-validation

k <- 5
n <- length(datag$INDIV)
sel <- rep(1:k, length.out = n) 
group <- sample(sel, n)
table(group)  

ypred_cv <- matrix(data=NA, nrow=n, ncol=1)
h2_cv <- matrix(data=NA, nrow=k, ncol=1)

y <- datag$WT_HB   # Response variable

for (g in 1:k)
{ 
  ycv <- y
  
  for (j in 1:n) {
    if(group[j] == g) { ycv[j]<-NA } 
  }

  modelGBLUPcv <- asreml(fixed = ycv ~ 1,
                     random = ~vm(INDIV, Ginv.bend),
                     residual = ~idv(units),
                     workspace = 128e06,
                     na.action = na.method(y = "include"),
                     data = datag) 
  #modelGBLUPcv <- update.asreml(modelGBLUPcv)
  #BLUPcv <- as.data.frame(summary(modelGBLUPcv, coef=TRUE)$coef.random)
  #BLUPcv <- BLUPcv[order(rownames(BLUPcv)),]
  #BLUPcv <- BLUPcv[,1]
  predGBLUPcv <- predict(modelGBLUPcv, classify="INDIV", sed=T)$pvals
  predGBLUPcv <- predGBLUPcv[,2]
  
  (h2_GBLUPcv <- vpredict(modelGBLUPcv, h2~V1/(V1+V2)))
  (h2_cv[g] <- as.numeric(h2_GBLUPcv[1]))
  
  for (j in 1:n) {
    #if(group[j] == g) { ypred_cv[j] <- BLUPcv[j] } 
    if(group[j] == g) { ypred_cv[j] <- predGBLUPcv[j] } 
  }
}

head(cbind(group,y, ypred_cv))
plot(y, ypred_cv)
(PA <- cor(y, ypred_cv, method='pearson', use="complete.obs"))

# Heritability (GS cross-validation)
h2_cv
mean(h2_cv)

# Predictive Ability  corr(yadj,ghat)
PA

# Predictive Accuracy corr(greal,ghat)
(ACC <- PA/sqrt(mean(h2_cv)))  # Approximation

