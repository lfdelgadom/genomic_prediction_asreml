################################
##  Genomic Selection: GBLUP  ##
################################

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

## GBLUP ###
# y = mu + vm(INDIV,Ginv) + e

modelGBLUP <- asreml(fixed = WT_HB ~ 1,
                  random = ~vm(INDIV, Ginv.bend),
                  residual = ~idv(units),
                  workspace = 128e06,
                  na.action = na.method(y = "include"),
                  data = datag)
plot(modelGBLUP)
summary(modelGBLUP)$varcomp
vpredict(modelGBLUP, h2~V1/(V1+V2))

## Obtaining GEBVs and Predictions - BLUP
BLUP <- as.data.frame(summary(modelGBLUP, coef = TRUE)$coef.random)
head(BLUP)
vc <- summary(modelGBLUP)$varcomp
BLUP$reliab <- 1 - BLUP$std.error^2/vc[1,1]  # Reliab = 1 - PEV/s2a
head(BLUP)
View(BLUP)

# Yhat, Ypred 
predGBLUP <- predict(modelGBLUP, classify = "INDIV")$pvals
head(predGBLUP)
#View(predGBLUP)

preds <- as.matrix(predGBLUP[,2])
(corr_pearson <- cor(datag$WT_HB,preds, method='pearson', use="complete.obs"))
plot(datag$WT_HB, preds)

# Predicted Individuals (not measured)
dataP <- cbind(datag[,1:5], preds)
(dataP[which(is.na(datag$WT_HB)),])

# Saving predictions
# write.table(predGBLUP, file = "PredGS.txt")

#######################################################
# Calculating goodness-of-fit statistics

# Heritability (GBLUP)
summary(modelGBLUP)$varcomp
(h2_GBLUP <- vpredict(modelGBLUP, h2_ind~V1/(V1+V2)))

# Predictive Ability  corr(yadj,ghat)
(PA <- cor(datag$WT_HB, preds, method='pearson', use="complete.obs"))

# Predictive Accuracy corr(greal,ghat)
# Need a vector with the greal
(AC <- PA/sqrt(as.numeric(h2_GBLUP[1])))


###############################
# Extended GBLUP model - Family

## GBLUP ###
# y = mu + INDIV + Family + e

modelGBLUPf <- asreml(fixed = WT_HB~1,
                    random = ~vm(INDIV, Ginv.bend) + Family,
                    residual = ~idv(units),
                    workspace = 128e06,
                    na.action = na.method(y = "include"),
                    data = datag)
modelGBLUPf <- update.asreml(modelGBLUPf)
plot(modelGBLUPf)
summary(modelGBLUPf)$varcomp
vpredict(modelGBLUPf, h2~V2/(V1+V2+V3))
vpredict(modelGBLUPf, d2~4*V1/(V1+V2+V3))


#######################################
# Linking to get u's

BLUP <- as.data.frame(summary(modelGBLUP, coef=TRUE)$coef.random)
head(BLUP) # GEBVs

# a = M*u
# M'*a = M'*M*u / inv(M'*M)
# u     =  inv(M'*M) * M' * a
#  p*1  =     p*p          p*n    n*1  

# Reading full marker data
markerFULL <- read.table(file='../Mclean.txt', header=TRUE)
markerFULL[1:5,1:5]
dim(markerFULL) # 247 x 2732

library(MASS)
? ginv
a <- as.matrix(BLUP$solution)
M <- as.matrix(markerFULL) # 247 x 2828
CM <- t(M) %*% M
invCM <- ginv(CM)  # Using ginverse - it takes a long time!
#invCM[1:5,1:5]
u <- invCM %*% t(M) %*% a
head(u)
boxplot(u)
dim(u)

# Predictions (new individuals, using Mnew!)
yhatGP <- mean(datag$WT_HB, na.rm=TRUE) + M %*% u
head(yhatGP)

plot(yhatGP, predGBLUP$predicted.value)
cor(yhatGP, predGBLUP$predicted.value)

#save(u, file='uGBLUP.RData')