################################
##  Genomic Selection: GBLUP  ##
##    Multi-trait Model       ##
################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/TestData/SalmonGIll/GBLUP")
library(asreml)
library(ASRgenomics)

# Reading Phenotypic data
data(pheno.salmon)
datag <- pheno.salmon
head(datag)

# Reading AHAT inverse (sparse format) 
load(file='../GAmat.RData')

# Checking Ginverse
head(Ginv.blend)
head(attr(Ginv.blend, "rowNames"))
head(attr(Ginv.blend, "colNames"))
attr(Ginv.blend, "INVERSE")

# Matching data
check <- match.kinship2pheno(K=G_blend, pheno.data=datag,
                             indiv='ID', clean=FALSE, mism=TRUE)

# Some EDA
head(datag)
hist(datag$mean_gill_score)
plot(datag$mean_gill_score, datag$amoebic_load)
cor(datag$mean_gill_score, datag$amoebic_load)

########################
# Fitting GBLUP Model in ASReml-R

# Defining factors
head(datag)
datag$ID <- as.factor(datag$ID)
datag <- datag[order(datag$ID),]
str(datag)

## GBLUP ###
# y = mu + vm(INDIV,Ginv) + e

# Fitting gill_score
modelGBLUP1 <- asreml(fixed=mean_gill_score~1,
                    random=~vm(ID, Ginv.blend),
                    residual=~idv(units),
                    workspace=128e06,
                    na.action=na.method(y="include"),
                    data=datag)
plot(modelGBLUP1)
summary(modelGBLUP1)$varcomp
vpredict(modelGBLUP1, h2~V1/(V1+V2))

# Fitting amoebic_load
modelGBLUP2 <- asreml(fixed=amoebic_load~1,
                     random=~vm(ID, Ginv.blend),
                     residual=~idv(units),
                     workspace=128e06,
                     na.action=na.method(y="include"),
                     data=datag)
plot(modelGBLUP2)
summary(modelGBLUP2)$varcomp
vpredict(modelGBLUP2, h2~V1/(V1+V2))

# Fitting Bivariate Model
initG <- c(-0.8, 0.18, 2.60)
initR <- c(-0.8, 0.56, 7.91)
modelBV <- asreml(fixed=cbind(mean_gill_score,amoebic_load)~trait,
                      random=~vm(ID, Ginv.blend):corgh(trait, init=initG),
                      residual=~id(units):corgh(trait, init=initR),
                      workspace=128e06,
                      na.action=na.method(y="include"),
                      data=datag)
plot(modelBV)
summary(modelBV)$varcomp

## Obtaining GEBVs and Predictions - BLUP
BLUP <- as.data.frame(summary(modelBV, coef=TRUE)$coef.random)
head(BLUP)
