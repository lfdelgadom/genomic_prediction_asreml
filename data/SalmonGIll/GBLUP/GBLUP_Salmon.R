################################
##  Genomic Selection: GBLUP  ##
##  Library ASReml-R
################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/SalmonGIll/GBLUP")
library(asreml)
library(ASRgenomics)

# Reading Phenotypic data
data(pheno.salmon)
datag <- pheno.salmon
head(datag)

# Reading AHAT inverse (sparse format) 
load(file = '../GAmat.RData')

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
summary(datag$mean_gill_score)

########################
# Fitting GBLUP Model in ASReml-R

# Defining factors
head(datag)
datag$ID <- as.factor(datag$ID)
datag <- datag[order(datag$ID),]
str(datag)

## GBLUP ###
# y = mu + vm(INDIV,Ginv) + e

modelGBLUP <- asreml(fixed = mean_gill_score ~ 1,
                    random = ~vm(ID, Ginv.blend),
                    residual = ~idv(units),
                    workspace = 128e06,
                    na.action = na.method(y = "include"),
                    data = datag)
plot(modelGBLUP)
summary(modelGBLUP)$varcomp
vpredict(modelGBLUP, h2~V1/(V1+V2))

## Obtaining GEBVs and Predictions - BLUP
BLUP <- as.data.frame(summary(modelGBLUP, coef=TRUE)$coef.random)
BLUP <- BLUP[order(rownames(BLUP)),]
head(BLUP)
vc <- summary(modelGBLUP)$varcomp
BLUP$reliab <- 1 - BLUP$std.error^2/vc[1,1]  # Reliab = 1 - PEV/s2a
head(BLUP)
View(BLUP)

# Saving BLUPs (or GEBVs)
#write.table(BLUP, file='../BLUPGS.txt')

# Yhat, Ypred
predGBLUP <- predict(modelGBLUP, classify = "ID")$pvals
predGBLUP <- as.data.frame(predGBLUP)
predGBLUP <- predGBLUP[order(as.character(predGBLUP$ID)),]
head(predGBLUP)

(corr_pearson <- cor(datag$mean_gill_score, predGBLUP[,2], method='pearson', use="complete.obs"))
plot(datag$mean_gill_score, predGBLUP[,2])

# Saving predictions
# write.table(predGBLUP, file='../PredGS.txt', sep = " ", quote = F)

#######################################################
# Calculating goodness-of-fit statistics

# Heritability (GBLUP)
summary(modelGBLUP)$varcomp
(h2_GBLUP <- vpredict(modelGBLUP, h2~V1/(V1+V2)))

# Predictive Ability  corr(yadj,ghat)
(PA <- cor(datag$mean_gill_score, predGBLUP[,2], method='pearson', use="complete.obs"))

# Predictive Accuracy corr(greal,ghat)
# Need a vector with the greal
(AC <- PA/sqrt(as.numeric(h2_GBLUP[1])))


#######################################
# Linking to get u's - IT TAKES TIME

BLUP <- as.data.frame(summary(modelGBLUP, coef=TRUE)$coef.random)
#BLUP <- BLUP[order(rownames(BLUP)),]
head(BLUP)  # GEBV

# a = M*u
# M'*a = M'*M*u / inv(M'*M)
# u     =  inv(M'*M) * M' * a
#  p*1  =     p*p          p*n    n*1  

# Reading full marker data
markerFULL <- read.table(file='../Mclean.txt', header=TRUE)
markerFULL[1:5,1:5]
dim(markerFULL) # 1481 x 11146

library(MASS)
? ginv
a <- as.matrix(BLUP$solution)
M <- as.matrix(markerFULL) 
CM <- t(M) %*% M
invCM <- ginv(CM)  # Using ginverse - it takes a long time!
#invCM[1:5,1:5]
u <- invCM %*% t(M) %*% a
head(u)
boxplot(u)
dim(u)

# Predictions (new individuals, using Mnew!)
yhatGP <- mean(datag$mean_gill_score, na.rm=TRUE) + M %*% u
head(yhatGP)

plot(yhatGP, predGBLUP$predicted.value)
cor(yhatGP, predGBLUP$predicted.value)

#save(u, file='uGBLUP.RData')