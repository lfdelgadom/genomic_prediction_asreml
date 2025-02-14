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
datag$INDIV <- as.factor(datag$INDIV)
str(datag)

## GBLUP ###
# y = mu + INDIV + e

modelGBLUP <- asreml(fixed = FF_HB ~1,
                  random = ~vm(INDIV,Ginv.bend),
                  residual = ~idv(units),
                  workspace = 128e06,
                  na.action = na.method(y = "include"),
                  data = datag)
plot(modelGBLUP)
summary(modelGBLUP)$varcomp
vpredict(modelGBLUP, h2~V1/(V1+V2))

##################
# Bivariate Model

# [y1 y2] = trait + INDIV:trait + e:trait
modelBV <- asreml(fixed= cbind(JUI_HB, FF_HB) ~trait,
                 random = ~vm(INDIV, Ginv.bend):corgh(trait),
                 residual = ~id(units):corgh(trait),
                 workspace = 128e06,
                 na.action = na.method(y = "include"),
                 data = datag)
plot(modelBV)
summary(modelBV)$varcomp

vpredict(modelBV, h2_1~V2/(V2+V6))
vpredict(modelBV, h2_2~V3/(V3+V7))

###################
# MET

# Reading Phenotypic data
datal <- read.table(file = 'phenotAppleMET.txt', header=TRUE)
head(datal)

# Defining factors
head(datal)
datal$INDIV <- as.factor(datal$INDIV)
datal$SITE <- as.factor(datal$SITE)
str(datal)

# y = mu + SITE + INDIV:SITE + e
modelMET <- asreml(fixed = CRI ~SITE,
                 random = ~vm(INDIV, Ginv.bend):corgh(SITE),
                 residual = ~dsum(~units|SITE),
                 workspace = 128e06,
                 na.action = na.method(y = "include"),
                 data = datal)
modelMET <- update.asreml(modelMET)
summary(modelMET)$varcomp
wald(modelMET, denDF="numeric")
plot(modelMET)

vpredict(modelMET, h2_HB~V2/(V2+V4))
vpredict(modelMET, h2_MOT~V3/(V3+V5))

BLUP <- summary(modelMET, coef=TRUE)$coef.random
View(BLUP)

preds <- predict.asreml(modelMET, classify='INDIV:SITE')$pvals
head(preds)
predsA <- predict.asreml(modelMET, classify='INDIV')$pvals
head(predsA)
