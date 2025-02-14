#####################################
##  Genomic Selection: GBLUP       ##
##    Additive + Dominance         ##
#####################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/Apple/GBLUP")

# Reading Phenotypic data
data(pheno.apple)
datag <- pheno.apple
head(datag)

# Reading AHAT inverse (sparse format) 
load(file = '../Gmat.RData')
load(file = '../GDmat.RData')

# Checking Ginverse
head(Ginv.bend)
head(attr(Ginv.bend, "rowNames"))
head(attr(Ginv.bend, "colNames"))
attr(Ginv.bend, "INVERSE")

head(GDinv)
head(attr(GDinv, "rowNames"))
head(attr(GDinv, "colNames"))
attr(GDinv, "INVERSE")

# Defining factors
head(datag)
datag$INDIVA <- as.factor(datag$INDIV)
datag$INDIVD <- as.factor(datag$INDIV)
datag$Family <- as.factor(datag$Family)
str(datag)

## GBLUP ###

# Only additive effects
modelA <- asreml(fixed = JUI_MOT ~1,
                  random = ~vm(INDIVA, Ginv.bend),
                  residual = ~idv(units),
                  workspace = 128e06,
                  na.action = na.method(y = "include"),
                  data = datag)
plot(modelA)
summary(modelA)$varcomp
(h2<-vpredict(modelA, h2~V1/(V1+V2)))

## Obtaining Predictions - BLUP - 
BLUP <- summary(modelA, coef=TRUE)$coef.random
a.effects <- BLUP[1:247, 1]
(corr_pearson_A <- cor(datag$JUI_MOT, a.effects, method='pearson', use="complete.obs"))

# Both additive and dominant effects
modelAD <- asreml(fixed = JUI_MOT ~1,
                   random = ~vm(INDIVA, Ginv.bend) + vm(INDIVD, GDinv),
                   residual = ~idv(units),
                   workspace = 128e06,
                   na.action = na.method(y = "include"),
                   data = datag)
plot(modelAD)
summary(modelAD)$varcomp
(h2 <- vpredict(modelAD, h2~V1/(V1+V2+V3)))
(d2 <- vpredict(modelAD, d2~V2/(V1+V2+V3)))

lrt.asreml(modelA, modelAD, boundary=TRUE)

## Obtaining Predictions - BLUP 
BLUP <- summary(modelAD, coef=TRUE)$coef.random
a.effects <- BLUP[1:247, 1]
d.effects <- BLUP[248:494, 1]
g.effects <- a.effects + d.effects

(corr_pearson_A  <- cor(datag$JUI_MOT, a.effects, method='pearson', use="complete.obs"))
(corr_pearson_AD <- cor(datag$JUI_MOT, g.effects, method='pearson', use="complete.obs"))

