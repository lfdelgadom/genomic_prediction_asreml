#####################################
##  Genomic Selection: GBLUP       ##
##    Additive + Dominance         ##
#####################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/SalmonGIll/GBLUP")
library(ASRgenomics)
library(asreml)

# Reading Phenotypic data
data(pheno.salmon)
datag <- pheno.salmon
head(datag)

# Reading AHAT inverse (sparse format) 
load(file='../GAmat.RData')
load(file='../GDmat.RData')

# Checking Ginverse
head(Ginv.blend)
head(attr(Ginv.blend, "rowNames"))
head(attr(Ginv.blend, "colNames"))
attr(Ginv.blend, "INVERSE")

head(GDinv)
head(attr(GDinv, "rowNames"))
head(attr(GDinv, "colNames"))
attr(GDinv, "INVERSE")

# Defining factors
head(datag)
datag$IDA <- as.factor(datag$ID)
datag$IDD <- as.factor(datag$ID)
str(datag)

## GBLUP ###

# Only additive effects
modelA <- asreml(fixed = mean_gill_score ~1,
                 random = ~vm(IDA, Ginv.blend),
                 residual = ~idv(units),
                 workspace = 128e06,
                 na.action = na.method(y = "include"),
                 data = datag)
plot(modelA)
summary(modelA)$varcomp
(h2 <- vpredict(modelA, h2~V1/(V1+V2)))

## Obtaining Predictions - BLUP - 
BLUP <- summary(modelA, coef=TRUE)$coef.random
BLUP <- BLUP[order(rownames(BLUP)),]
a.effects <- BLUP[,1]
(corr_pearson_A <- cor(datag$mean_gill_score, a.effects, method='pearson', use="complete.obs"))

# Both additive and dominant effects
modelAD <- asreml(fixed = mean_gill_score ~1,
                random = ~vm(IDA, Ginv.blend) + vm(IDD, GDinv),
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
BLUP <- BLUP[order(rownames(BLUP)),]
a.effects <- BLUP[1:1481,1]
d.effects <- BLUP[1482:2962,1]
g.effects <- a.effects + d.effects

(corr_pearson_A  <- cor(datag$mean_gill_score, a.effects, method='pearson', use="complete.obs"))
(corr_pearson_AD <- cor(datag$mean_gill_score, g.effects, method='pearson', use="complete.obs"))
