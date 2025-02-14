##################################
## Spatial Analysis - WHEAT     ##
##################################

rm(list=ls()) 

library(asreml)
# setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/Phenotypic/Spat")

# Reading data
spatial <- read.table(file='data/Phenotypic/Spat/wheat.txt', header=TRUE)
head(spatial)

# Generating factors 
spatial$rep <- as.factor(spatial$rep)
spatial$ibk <- as.factor(spatial$ibk)
spatial$check <- as.factor(spatial$check)
spatial$gen <- as.factor(spatial$gen)
spatial$col <- as.factor(spatial$col)  
spatial$row <- as.factor(spatial$row)  
str(spatial)

# Some plotting
library(desplot)
desplot(spatial,yield~row+col)
desplot(spatial,yield~row+col, out1=ibk)
desplot(spatial,yield~row+col, out1=rep)
head(table(spatial$gen,spatial$rep))

#########################
# Spatial - WITH many things

spatial1 <- asreml(fixed=yield~rep + at(check,'1'):gen,
                  random= ~ rep:ibk + at(check,'0'):gen,
                  residual=~ar1(col):ar1v(row),
                  data=spatial)
spatial1 <- update.asreml(spatial1)
plot(varioGram(spatial1))
plot(spatial1)
summary(spatial1)$varcomp # variance componentes

wald.asreml(spatial1,denDF='numeric')

(H2<-vpredict(spatial1,H2~V2/(V1+V2+V6)))   # Traditional definition

#########################
# Getting Predictions

spatialp <- asreml(fixed=yield~rep + gen,
                   random= ~ rep:ibk,
                   residual=~ar1(col):ar1v(row),
                   data=spatial)
summary(spatialp)$varcomp

preds <- predict.asreml(spatialp, classify='gen')$pvals
head(preds)

# Comparing with means
mgen <- aggregate(yield~gen, FUN=mean, data=spatial)
head(mgen)

comb <- cbind(preds,mgen)
head(comb)
plot(comb$predicted.value, comb$yield)
cor(comb$predicted.value, comb$yield)

