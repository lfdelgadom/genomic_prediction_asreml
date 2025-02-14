###############################
## Sire Model (Female) 
## Data Open Polinization
###############################

rm(list=ls()) 
# setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/Phenotypic/Tree")
library(asreml)

# Open Pollinated Example
openp <- read.table("data/Phenotypic/Tree/OPENPOL.txt", h=T)
head(openp,20)
openp$REP <- as.factor(openp$REP)
openp$FEMALE <- as.factor(openp$FEMALE)
openp$PLOT <- as.factor(openp$PLOT)
openp$TYPE <- as.factor(openp$TYPE)
openp$INDIV <- as.factor(openp$INDIV)
str(openp)

table(openp$INDIV, openp$PLOT)

# Fitting Individual model
# y = mu + REP + REP:PLOT + INDIV + e 
mi <- asreml(fixed=HT~REP,
             random=~INDIV+REP:PLOT,
             #residual=~ar1v(row):ar1(col),
             na.action=na.method(y="include",x="include"),
             data=openp) # model did not converge bc we only have one obs per indv

#################
# Use residuals
# res = y - (mu + REP + REP:PLOT + INDIV)
# res* + INDIV = y - (mu + REP + REP:PLOT) 

mres <- asreml(fixed=HT~-1+REP,
             random=~REP:PLOT,
             na.action=na.method(y="include",x="include"),
             data=openp)
summary(mres)$varcomp
wald.asreml(mres)
BLUE <- as.data.frame(summary(mres, coef=TRUE)$coef.fixed)
BLUE
mean(BLUE$solution)

res <- resid(mres)
openp$yadj_ht <- res + mean(BLUE$solution)
head(openp)
plot(openp$HT, openp$yadj_ht)
cor(openp$HT, openp$yadj_ht, use='complete.obs')

