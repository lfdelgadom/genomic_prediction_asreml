#################################
## Single-step GBLUP           ##
##   Pine genomic Data         ##
#################################
  
rm(list=ls())
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/ssGBLUP")

library(ASRgenomics) 
library(asreml)      

######################################################
#  Genotypic dataset (655 genotypes) - Loblolly pine

? geno.pine655
? pheno.pine

# Phenotypic data # 861
data(pheno.pine)
head(pheno.pine)
tail(pheno.pine)
dim(pheno.pine)

# Pedigree data # 2034
data(ped.pine)
head(ped.pine)
tail(ped.pine)
dim(ped.pine)

# Molecular data # 655
data(geno.pine655)
geno.pine655[1:5, 1:5]
dim(geno.pine655)

######################################################
# Reading and Filtering a Molecular Dataset

geno.pine655[1:5, 1:5]
dim(geno.pine655)
? ASRgenomics::qc.filtering

M_filter <- qc.filtering(M = geno.pine655, base = FALSE, ref = NULL,
                         marker.callrate = 0.20, ind.callrate = 0.20,
                         maf = 0.05, heterozygosity = 0.95, Fis = 1,
                         impute = FALSE, na.string = '-9', plots = TRUE)
M_filter$M.clean[1:5, 1:5]
dim(M_filter$M.clean)
M_filter$plot.missing.ind
M_filter$plot.missing.SNP
M_filter$plot.maf
M_filter$plot.heteroz
M_filter$plot.Fis

######################################################
# Generating a Kinship Matrix

# Pedigree-based relationship matrix (A)
head(ped.pine)
A <- AGHmatrix::Amatrix(data = ped.pine)
A[601:605, 601:605]

# Genomic-based relationship matrix (GRM)
? G.matrix
G <- G.matrix(M = M_filter$M.clean, method = 'VanRaden', na.string = NA)$G
dim(G)
G[1:5, 1:5]

######################################################
# Diagnostics on the Kinship Matrix

check_G <- kinship.diagnostics(K = G,  diagonal.thr.small = 0.8,
                               diagonal.thr.large = 1.2, duplicate.thr = 0.95)
check_G$list.diagonal
check_G$list.duplicate
check_G$plot.diag
check_G$plot.offdiag

######################################################
# Tuning a Genomic-based Kinship Matrix

? G.tuneup
# Options: blend, bend, align

# Matching GRM to A matrix
dim(A)
dim(G)
G2A <- match.G2A(A = A, G = G, clean = TRUE, ord = TRUE, mism = TRUE, RMdiff = TRUE)
ls(G2A)
dim(G2A$Aclean)
dim(G2A$Gclean)
G2A$Aclean[343:347, 343:347]
G2A$Gclean[343:347, 343:347]

# Some comparison of GRM to A values
G2A$plotG2A
head(G2A$RM)
head(G2A$RM[G2A$RM$absdiff > 0.20,])

# Aligning the G matrix to the A matrix
G_align <- G.tuneup(G = G2A$Gclean, A = G2A$Aclean, align = TRUE)$Gb
G2A$Gclean[343:347, 343:347]
G_align[343:347, 343:347]

# Checking again RM to see if there still issues
Ga2A <- match.G2A(A = A, G = G_align, clean = TRUE, ord = TRUE, mism = TRUE,  RMdiff = TRUE)
Ga2A$plotG2A
dim(Ga2A$RM[Ga2A$RM$absdiff > 0.20,])
head(Ga2A$RM[Ga2A$RM$absdiff > 0.25,])

# Diagnostics of the align matrix
check_G$plot.diag      # original
check_G$plot.offdiag   # original

check_G_align <- kinship.diagnostics(K = G_align,  diagonal.thr.small = 0.8,
                               diagonal.thr.large = 1.2, duplicate.thr = 0.95)
check_G_align$plot.diag     # aligned
check_G_align$plot.offdiag  # aligned

######################################################
# Preparing to fit a Genomic-BLUP (GBLUP) model with ASReml-R

head(pheno.pine)
dim(pheno.pine)
dim(G_align)

pheno.G <- match.kinship2pheno(K = G_align, pheno.data = pheno.pine,
                               indiv = 'Genotype', clean = FALSE, mism = TRUE)
head(pheno.G$matchesP)
pheno.subset <- pheno.pine[pheno.G$matchesP, ]
dim(pheno.subset)
head(pheno.subset)

# Obtaining inverse sparse matrix for ASReml-R
? G.inverse
Ginv.sparse <- G.inverse(G = G_align, sparseform = TRUE)$Ginv
head(Ginv.sparse)
head(attr(Ginv.sparse, "rowNames"))
head(attr(Ginv.sparse, "colNames"))
attr(Ginv.sparse, "INVERSE")


######################################################
######################################################
# Fitting a Genomic-BLUP (GBLUP) model with ASReml-R

library(asreml)
dim(pheno.subset)

pheno.subset$Genotype <- as.factor(pheno.subset$Genotype)
head(pheno.subset)
str(pheno.subset)

# y = mu + Genotype + e
GBLUP <- asreml(fixed = DBH_Adj ~ 1,
                random = ~vm(Genotype, Ginv.sparse),
                residual = ~idv(units),
                na.action = na.method(y = 'include'),
                workspace = 1e07,
                data = pheno.subset)
summary(GBLUP)$varcomp
plot(GBLUP)

# Heritability and other output (BLUP, GEBV)
vpredict(GBLUP, h2~V1/(V1+V2))
GEBV <- summary(GBLUP, coef = TRUE)$coef.random # Genomic breeding values
head(GEBV)
#preds <- predict.asreml(GBLUP, classify = 'Genotype')$pvals  # Predictions
#head(preds)


######################################################
######################################################
# Generating/Exploring a Hybrid Genomic Matrix (H)

? H.inverse

Ginv <- G.inverse(G = G_align, sparseform = FALSE)$Ginv
Ginv[1:5, 1:5] # Note we are using G_align ...
Hinv.sparse <- H.inverse(A = A, G = Ginv, lambda = 0.90, sparseform = TRUE)
head(Hinv.sparse)

H <- H.matrix(A = A, Ginv = Ginv, lambda = 0.90, sparseform = FALSE)
A[25:30, 25:30]
H[25:30, 25:30]
A[44:49, 44:49]
H[44:49, 44:49]
A[605:610, 605:610]
H[605:610, 605:610]
dim(H)

check_H <- kinship.diagnostics(K = H,  diagonal.thr.small = 0.8,
                               diagonal.thr.large = 1.2, duplicate.thr = 0.95,
                               sample.plot = 0.2)
check_H$plot.diag
check_H$plot.offdiag

######################################################
# Fitting a Single-Step Genomic-BLUP model (ssGBLUP) with ASReml-R

Hinv.sparse <- H.inverse(A = A, G = Ginv, lambda = 0.90, sparseform = TRUE)
pheno.pine$Genotype <- as.factor(pheno.pine$Genotype)
dim(pheno.pine)

ssGBLUP <- asreml(fixed = DBH_Adj ~ 1,
                  random = ~vm(Genotype, Hinv.sparse),
                  residual = ~idv(units),
                  na.action = na.method(y = 'include'),
                  workspace = 1e07,
                  data = pheno.pine)
plot(ssGBLUP)
summary(ssGBLUP)$varcomp
vpredict(ssGBLUP, h2~V1/(V1+V2))

HGEBV <- summary(ssGBLUP, coef = TRUE)$coef.random
head(HGEBV)
head(GEBV)
HGEBV <- data.frame(HGEBV)
GEBV <- data.frame(GEBV)

# Some comparison between HGEBV and GEBV
GEBV$gen <- substr(rownames(GEBV), 27, 40) 
HGEBV$gen <- substr(rownames(HGEBV), 27, 40) 

TT <- merge.data.frame(GEBV, HGEBV, by='gen')
head(TT) # In general, we see a reduction of std.error
dim(TT)
summary(TT)

plot(TT$solution.x, TT$solution.y)
cor(TT$solution.x, TT$solution.y)

mean(TT$solution.x)
mean(TT$solution.y)
