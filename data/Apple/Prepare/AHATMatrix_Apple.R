#####################################
##  Genomic Selection              ##
##  Obtaining AHAT and its inverse ##
##                                 ##
##    Apple                        ##
#####################################

rm(list=ls()) 
setwd("D:\\OneDrive - CGIAR\\Data Analysis\\genomic_seelction_training_vsni\\data\\Apple\\Prepare")
library(ASRgenomics)

# Apple Dataset
? pheno.apple
? geno.apple

##############################
# Obtaining APED

library(nadiv)
pedigree<-read.table("Pedigree_Apple.txt",h=TRUE) 
head(pedigree)
tail(pedigree)

inv <- ainverse(pedigree[c(1:3)])

APED<-makeA(pedigree[,1:3])
APED[1:5,1:5]
APED<-as.matrix(APED[48:294,48:294])  # Only offspring
APED[1:5,1:5]

ASRgenomics::kinship.heatmap(APED)

##############################
# Obtaining M (under QC)

# Reading full marker data
data(geno.apple)
markerFULL <- geno.apple
markerFULL[1:5,1:5]
dim(markerFULL)

M_filter <- qc.filtering(M=markerFULL, base=FALSE, ref=NULL, maf=0.05, marker.callrate=0.2,
                         ind.callrate=0.20, heterozygosity = 0.7, Fis = 0.5,
                         impute=FALSE, na.string='-9', plots=TRUE)
M_filter$M.clean[1:5,1:5]
dim(M_filter$M.clean) # 247 x 2732

M_filter$plot.heteroz
M_filter$plot.maf
M_filter$plot.missing.ind
M_filter$plot.missing.SNP
M_filter$plot.heteroz
M_filter$plot.Fis

# Saving imputed, filtered M matrix
write.table(M_filter$M.clean, file='data/Apple/Prepare/Mclean.txt')

##############################
# Obtaining AHAT

G <- G.matrix(M=M_filter$M.clean, method='VanRaden', na.string=NA)$G
G[1:5, 1:5]

#################################
# Diagnostics on the Kinship Matrix

check_G <- kinship.diagnostics(K=G)
ls(check_G)
check_G$list.diagonal
check_G$plot.diag
check_G$plot.offdiag  
check_G$list.duplicate

#################################
## Checking Pedigree vs Genomic-base Matrices

G2A <- match.G2A(A=APED, G=G, clean=TRUE, ord=TRUE, mism=TRUE, RMdiff=TRUE)

dim(G2A$Aclean)
dim(G2A$Gclean)
G2A$Aclean[1:8, 1:8]
G2A$Gclean[1:8, 1:8]

G2A$plotG2A
head(G2A$RM)

boxplot(G2A$RM$absdiff)
G2A$RM[G2A$RM$absdiff > 0.75,]

#################################
## Tunning a Genomic-based Kinship Matrix

G_blend <- G.tuneup(G=G, blend=TRUE, pblend=0.05)$Gb
G_bend <- G.tuneup(G=G, bend=TRUE)$Gb

G[1:8, 1:8]
G_blend[1:8, 1:8]
G_bend[1:8, 1:8]

check_G_blend <- kinship.diagnostics(K=G_blend)
check_G_blend$plot.diag
check_G_blend$plot.offdiag

check_G_bend <- kinship.diagnostics(K=G_bend)
check_G_bend$plot.diag
check_G_bend$plot.offdiag

#################################
## Verification with the G-inverse

Ginv.blend <- G.inverse(G=G_blend, sparseform=TRUE)$Ginv
Ginv.bend <- G.inverse(G=G_bend, sparseform=TRUE)$Ginv

head(Ginv.bend)
head(attr(Ginv.bend, "rowNames"))
head(attr(Ginv.bend, "colNames"))
attr(Ginv.bend, "INVERSE")

# Plotting G_bend Matrix
kinship.heatmap(G_bend)

# Saving Matrices
save(G_bend, Ginv.bend, file = "../Gmat.RData")

##################################
## Generating a Dominance Kinship Matrix

GD <- G.matrix(M=M_filter$M.clean, method='Su', na.string=NA)$G
GD[1:5, 1:5]
check_GD <- kinship.diagnostics(K=GD)

GD_bend <- G.tuneup(G=GD, bend=TRUE)$Gb
GDinv <- G.inverse(G=GD_bend, sparseform=TRUE)$Ginv
head(GDinv)

save(GD_bend, GDinv, file = "../GDmat.RData")
