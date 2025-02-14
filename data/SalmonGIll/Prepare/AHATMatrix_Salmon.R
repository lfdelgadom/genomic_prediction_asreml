#########################################
##  Genomic Selection                  ##
##  Obtaining AHAT and its inverse     ##
##                                     ##
##    Salmon - Gill Score              ##
#########################################

rm(list=ls()) 
setwd("C:/Users/sgeza/OneDrive/Desktop/Online_Training/GS_GWAS/Datasets/SalmonGIll/Prepare")
library(ASRgenomics)

# Reading full marker data
? pheno.salmon
? geno.salmon
data(geno.salmon)
markerFULL <- geno.salmon
markerFULL[1:5,1:5]
dim(markerFULL)  # 1481 x 17156

# Reading pedigree data
data(ped.salmon)
ped <- ped.salmon
head(ped)
tail(ped)
dim(ped)  # 1561 x 3

# Reading Phenotypic data
data(pheno.salmon)
datag <- pheno.salmon
head(datag)

##############################
# Obtaining M (under QC)

# Filtering Marker Data (not imputed)
M_filter <- qc.filtering(M=markerFULL, base=FALSE, ref=NULL, maf=0.05, marker.callrate=0.2,
                         ind.callrate=0.20, heterozygosity = 0.60, Fis = 0.75, 
                         impute=FALSE, na.string=NA, plots=TRUE)
M_filter$M.clean[1:5,1:5]
dim(M_filter$M.clean) # 1481 x 11146

M_filter$plot.heteroz
M_filter$plot.maf
M_filter$plot.missing.ind
M_filter$plot.missing.SNP
M_filter$plot.heteroz
M_filter$plot.Fis

# Obtaining M_filter IMPUTED
M_filter <- qc.filtering(M=markerFULL, base=FALSE, ref=NULL, maf=0.05, marker.callrate=0.2,
                         ind.callrate=0.20, heterozygosity = 0.60, Fis = 0.75, 
                         impute=TRUE, na.string=NA, plots=TRUE)
M_filter$M.clean[1:5,1:5]

# Saving imputed, filtered M matrix
write.table(M_filter$M.clean, file='../Mclean.txt')

##############################
# Obtaining APED

head(ped.salmon)
tail(ped.salmon)

APED <- AGHmatrix::Amatrix(ped[,1:3])
APED[1:5,1:5]
APED <- as.matrix(APED[81:1561,81:1561])  # Only offspring
APED[1:5,1:5]

ASRgenomics::kinship.heatmap(APED)

##############################
# Obtaining AHAT

# Methods: VanRaden, Yang
G <- G.matrix(M=round(M_filter$M.clean,0), method='VanRaden', na.string=NA)$G
#G <- G.matrix(M=M_filter$M.clean, method='VanRaden', na.string=NA)$G
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
# Checking Pedigree vs Genomic-based Matrices

G2A <- match.G2A(A=APED, G=G, clean=TRUE, ord=TRUE, mism=TRUE, RMdiff=TRUE)

dim(G2A$Aclean)
dim(G2A$Gclean)
G2A$Aclean[1:8, 1:8]
G2A$Gclean[1:8, 1:8]

G2A$plotG2A
head(G2A$RM)
boxplot(G2A$RM$absdiff)
G2A$RM[G2A$RM$absdiff > 0.6,]

#################################
# Tunning a Genomic-based Kinship Matrix

# Trying Blending
G_blend <- G.tuneup(G=G, blend=TRUE, pblend=0.05)$Gb

check_G_blend <- kinship.diagnostics(K=G_blend)
check_G_blend$plot.diag
check_G_blend$plot.offdiag

# Trying Bending
G_bend <- G.tuneup(G=G, bend=TRUE)$Gb

G[1:8, 1:8]
G_blend[1:8, 1:8]
G_bend[1:8, 1:8]

check_G_bend <- kinship.diagnostics(K=G_bend)
check_G_bend$plot.diag
check_G_bend$plot.offdiag

#################################
# Verification with the G-inverse

# Using G_blend
Ginv.blend <- G.inverse(G=G_blend, sparseform=TRUE)$Ginv

head(Ginv.blend)
head(attr(Ginv.blend, "rowNames"))
head(attr(Ginv.blend, "colNames"))
attr(Ginv.blend, "INVERSE")

# Plotting G_bend Matrix
kinship.heatmap(G_blend)

# Saving Matrices
save(G_blend, Ginv.blend, file = "../GAmat.RData")

##################################
# Generating a Dominance Kinship Matrix

# Methods: Su, Vitezica
GD <- G.matrix(M=round(M_filter$M.clean,0), method='Su', na.string=NA)$G
GD[1:5, 1:5]

check_GD <- kinship.diagnostics(K=GD)
check_G_bend$plot.diag
check_G_bend$plot.offdiag

#GD_bend <- G.tuneup(G=GD, bend=TRUE)$Gb
GDinv <- G.inverse(G=GD, sparseform=TRUE)$Ginv
head(GDinv)

save(GD, GDinv, file = "../GDmat.RData")
