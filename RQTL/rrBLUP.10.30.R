# Simplest Version 11.3.14
require(rrBLUP)

setwd("/Users/JKC/Documents/M.S. Project/Thesis Material/Thesis/RQTL")

pheno<-as.data.frame(read.table("Pheno.10.30.txt", header=T, sep="\t"))

geno<-as.data.frame(read.table("Geno.10.30.txt",header=T,sep="\t"))

GWAS(pheno,geno,P3D=T)
#end

Markers<-as.matrix(read.table("markers.10.30.txt", header=F))

head(Markers)

pheno<-as.data.frame(read.table("Pheno.10.30.txt", header=T, sep="\t"))

geno<-as.data.frame(read.table("Geno.10.30.txt",header=T,sep="\t"))

attach(geno)

head(pheno)

dim(Markers)
dim(pheno)

impute=A.mat(Markers,max.missing=0.5, impute.method="mean", return.imputed=T)

Markers_impute=impute$imputed

impute$imputed

write.table(Markers_impute, "Markers_impute.10.30.txt", sep="\t", row.names=F, col.names=F)

Markers_impute<-as.matrix(read.table("Markers_impute.10.30.txt",header=T, sep="\t"))

GWAS(pheno,geno,K=Markers_impute)
par(mfrow=c(1,2))
gwas.10.30<-GWAS(pheno,geno, P3D=F)
gwas.10.31<-GWAS(pheno,geno)
GWAS(pheno,geno,P3D=T)
write.table(gwas.10.31,"GWAS.10.31.txt", sep="\t")

write.table(geno, "Geno.txt", sep="\t")

require(knitr)
install.packages("qqman")
library(qqman)

x<-read.table("GWAS.10.31.txt", header=T, sep="\t")

str(x)

manhattan(x, ymax = 6)

as.data.frame(table(x$CHR))

qq(x$P)

# 11.2 with all mapped SNPs

require(rrBLUP)

setwd("/Users/JKC/Documents/M.S. Project/Thesis Material/Thesis/RQTL")

Markers1<-as.matrix(read.table("markers.11.2.txt", header=F))

pheno<-as.data.frame(read.table("Pheno.10.30.txt", header=T, sep="\t"))

impute=A.mat(Markers1,max.missing=0.5, impute.method="mean", return.imputed=T)

Markers1_impute=impute$imputed

impute$imputed

write.table(Markers1_impute, "Markers_impute.11.2.txt", sep="\t", row.names=F, col.names=F)

