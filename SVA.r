#######################################################################################################################
#Surrogate variable analysis (SVA) estimates potential confounding factors based on singular value decompositions 
#(SVD) of residuals.
#######################################################################################################################

#(a) Install the bioconductor package “sva“, “limma” and CRAN package “MASS” in R using: 

source("https://bioconductor.org/biocLite.R");
biocLite(c("sva","limma"));
install.packages("MASS");

#(b) Use the following R code to implement SVA

lib = c("MASS","sva","limma");
lapply(lib, require, character.only = TRUE);
mod1<-model.matrix(~x1+x2+x3+x4+x5+x6)
mod01<-model.matrix(~x2+x3+x4+x5+x6)
svobj1= sva(edata,mod1,mod01,n.sv=NULL,method="two-step");
modSv1 = cbind(mod1,svobj1$sv);
fit1 = lmFit(edata,modSv1,method="robust");
fite1 = eBayes(fit1);
tab1 = topTable(fite1, coef = "x1",number=length(y1[,1]), p.val=0.05,adjust = "fdr");
write.table(tab1,file="...\filename.csv",sep=",")

################################################################################################################################
#filname.csv will have list of CpG found to be significantly associated with the primary covariate 
#(contained in file “input_phenotype.txt”) after adjusting for cellular heterogeneity using SVA.
################################################################################################################################
