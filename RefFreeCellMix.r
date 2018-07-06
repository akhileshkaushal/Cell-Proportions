##################################################################################################################################
#RefFreeCellMix is a function available in R package RefFreeEWAS. It uses a variant of non-negative matrix 
#factorization to decompose the total methylation sites into CpG-specific methylation states for a pre-specified number 
#of cell types and subject-specific cell-type distributions.
##################################################################################################################################

#(a) Install R-package RefFreeEWAS version 2.0

#(b) Use the code below:

library(RefFreeEWAS)
cell<-RefFreeCellMix(y2,mu0=NULL,K=7,iters=5,Yfinal=NULL,verbose=TRUE)
mod1<-model.matrix(~x1+x2+x3+x4+x5+x6+cell$Omega)

fit1<-lmFit(y2,mod1,method="robust")
fite1<-eBayes(fit1)
tab1 <- topTable(fite1, coef = "x1",number=length(O[,1]), p.val=0.05,adjust = "fdr")

write.table(tab1,file="Ewasher_data_RefFreeCellMix_results.csv", sep=",",row.names=T)
