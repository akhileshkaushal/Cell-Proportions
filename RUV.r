################################################################################################################################
#The method of removing unwanted variation (RUV) uses information from reference database, but it does not estimate 
#cell type proportions. Instead, this approach bases on the information of negative control probes and performs factor 
#analysis on these probes to identify factors due to unmeasured confounders. These factors are then included in subsequent 
#analyses to adjust for cell type effects.
################################################################################################################################

#(a) Top 600 CpG sites associated with the blood cell types for 27k example data was obtained from the link: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-86 under the heading Electronic supplementary material. It can also be downloaded from google drive link mentioned under the Houseman’s method, the name of file is “WBC-Analysis-110319.RData”.

#(b) Code below identifies principal components related to the dna-methylation of 600 CpG sites.

library(data.table)
library(psych)
library(pracma)

beta<-read.csv("Top_600_Ewasher_27k.csv",header=T)
beta1=beta[,-1]
rownames(beta1)=beta[,1]
beta1_t<-t(beta1)

pca_mval<-prcomp(beta1_t,center=T,scale.=T)
plot(pca_mval, type = "l")

###Based upon scree plot we can choose number of PC's##
rawLoadings <- pca_mval$rotation[,1:7] %*% diag(pca_mval$sdev, 7, 7)
rotatedLoadings <- varimax(rawLoadings,normalize = TRUE, eps = 1e-5)$loadings
invLoadings <- t(pracma::pinv(rotatedLoadings))
scores <- scale(beta1_t) %*% invLoadings

###x1,x2,x3,x4,x5 and x6 are primary and secondary covariates,
###PC1,PC2, PC3 and PC4 are principal components
##############################################################
mod<-model.matrix(~x1+x2+x3+x4+x5+x6+PC1+PC2+PC3+PC4)

fit<-lmFit(edata,mod,method="robust")
fite<-eBayes(fit)

tab <- topTable(fite, coef = "x1",number=length(edata[,1]), p.val=0.05,adjust = "fdr")
