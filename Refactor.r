################################################################################################################################
#ReFACTor, uses a variant form of principal component analysis (PCA) to adjust for the cell type effects. 
#Principal components are based on top 500 most informative CpG sites. The identified principal components can be used 
#in the model to adjust for the confounding effects of cell types.
################################################################################################################################

#(a) Download the Refactor source code from following link: http://www.cs.tau.ac.il/~heran/cozygene/software/refactor.html

#(b) We modified the code obtained from above to exclude the CpG sites with standard deviation in the lowest 5th percentile.


sd1=apply(y1, 1, sd, na.rm = F)
include = which(sd1>=quantile(sd1,0.05))
O = y1[include,]
edata1 = edata[include,]

cpgnames <- rownames(O)
for (site in 1:nrow(O))
 {
 model <- lm(O[site,] ~ x1+x2+x3+x4+x5+x6)
 O_adj[site,] = residuals(model)
 }
 O = O_adj
 print('Running a standard PCA...')
 pcs = prcomp(scale(t(O)));

 coeff = pcs$rotation
 score = pcs$x

 print('Compute a low rank approximation of input data and rank sites...')
 x = score[,1:7]%*%t(coeff[,1:7]);
 An = scale(t(O),center=T,scale=F)
 Bn = scale(x,center=T,scale=F)
 An = t(t(An)*(1/sqrt(apply(An^2,2,sum))))
 Bn = t(t(Bn)*(1/sqrt(apply(Bn^2,2,sum))))

# Find the distance of each site from its low rank approximation.
 distances = apply((An-Bn)^2,2,sum)^0.5 ;
 dsort = sort(distances,index.return=T);
 ranked_list = dsort$ix

 print('Compute ReFACTor components...')
 sites = ranked_list[1:500];
 pcs = prcomp(scale(t(O[sites,])));
 first_score <- score[,1:7];
 score = pcs$x


mod=model.matrix(~x1+x2+x3+x4+x5+x6+first_score)
fit1<-lmFit(edata1,mod,method="ls")
fite1<-eBayes(fit1)
tab1 <- topTable(fite1, coef = "x1",number=length(y1[,1]), p.val=0.05,adjust = "fdr")

write.table(tab1,file="Ewasher_data_Refactor_results.csv",sep=",",row.names=T)

