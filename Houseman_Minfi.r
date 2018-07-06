###############################################################################################################################
#The method by Jaffe and Irizarry [7] was adapted from the Houseman et al. [5] method and is tailored for 
#Illumina450k along with 27k array. The algorithm in Houseman et al. identifies 500 CpG sites used to estimate cell 
#mixture proportions from the Illumina 27k array. The modification of Jaffe and Irizarry was motivated because of 
#the existence of probe SNPs in the 500 CpG sites and the inconsistency of CpG sites between the 27k and 450k arrays. 
#In addition, the flow-sorted data of the six adult male subjects were used as references [8] .
###############################################################################################################################

#(a) Install the required packages in R using:

source("https://bioconductor.org/biocLite.R")
biocLite(c("minfi","quadprog","FlowSorted.Blood.450k",
"IlluminaHumanMethylation450kmanifest",
"IlluminaHumanMethylation450kanno.ilmn12.hg19"));

#(b) Below is the R code to obtain cell estimates for the example data set:


lib = c("minfi","quadprog","FlowSorted.Blood.450k",
"IlluminaHumanMethylation450kmanifest",
"IlluminaHumanMethylation450kanno.ilmn12.hg19");
lapply(lib, require, character.only = TRUE);

grSet1=read.table("input_data.txt",header=T); grSet=grSet1[,-1];
rownames(grSet)=grSet1[,1]; grSet=data.matrix(grSet);
referenceMset = get('FlowSorted.Blood.450k.compTable');

cell = c("CD8T","CD4T", "NK","Bcell","Mono","Gran","Eos");
compData = minfi:::pickCompProbes(referenceMset, cellTypes=cell);
coefs = compData$coefEsts;
coefs = coefs[ intersect(rownames(grSet), rownames(coefs)), ];
rm(referenceMset);

counts = minfi:::projectCellType(grSet[rownames(coefs), ], coefs);
rownames(counts) = colnames(grSet);
write.table(counts,file="Ewasher_data_minficell.csv",sep=",");


###############################################################################################################################
#The file by the name Ewasher_data_minficell.csv will contain the cell proportions estimated using the approach extended 
#from the Houseman et al. method and implemented in the R package minfi.
###############################################################################################################################

#(c) The estimated cell proportions will be used as adjusting co-variate in a linear model described below:

cell = read.csv("Ewasher_data_minficell.csv",header=T);

cell1 = cell[,1];cell2 = cell[,3];cell3 = cell[,4];cell5 = cell[,6];
cell6 = cell[,7];cell7 = cell[,8];

mod11 = model.matrix(~x1+x2+x3+x4+x5+x6+cell2+cell3+cell4+cell5+cell6+cell7);
fit1  = lmFit(y12,mod11,method="robust");
fite1 = eBayes(fit1);
tab1  = topTable(fite1, coef = "x1",number=length(y1[,1]), p.val=0.05,adjust = "fdr");
write.table(tab1,file="Ewasher_data_minfi_CpG.csv",sep=",");

###############################################################################################################################
#The file Ewasher_data_minfi_CpG.csv will have a list of CpG sites whose association with primary covariate x1 
#(input_phenotye.txt) is adjusted for underlying cell proportions.
###############################################################################################################################
