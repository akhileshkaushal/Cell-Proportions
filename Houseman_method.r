####################################################################################################################################
#Houseman et al. [5] developed a method for cell type correction that capitalizes on the idea that differentially 
#methylated regions (DMRs) can serve as a signature for the distribution of different types of white blood cells. 
#It uses these DMRs as a surrogate in a regression calibration based technique to identify the cell mixture distribution. 
#Regression calibration technique can lead to bias estimate, thus external validation data is used to calibrate the model 
#and to correct for the bias [6]. Validation data set used by Houseman et al. method consists of DNA-methylation data and 
#sorted cell types of 46 white blood cell samples from Infinium HumanMethylation27 Beadchip  and AllCells Ⓡ , LLC (Emeryville, CA) respectively. 
#Their method was specifically for the Illumina 27k beadchip array.
###################################################################################################################################

#(a) Visit this website
  https://drive.google.com/folderview?id=0B6IN3-9RV-LBeUt5bENMSXlsY2s&usp= sharing
#to download all the relevant files

#(b) Specify the input data in the R code by the name “Rcodes_Cell_mixture.R” on the line that reads following:

beta=t(as.matrix(read.table(".../input_data.txt",sep="\t",header=T)))

#(c) Run the above R code by using

source("Rcodes_Cell_mixture.R")

#inside R console (make sure the working directory contains all of the files downloaded in step (a).
#(d) This program will create a file by the name “Ewasher_Houseman_cell.csv”, 
#which will have cell proportions for the 204 subjects contained in the example data set obtained above.
#(e) The estimated cell proportions will be used as adjusting co-variate in a linear model described below:

library(limma);
library(MASS);
cell = read.csv("Ewasher_data_housecell.csv",header=T);
cell1 = cell[,2];cell2 = cell[,3];cell3 = cell[,4];cell4 = cell[,5];
cell5 = cell[,6];cell6 = cell[,7]

mod11 = model.matrix(~x1+x2+x3+x4+x5+x6+cell2+cell3+cell4+cell5+cell6);
fit1  = lmFit(y12,mod11,method="robust");
fite1 = eBayes(fit1);
tab1  = topTable(fite1, coef = "x1",number=length(y1[,1]), p.val=0.05,adjust = "fdr");
write.table(tab1,file="Ewasher_data_houseman_CpG.csv",sep=",");

#The file Ewasher_data_houseman_CpG.csv will have a list of CpG sites whose association with primary covariate x1 
#(input_phenotye.txt) is adjusted for underlying cell proportions.
