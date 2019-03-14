setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/human_panel/human_panel/aa_0927")

library(foreach)
library(doParallel)
library(LAMatrix)


args = commandArgs(trailingOnly=TRUE)
print("effect, model")
effect = as.numeric(args[1])
model = args[2]





pc = read.table("./data/hgcr_0927_covariate.txt",stringsAsFactors = FALSE,header = T)
pc1 = as.numeric(pc[2,])

gene_location_file_name="./data/gene_name_gene_position_v19.txt"
snps_Location_File_name="./data/hgcr_0927_81_qc_chr10.removedup.snp_chr10.snppos"
snpspos=read.table(snps_Location_File_name,header=T,stringsAsFactors = FALSE)
genepos=read.table(gene_location_file_name,header = F,stringsAsFactors = FALSE)
average_ancestry_by_gene = read.csv("./result/averaged_la_by_gene_1005.csv",header = T,stringsAsFactors = F) #average local ancestry by gene


p_range = seq(1e-6,1,length.out = 100)
pvalue = seq(0,6,length.out = 7)
pvalue2 = seq(0,6,length.out = 20)
pvalue = 10^(-pvalue)
pvalue2 = 10^(-pvalue2)
p=c(pvalue,pvalue2)


false_positive <- function(out, p){
  fp_r = sapply(p,function(x) length(which(out$pvalue<x))/nrow(out))
  return(fp_r)
}




snps = SlicedData$new();
snps$fileDelimiter = " ";      
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile("data/hgcr_0927_81_qc_chr22.dosage");

cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 0;       # one column of row labels
cvrt$LoadFile("data/hgcr_0927_covariate.txt");

cvrt1 = SlicedData$new();
cvrt1$fileDelimiter = " ";      
cvrt1$fileOmitCharacters = "NA"; # denote missing values;
cvrt1$fileSkipRows = 1;          # one row of column labels
cvrt1$fileSkipColumns = 0;       # one column of row labels
cvrt1$LoadFile("data/hgcr_0927_covariate_nopc.txt");

local = SlicedData$new();
local$fileDelimiter = " ";      
local$fileOmitCharacters = "NA"; # denote missing values;
local$fileSkipRows = 1;          # one row of column labels
local$fileSkipColumns = 1;       # one column of row labels
local$fileSliceSize = 2000;      # read file in slices of 2,000 rows
local$LoadFile("data/hgcr_0927_81_qc_chr10.removedup.snp.la.input");

cisDist=1e6

cl <- makeCluster(4) # create a cluster with 2 cores
registerDoParallel(cl) # register the cluster
res = foreach(i = 1:10, 
              .combine = "rbind",
              .packages = "LAMatrix") %dopar% {
              
                ########### 
                #no stratification
                if(model=="local"){
                  gene = average_ancestry_by_gene[sample(1:nrow(average_ancestry_by_gene) ,100, replace = FALSE, prob = NULL),]
                  gene_exp = apply(gene,1, function(x) rnorm(81,mean=0,sd =1)+as.numeric(x[-1]) * effect)
                }else{
                  gene_exp = sapply(1:100, function(x) rnorm(81,mean=0,sd = 1)+pc1 *effect)
                  gene = gene[sample(1:nrow(gene) ,100, replace = FALSE, prob = NULL)]
                }
                gene_exp = t(gene_exp)
                rownames(gene_exp) = gene[,1]
                colnames(gene_exp) = colnames(gene)[-1]
                gene_test = SlicedData$new( gene_exp );
                #gene_test$rowNameSlices[[1]] = gene$geneid[sample(1:nrow(gene) ,100, replace = FALSE, prob = NULL)]
              
                errorCovariance = numeric();
                output_file_name = tempfile();
                output_File_name_cis = tempfile()
                
                me = LAMatrix_main(
                  snps = snps,
                  gene = gene_test,
                  cvrt = cvrt1,
                  output_file_name = output_file_name,
                  output_file_name.cis = output_File_name_cis,
                  pvOutputThreshold = 0,
                  pvOutputThreshold.cis = 1,
                  useModel = modelLINEAR, 
                  errorCovariance = errorCovariance, 
                  verbose = TRUE,
                  snpspos = snpspos,
                  genepos = genepos,
                  cisDist = cisDist,
                  pvalue.hist = TRUE,
                  min.pv.by.genesnp = FALSE,
                  noFDRsaveMemory = FALSE);
            
                
                me_cov = LAMatrix_main(
                  snps = snps,
                  gene = gene_test,
                  cvrt = cvrt,
                  output_file_name = output_file_name,
                  output_file_name.cis = output_File_name_cis,
                  pvOutputThreshold = 0,
                  pvOutputThreshold.cis = 1,
                  useModel = modelLINEAR, 
                  errorCovariance = errorCovariance, 
                  verbose = TRUE,
                  snpspos = snpspos,
                  genepos = genepos,
                  cisDist = cisDist,
                  pvalue.hist = TRUE,
                  min.pv.by.genesnp = FALSE,
                  noFDRsaveMemory = FALSE);
                out = me$cis$eqtls
                
                me_local = LAMatrix_main(
                  snps = snps,
                  gene = gene_test,
                  cvrt = cvrt1,
                  local = local,
                  output_file_name = output_file_name,
                  output_file_name.cis = output_File_name_cis,
                  pvOutputThreshold = 0,
                  pvOutputThreshold.cis = 1,
                  useModel = modelLOCAL, 
                  errorCovariance = errorCovariance, 
                  verbose = TRUE,
                  snpspos = snpspos,
                  genepos = genepos,
                  cisDist = cisDist,
                  pvalue.hist = TRUE,
                  min.pv.by.genesnp = FALSE,
                  noFDRsaveMemory = FALSE);

              
                
                out = me$cis$eqtls
                out2 = me_cov$cis$eqtls
                out3 = me_local$cis$eqtls
                c(false_positive(out, p),false_positive(out2, p),false_positive(out3, p) )        
}
stopCluster(cl) # shut down the cluster

res = as.data.frame(res)
colnames(res) = c(sapply(p,function(x) paste0("no_",x)),sapply(p,function(x) paste0("pc_",x)),
  sapply(p,function(x) paste0("la_",x)))

write.csv(res,paste0("./result/res_simu_type1_",model,"_effect_cis_gene_beta",effect,".csv"),quote = F)
