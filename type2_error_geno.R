#setwd("R:/Basic_Sciences/Pharm/Perera_Lab/Yizhen/human_panel/human_panel/aa_0927/")
#setwd("/Volumes/fsmresfiles/Basic_Sciences/Pharm/Perera_Lab/Yizhen/human_panel/human_panel/aa_0927/")
setwd("/projects/b1047/zhong/aa_0927")
ancestry = read.table("data/hgcr_0927_81_qc_chr10.removedup.snp.la.input",header = TRUE, stringsAsFactors = FALSE)
snp_fn = "data/hgcr_0927_81_qc_chr22.dosage"

library(LAMatrix)
library(pROC)
args = commandArgs(trailingOnly=TRUE)
print("index, snp_fn, out_prefix, model")
index = as.numeric(args[1])
snp_fn = args[2]
out_prefix = args[3]
model = args[4]


snp = read.table(snp_fn, header = T,stringsAsFactors = F)
use_proc <- function(out, gt){
  res = list()
  out$eqtl = apply(out,1,function(x) paste0(x[1],"_",x[2]))
  label = rep(0,nrow(out))
  label[which(out$eqtl %in% gt)] = 1
  out$p = -log10(out$pvalue)
  f_roc = roc(label,out$p,plot=F)
  p_roc = roc(label,out$p,plot=F,algorithm = 2,partial.auc=c(100,80),partial.acu.foucs="sens",percent=T)
  res[["f_roc"]] = f_roc
  res[["p_roc"]] = p_roc
  return(res)
}

gt = sapply(seq(1,50,1),function(x) paste0(x,"_",x))

#cl <- makeCluster(4) # create a cluster with 2 cores
#registerDoParallel(cl) # register the cluster

#res = foreach(i = 1:2, 
#              .combine = "list",
#.packages = c("LAMatrix", "pROC")) %dopar% 
cvrt = SlicedData$new();
cvrt$fileDelimiter = " ";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 0;       # one column of row labels
cvrt$LoadFile("data/hgcr_0927_covariate.txt");


cvrt1 = SlicedData$new();
cvrt1$fileDelimiter = " ";      # the TAB character
cvrt1$fileOmitCharacters = "NA"; # denote missing values;
cvrt1$fileSkipRows = 1;          # one row of column labels
cvrt1$fileSkipColumns = 0;       # one column of row labels
cvrt1$LoadFile("data/hgcr_0927_covariate_nopc.txt");


set.seed(index)
for(i in index:(index+5000)){
  print(i)
  snp_id = sample(1:nrow(snp), 1000, replace = FALSE, prob = NULL)
  snp_mat = snp[snp_id,-1] #remove the SNP_ID in the first column
  rownames(snp_mat) = c(1:nrow(snp_mat))
  
  local_mat = ancestry[snp_id,-1] #remove the SNP_ID in the first column
  rownames(local_mat) = c(1:nrow(local_mat))
  
  snp_local = cbind(snp_mat[1:50,],local_mat[1:50,]) 

  if(model == "both"){
        gene_exp = apply(snp_local,1, function(x) rnorm(81,mean=0,sd =1)+
                  as.numeric(x[1:81])*0.9 + as.numeric(x[82:length(x)])*0.8) #gene expression with both genotype and local ancestry effect
  }else if(model == "geno"){
        gene_exp = apply(snp_mat[1:50,],1, function(x) rnorm(81,mean=0,sd =1)+as.numeric(x)*0.9)
  }else{
        print("chooose model from both and geno")
  }
  
  gene_exp_2 = sapply(1:450, function(x) rnorm(81,mean=0,sd =1))
  gene_exp = cbind(gene_exp,gene_exp_2)
  gene_exp = t(gene_exp)
  rownames(gene_exp) = c(1:nrow(gene_exp))
  
 
  
  gene = SlicedData$new();
  gene$CreateFromMatrix( gene_exp ) 
  
  snps = SlicedData$new();
  snps$CreateFromMatrix(as.matrix(snp_mat))


  local_mat = ancestry[snp_id,-1]
  rownames(local_mat) = c(1:nrow(local_mat))
  local = SlicedData$new();
  local$CreateFromMatrix(as.matrix(local_mat))
  
  
  errorCovariance = numeric();
  output_file_name = tempfile();
  
  
  
  me = LAMatrix_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt1,
    output_file_name = output_file_name,
    pvOutputThreshold = 1,
    useModel = modelLINEAR, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  me_cov = LAMatrix_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = 1,
    useModel = modelLINEAR, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  me_local = LAMatrix_main(
    snps = snps,
    gene = gene,
    local = local,
    cvrt = cvrt1,
    output_file_name = output_file_name,
    pvOutputThreshold = 1,
    useModel = modelLOCAL, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE);
  
  
  out1 = me$all$eqtls
  out2 = me_cov$all$eqtls
  out3 = me_local$all$eqtls
  
  out1 = out1[order(out1[,1], out1[,2]),]
  out2 = out2[order(out2[,1], out2[,2]),]
  out3 = out3[order(out3[,1], out3[,2]),]
  
  res = list()
  roc1= use_proc(out1, gt)
  roc2 = use_proc(out2, gt)
  roc3 = use_proc(out3, gt) #record the roc
  res[["roc1"]] = roc1
  res[["roc2"]] = roc2
  res[["roc3"]] = roc3
  
  #roc_p1p2 = roc.test(roc1[["f_roc"]], roc2[["f_roc"]], method = "delong", paired = T)
  #roc_p2p3 = roc.test(roc2[["f_roc"]], roc3[["f_roc"]], method = "delong", paired = T)
  #roc_p1p3 = roc.test(roc1[["f_roc"]], roc3[["f_roc"]], method = "delong", paired = T)
 
  
  #proc_p1p2 = roc.test(roc1[["p_roc"]], roc2[["p_roc"]], method = "bootstrap", paired = T, boot.n = 2)
  #proc_p2p3 = roc.test(roc2[["p_roc"]], roc3[["p_roc"]], method = "bootstrap", paired = T, boot.n = 2)
  #proc_p1p3 = roc.test(roc1[["p_roc"]], roc3[["p_roc"]], method = "bootstrap", paired = T,boot.n = 2)

  ps = c(as.numeric(roc1[["f_roc"]]$auc),as.numeric(roc2[["f_roc"]]$auc), as.numeric(roc3[["f_roc"]]$auc), as.numeric(roc1[["p_roc"]]$auc),as.numeric(roc2[["p_roc"]]$auc), as.numeric(roc3[["p_roc"]]$auc))
  #ps = c(roc_p1p2$p.value,roc_p2p3$p.value,roc_p1p3$p.value, 
  #       as.numeric(roc1[["f_roc"]]$auc),as.numeric(roc2[["f_roc"]]$auc), as.numeric(roc3[["f_roc"]]$auc),
  #       proc_p1p2$p.value, proc_p2p3$p.value, proc_p1p3$p.value,
  #       as.numeric(roc1[["p_roc"]]$auc),as.numeric(roc2[["p_roc"]]$auc), as.numeric(roc3[["p_roc"]]$auc))
  print(ps)
  saveRDS(res, paste0(out_prefix,"_",i,".rds"))
  saveRDS(ps, paste0(out_prefix,"_p_",i,".rds"))
}








