library(LAMatrix)


# ################################################
# #test a table
# ################################################

base.dir = find.package('LAMatrix');
# Genotype file name
SNP_file_name = paste(base.dir, "/data/SNP.txt", sep="");

#local file name
LOCAL_file_name = paste(base.dir, "/data/local.txt", sep="");

#Gene expression file name
expression_file_name = paste(base.dir, "/data/GE.txt", sep="");

#Covariates file name
# # Set to character() for no covariates
covariates_file_name = paste(base.dir, "/data/Covariates.txt", sep="");
#
# # Output file name
output_file_name = tempfile();
#
# # Only associations significant at this level will be saved
pvOutputThreshold = 1e-2;
#
# # Error covariance matrix
# # Set to numeric() for identity.
errorCovariance = numeric();
# # errorCovariance = read.table("Sample_Data/errorCovariance.txt");
#
#
# ## Load genotype data
#
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);
#
# ## Load gene expression data
#
genes = SlicedData$new();
genes$fileDelimiter = "\t";      # the TAB character
genes$fileOmitCharacters = "NA"; # denote missing values;
genes$fileSkipRows = 1;          # one row of column labels
genes$fileSkipColumns = 1;       # one column of row labels
genes$fileSliceSize = 2000;      # read file in slices of 2,000 rows
genes$LoadFile(expression_file_name);
#
# ## Load covariates
#
cvrts = SlicedData$new();
cvrts$fileDelimiter = "\t";      # the TAB character
cvrts$fileOmitCharacters = "NA"; # denote missing values;
cvrts$fileSkipRows = 1;          # one row of column labels
cvrts$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrts$LoadFile(covariates_file_name);
}
#
# ## Load local data
#
locals = SlicedData$new();
locals$fileDelimiter = "\t";      # the TAB character
locals$fileOmitCharacters = "NA"; # denote missing values;
locals$fileSkipRows = 1;          # one row of column labels
locals$fileSkipColumns = 1;       # one column of row labels
locals$fileSliceSize = 2000;      # read file in slices of 2,000 rows
locals$LoadFile(LOCAL_file_name);
#

## Run the analysis

me2 = LAMatrix_main(
  snps = snps,
  gene = genes,
  cvrt = cvrts,
  local = locals,
  shuffle = FALSE,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = modelLOCAL,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)

snp_con = read.table(SNP_file_name,header = TRUE,stringsAsFactors = FALSE)
snp_con = snp_con[,c(2:ncol(snp_con))]
gene_con = read.table(expression_file_name,header = TRUE,stringsAsFactors = FALSE)
gene_con = gene_con[,c(2:ncol(gene_con))]
local_con = read.table(LOCAL_file_name,header = TRUE)
local_con = local_con[,c(2:ncol(local_con))]
cvrt_con = read.table(covariates_file_name,header = TRUE)
cvrt_con = cvrt_con[,c(2:ncol(cvrt_con))]
table = NULL
for(i in 1:nrow(snp_con)){
  for(j in 1:nrow(gene_con)){
    temp = summary(lm(as.numeric(gene_con[j,]) ~ as.numeric(snp_con[i,])+as.numeric(local_con[i,])+as.numeric(cvrt_con[1,])+as.numeric(cvrt_con[2,])))$coefficients[2,]
    temp_table = matrix(c(i,j,as.numeric(temp)),nrow=1,ncol=6)
    table = rbind(table,temp_table)
  }
}

table = as.data.frame(table)
colnames(table) = c("snp","gene","estimate","std.error","t.value","p")
table = table[order(table$p),]

cat("\n\n Matrix eQTL: \n");
print(show(me$all$eqtls));
cat("\n R summary output: \n");
print(table[which(table$p <= pvOutputThreshold),])

