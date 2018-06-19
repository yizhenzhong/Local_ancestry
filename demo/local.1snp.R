library(LAMatrix)



n = 200;# Number of columns (samples)


nc = 10;# Number ofs covariates



cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);# Generate the covariates
snps.mat = floor(runif(n, min = 0, max = 3)); # Generate the genotype
local.mat = floor(runif(n, min = 0, max = 3)); # Generate the local ancestry


gene.mat = cvrt.mat %*% rnorm(nc) + rnorm(n) + 0.5 * snps.mat + 1 + 0.3 *local.mat; # Generate the expression vector

snps = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
genes = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrts = SlicedData$new( t(cvrt.mat) );
locals = SlicedData$new( matrix(local.mat,nrow=1))

# Produce no output files
filename = NULL; # tempfile()

me = LAMatrix_main(
  snps = snps,
  gene = genes,
  cvrt = cvrts,
  local = locals,
  output_file_name = filename,
  pvOutputThreshold = 1,
  useModel = modelLOCAL,
  verbose = TRUE,
  pvalue.hist = FALSE );


# Pull Matrix eQTL results - t-statistic and p-value
beta = me$all$eqtls$beta;
tstat = me$all$eqtls$statistic;
pvalue = me$all$eqtls$pvalue;
rez = c(beta = beta, tstat = tstat, pvalue = pvalue);


# And compare to those from the linear regression in R
cat("\n\n Matrix eQTL: \n");
print(rez);
cat("\n R summary(lm()) output: \n");
lmdl = lm( gene.mat ~ snps.mat + cvrt.mat + local.mat);
lmout = summary(lmdl)$coefficients[2,c("Estimate","t value","Pr(>|t|)")];
print( lmout );


# Results from Matrix eQTL and "lm" must agree
stopifnot(all.equal(lmout, rez, check.attributes=FALSE));
