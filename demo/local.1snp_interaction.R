library(LAMatrix)


# Number of columns (samples)
n = 200;

# Number ofs covariates
nc = 10;

# Generate the covariates

cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);
snps.mat = floor(runif(n, min = 0, max = 3)); 
local.mat = floor(runif(n, min = 0, max = 3));

# Generate the interaction
int.mat = snps.mat * local.mat

# Generate the vectors with genotype and expression variables

gene.mat = cvrt.mat %*% rnorm(nc) + rnorm(n) + 0.5 * snps.mat + 0.3 *local.mat + 0.6 * int.mat;

# Create 3 SlicedData objects for the analysis
snps = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
genes = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrts = SlicedData$new( t(cvrt.mat) );
locals = SlicedData$new( matrix(local.mat,nrow=1))
ints = SlicedData$new( matrix(int.mat,nrow=1))

# Produce no output files
filename = NULL; # tempfile()

me = LAMatrix_main(
  snps = ints,
  gene = genes,
  cvrt = cvrts,
  local = snps,
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
lmdl = lm( gene.mat ~ snps.mat + cvrt.mat + int.mat);
lmout = summary(lmdl)$coefficients;
print( lmout[nrow(lmout),c("Estimate","t value","Pr(>|t|)")] )


# Results from Matrix eQTL and "lm" must agree
stopifnot(all.equal(lmout, rez, check.attributes=FALSE));

