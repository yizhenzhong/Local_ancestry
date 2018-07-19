pkgname <- "LAMatrix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('LAMatrix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("LAMatrix_main")
### * LAMatrix_main

flush(stderr()); flush(stdout())

### Name: LAMatrix_main
### Title: Efficient eQTL Mapping with Local Ancestry
### Aliases: LAMatrix_main

### ** Examples

library("LAMatrix")

n = 200;# Number of columns (samples)
nc = 10;# Number ofs covariates

cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);# Generate the covariates
snps.mat = floor(runif(n, min = 0, max = 3)); # Generate the genotype
local.mat = floor(runif(n, min = 0, max = 3)); # Generate the local ancestry

# Generate the expression vector
gene.mat = cvrt.mat %*% rnorm(nc) + rnorm(n) + 0.5 * snps.mat + 1 + 0.3 *local.mat;


#Create SlicedData objects
snps = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
genes = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrts = SlicedData$new( t(cvrt.mat) );
locals = SlicedData$new( matrix(local.mat,nrow=1))

# Produce no output files
filename = NULL; # tempfile()

modelLOCAL = 930507L;
me = LAMatrix_main(
 snps = snps,
 gene = genes,
 cvrt = cvrts,
 local = locals,
 output_file_name = filename,
 pvOutputThreshold = 1,
 useModel = modelLOCAL,
 verbose = TRUE,
 pvalue.hist = FALSE);


# Pull results - t-statistic and p-value
beta = me$all$eqtls$beta;
tstat = me$all$eqtls$statistic;
pvalue = me$all$eqtls$pvalue;
rez = c(beta = beta, tstat = tstat, pvalue = pvalue);
print(rez)

# Results from linear
lmdl = lm( gene.mat ~ snps.mat + cvrt.mat + local.mat);
lmout = summary(lmdl)$coefficients[2,c("Estimate","t value","Pr(>|t|)")];
print( lmout );



cleanEx()
nameEx("modelLOCAL")
### * modelLOCAL

flush(stderr()); flush(stdout())

### Name: modelLOCAL
### Title: Constant for 'LAMatrix_main'
### Aliases: modelLOCAL
### Keywords: datasets

### ** Examples

library("LAMatrix")

n = 200;# Number of columns (samples)
nc = 10;# Number ofs covariates

cvrt.mat = 2 + matrix(rnorm(n*nc), ncol = nc);# Generate the covariates
snps.mat = floor(runif(n, min = 0, max = 3)); # Generate the genotype
local.mat = floor(runif(n, min = 0, max = 3)); # Generate the local ancestry

# Generate the expression vector
gene.mat = cvrt.mat %*% rnorm(nc) + rnorm(n) + 0.5 * snps.mat + 1 + 0.3 *local.mat;


#Create SlicedData objects
snps = SlicedData$new( matrix( snps.mat, nrow = 1 ) );
genes = SlicedData$new( matrix( gene.mat, nrow = 1 ) );
cvrts = SlicedData$new( t(cvrt.mat) );
locals = SlicedData$new( matrix(local.mat,nrow=1))

# Produce no output files
filename = NULL; # tempfile()

modelLOCAL = 930507L;
me = LAMatrix_main(
 snps = snps,
 gene = genes,
 cvrt = cvrts,
 local = locals,
 output_file_name = filename,
 pvOutputThreshold = 1,
 useModel = modelLOCAL,
 verbose = TRUE,
 pvalue.hist = FALSE);


# Pull results - t-statistic and p-value
beta = me$all$eqtls$beta;
tstat = me$all$eqtls$statistic;
pvalue = me$all$eqtls$pvalue;
rez = c(beta = beta, tstat = tstat, pvalue = pvalue);
print(rez)

# Results from linear
lmdl = lm( gene.mat ~ snps.mat + cvrt.mat + local.mat);
lmout = summary(lmdl)$coefficients[2,c("Estimate","t value","Pr(>|t|)")];
print( lmout );




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
