estim.regul.v3(X,Y, )
estim.regul.v2(X, Y, plt = FALSE, grid1 = seq(0.01, 0.95, l=5), grid2 = seq(0.01, 0.27, l=5))


save(nm,sm,file = "/host/yeatman/local_raid/rcruces/git_here/cognition_connectomics_TLE/databases/cca-XYmatrices.RData")
# Repository path
home=paste0(Sys.getenv("KOTI"),'git_here/cognition_connectomics_TLE')
setwd(home)

# Loads the matrices
load('databases/cca-XYmatrices.RData')

# Loads the functions
source('code/functions_cca.R')

# ------------------------------------------------------------- #
#  Regularized Canonical Correlation Analysis 
# sm subject metrics, only patients, excluding IQ
Y <- sm[sm$clust!=0,9:17]
# nm network metrics, only patients
X <- nm[sm$clust!=0,]

# Estimates the regularizatoin parameters (lambdas)

res.rcc <- rcc(X, Y, lambda1 = 0.4122449, lambda2 = 0.2130204)
