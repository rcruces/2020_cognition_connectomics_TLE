# ------------------------------------------------------------- #
#  Regularized Canonical Correlation Analysis 

# Repository path
home=paste0(Sys.getenv("KOTI"),'git_here/cognition_connectomics_TLE')
setwd(home)

# Loads the matrices
load('databases/rcca-XYmatrices.RData')

# Loads the functions and libraries
source('code/functions_cca.R')
library(CCA)

# ------------------------------------------------------------- #
#  Loads all the data 
# sm subject metrics, only patients, excluding IQ
Y <- sm[sm$clust!=0,9:17]
# nm network metrics, only patients
X <- nm[sm$clust!=0,]

# ------------------------------------------------------------- #
# Estimates the regularizatoin parameters (lambdas)
# WARNING the exahustive search might take a while
# regul <- estim.regul.v3(X, Y, plt = FALSE, grid1 = seq(0.1, 0.95, l=20), grid2 = seq(0.01, 0.25, l=20))
# I recommend to load the already optimized regularization parameters: l1=0.4122449, l2=0.2130204
regul <- readRDS(paste0(Sys.getenv("KOTI"),'git_here/epilepsia_processing/code/CCA/crossValidation/cca-cv-criterion_gridfine_ultra.rds'))
# plots the Cross Validation grid
img.estim.regul(regul)
L1=regul$lambda1
L2=regul$lambda2

# Estimates de confidence intervals and significance of the rCCA
Nrep=10000
# res.stats <- rcca.stats(X, Y, l1=0.4122449, l2=0.2130204, Nrep)
load('databases/rcca_stats.RData')

# Run the actual model of Regularized CCA
reg.cca <- rcc(X, Y, lambda1 = 0.4122449, lambda2 = 0.2130204)

# ------------------------------------------------------------- #
# Plots the canonical correlation
plot.cor(IC = res.stats, Title = 'Model')

# Plots V1 vs U1
plot.scores(reg.cca, sm$PR, Title = '', Clust = sm$clust!=0, sm = sm)

# Plots Y-loadings
sort.Yloads(rCCA = reg.cca, IC = res.stats, Xlim = c(-1, 0), Title = '')

# Plots Canonical variate 1 vs Canonical Variate 2 for X and Y loadings
Col <- rep(c('purple3', 'chartreuse3', 'orange'), each=dim(X)[2]/3)
plt.var.v3(res = reg.cca, d1 = 1, d2 = 2, Ycol = '#00CCCC', Xcol = Col, var.label = TRUE)

# Plots the subcortical loadings

# Plots the cortical loadings (as a Network)

