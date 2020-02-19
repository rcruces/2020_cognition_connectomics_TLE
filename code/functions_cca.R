# Functions additional to the CCA package
# Please see:
# https://cran.r-project.org/web/packages/CCA/CCA.pdf
# https://core.ac.uk/download/pdf/6303071.pdf

# For the original rcc and CCA package cite:
#   González, I., Déjean, S., Martin, P. G., & Baccini, A. (2008). 
#   CCA: An R package to extend canonical correlation analysis. Journal of Statistical Software, 23(12), 1-14.

# For the usage of the implementations found here cite:
#   Multidimensional Associations between Cognition and Connectome Organization in Temporal Lobe Epilepsy
#   Raúl Rodríguez-Cruces, Boris C. Bernhardt, Luis Concha; bioRxiv 675884; doi: https://doi.org/10.1101/675884

# -------------------------------------------------------------------------------------------------------------------------- #
#### FUNCTIONS for CCA estimations ####
# -------------------------------------------------------------------------------------------------------------------------- #
estim.regul.v2 <- function (X, Y, grid1 = NULL, grid2 = NULL, plt = TRUE) {
  # loo function corrected for dot product calculation
  # absolute value of the cross validation score
  require("CCA")
  loo <- function (X, Y, lambda1, lambda2) {
    n = nrow(X)
    xscore = vector(mode = "numeric", length = n)
    yscore = vector(mode = "numeric", length = n)
    for (i in 1:n) {
      Xcv = X[-i, ]
      Ycv = Y[-i, ]
      res = rcc(Xcv, Ycv, lambda1, lambda2)
      xscore[i] = as.vector(t(X[i, ])) %*% res$xcoef[, 1]
      yscore[i] = as.vector(t(Y[i, ])) %*% res$ycoef[, 1]
    }
    cv = abs(cor(xscore, yscore, use = "pairwise"))
    print(paste("lambda1:",lambda1,"and lambda2:",lambda2," has a CV: ",round(cv,3)))
    return(invisible(cv))
  }
  
  if (is.null(grid1)) {
    grid1 = seq(0.001, 1, length = 5)
  }
  if (is.null(grid2)) {
    grid2 = seq(0.001, 1, length = 5)
  }
  grid = expand.grid(grid1, grid2)
  res = apply(grid, 1, function(x) {
    loo(X, Y, x[1], x[2])
  })
  res.grid = cbind(grid, res)
  mat = matrix(res, nc = length(grid2))
  if (isTRUE(plt)) {
    img.estim.regul(list(grid1 = grid1, grid2 = grid2, mat = mat))
  }
  opt = res.grid[res.grid[, 3] == max(res.grid[, 3]), ]
  cat("  lambda1 = ", opt[[1]], "\n", " lambda2 = ", opt[[2]],
      "\n", "CV-score = ", opt[[3]], "\n")
  return(invisible(list(lambda1 = opt[[1]], lambda2 = opt[[2]],
                        CVscore = opt[[3]], grid1 = grid1, grid2 = grid2, mat = mat)))
}

# -------------------------------------------------------------------------------------------------------------------------- #
estim.regul.v3 <- function (X, Y, grid1 = NULL, grid2 = NULL, plt = TRUE) {
  # Uses multicore for the loo procedure
  # Correct the dot product multiplication
  # Absolute value to find the maxmum cross validation score
  # Adds a progress bar to the apply function
  require("foreach")
  require("doMC")
  require("CCA")
  require("pbapply")
  # loo function corrected for dot product
  loo <- function (X, Y, lambda1, lambda2) {
    n = nrow(X)
    registerDoMC()
    CVs <- foreach(i=1:n,.verbose = FALSE) %dopar% {
      res<-rcc(X[-i, ], Y[-i, ], lambda1, lambda2)
      c(as.vector(t(X[i, ])) %*% res$xcoef[, 1], as.vector(t(Y[i, ])) %*% res$ycoef[, 1])}
    CVs <-matrix(unlist(CVs),n,2,byrow = TRUE)
    cv = abs(cor(as.vector(CVs[,1]), as.vector(CVs[,2]), use = "pairwise"))
    return(invisible(cv))
  }
  
  if (is.null(grid1)) {
    grid1 = seq(0.001, 1, length = 5)
  }
  if (is.null(grid2)) {
    grid2 = seq(0.001, 1, length = 5)
  }
  grid = expand.grid(grid1, grid2)
  print(paste('[INFO] ... Estimation of the optimal regularization parameters out of',dim(grid)[1],'combinations'))
  res = pbapply(grid, 1, function(x) {
    loo(X, Y, x[1], x[2])
  })
  res.grid = cbind(grid, res)
  mat = matrix(res, nc = length(grid2))
  if (isTRUE(plt)) {
    img.estim.regul(list(grid1 = grid1, grid2 = grid2, mat = mat))
  }
  opt = res.grid[res.grid[, 3] == max(res.grid[, 3]), ]
  cat("  lambda1 = ", opt[[1]], "\n", " lambda2 = ", opt[[2]],
      "\n", "CV-score = ", opt[[3]], "\n")
  return(invisible(list(lambda1 = opt[[1]], lambda2 = opt[[2]],
                        CVscore = opt[[3]], grid1 = grid1, grid2 = grid2, mat = mat)))
}

# -------------------------------------------------------------------------------------------------------------------------- #
# rCCA.stats - Statistics for regularized Canonical Correlations
rcca.stats <- function(X, Y, l1, l2, Nrep=1000){
  start.time <- Sys.time()
  require(CCA)
  require(candisc)
  require(boot)
  # -------------------------------------------------------------------------- #
  # Estimates the Confidence Interval with bootstrap
  # function for the rCCA
  boot.cca <- function(Dataset, indx) {
    n <- length(indx)
    # Run rCCA with bootstraped index
    rcca <- rcc(X[indx,], Y[indx,], lambda1 = l1, lambda2 = l2)
    # Cognitive loadings dimension 1 & Canonical correlations
    c(rcca$cor, abs(rcca$scores$corr.Y.xscores[,1]), abs(rcca$scores$corr.X.xscores[,1]))
  }
  if (Nrep<1000) {print('WARNING ... confidence interval and significance might be compromise with Nrep lower than 1000')}
  # Run bootstrap
  print(paste("[INFO] ... Running",Nrep,"bootstraps in rcca"))
  bcca <- boot(X, statistic = boot.cca, R=Nrep)
  # Gets the CI fot each variable
  ci <- matrix(NA,nrow = length(bcca$t0),ncol = 3)
  ci[,1] <- bcca$t0
  for (i in 1:length(bcca$t0)) {
    CI <- try(boot.ci(boot.out = bcca,conf = 0.95, type = c("bca"), index=i)$bca[4:5])
    if (class(CI)=="try-error") { ci[i,2:3] <- c(NA, NA) }
    ci[i,2:3] <- as.numeric(CI)
  }
  # set ci as dataframe
  ci <- data.frame(ci)
  # Conditional test
  if (dim(ci)[2] < 3) {print("Something is wrong with the bootstrap you might need more Nboot, minimum 1000")}
  
  print("[DONE] ... Bootstrap for Confidence Intervals")
  
  # --------------------------------------------------------------------------
  # Estimates the null distribution with Permutation for significance (p-value)
  print(paste("[INFO] ... Running",Nrep,"permutations in rcca"))
  ## CCA permutation testing
  Pcca <- list(ccor=matrix(NA,nrow = dim(Y)[2],ncol = Nrep),
               Yperm=matrix(NA,nrow = dim(Y)[2],ncol = Nrep),
               Xperm=matrix(NA,nrow = dim(X)[2],ncol = Nrep))
  rownames(Pcca$Yperm) <- colnames(Y)
  rownames(Pcca$Xperm) <- colnames(X)
  
  ## CCA permutation testing with iteration
  n <- dim(X)[1]
  # create progress bar
  pb <- txtProgressBar(min = 0, max = Nrep, style = 3)
  for (j in 1:Nrep){
    # We labeled Rows as random to create a Null distribution
    Rowi <- sample(1:n,replace = FALSE)
    cca.perm <- rcc(X, Y[Rowi,], lambda1 = l1, lambda2 = l2)
    Pcca$ccor[,j] <- cca.perm$cor
    Pcca$Yperm[,j] <- cca.perm$scores$corr.Y.xscores[,1]
    Pcca$Xperm[,j] <- cca.perm$scores$corr.X.xscores[,1]
    #print(paste("Running permutation",j))
    # update progress bar
    setTxtProgressBar(pb, j)
  }
  
  # Computational Time
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
  # get pvalues
  res.rcc <- rcc(X,Y,l1,l2)
  grotR <- res.rcc$cor
  pval.cor <- rep(0,length(grotR))
  for (i in 1:length(grotR)){ pval.cor[i] = (sum(Pcca$ccor[i,] > grotR[i]))/Nrep }
  
  Yload <- res.rcc$scores$corr.Y.xscores[,1]
  pval.Y <- rep(0,length(Yload))
  for (i in 1:length(Yload)){ pval.Y[i] <- (sum(abs(Pcca$Yperm[i,]) > abs(Yload[i]))/Nrep) }
  
  Xload <- res.rcc$scores$corr.X.xscores[,1]
  pval.X <- rep(0,length(Xload))
  for (i in 1:length(Xload)){ pval.X[i] <- (sum(abs(Pcca$Xperm[i,]) > abs(Xload[i]))/Nrep) }
  
  # data frame of the statistics
  rcca.stats <- data.frame(ci)
  # Colnames of the statistics
  colnames(rcca.stats) <- c("orig","ci.bot","ci.up")
  # p-values
  rcca.stats$pval <- c(pval.cor,pval.Y,pval.X)
  # Rownames of the statistics
  rownames(rcca.stats) <- c(paste0("cancor",1:length(grotR)),colnames(Y),colnames(X))
  
  return(rcca.stats)
}

# -------------------------------------------------------------------------------------------------------------------------- #
#### FUNCTIONS for CCA data handling ####
# -------------------------------------------------------------------------------------------------------------------------- #
#   rCCA.data - Get XY alls loading from a CCA model  
cca.getXY <- function(Net){
  P <- function(Path){koti=Sys.getenv("KOTI"); return(paste0(koti,Path))}
  # Renames the colums of the Path Lenght
  colnames(Net$K) <- paste0("k.",Net$roi.id)
  colnames(Net$CC) <- paste0("cc.",Net$roi.id)
  colnames(Net$L) <- paste0("L.",Net$roi.id)
  Net$id <- as.numeric(Net$id)
  
  # Network measurements matrices
  nm <- cbind(Net$K,Net$CC,Net$L)
  rownames(nm) <- Net$id
  
  # Subject Measures (Cognitive and Clinical)
  sm.data <- merge(read.csv(P("git_here/epilepsia_processing/code/CCA/subject_measures.csv")),
                   read.csv(P("git_here/epilepsia_processing/code/CCA/aparc2009_fs.csv")), by="id")
  rownames(sm.data) <- sm.data$id
  sm <- sm.data[,2:19]
  Indx <- match(sm.data$id, Net$id)
  nm <- na.exclude(nm[Indx,])
  sm <- sm[rownames(nm),]
  
  # Matrices for rCCA model
  X <-nm[sm$clust!=0,]
  Y <- sm[sm$clust!=0,9:17]
  return(list(X,Y,sm))
}

# -------------------------------------------------------------------------------------------------------------------------- #
#   rCCA.data - Obtains ordered loadings in a list from CCA
gral.loadings <- function(rCCA,IC,p.thr,D, LUT){
  d <- dim(IC)[1]
  Xn <- length(rCCA$cor)*2+1
  rho <- rCCA$scores$corr.X.xscores[which(IC[Xn:d,4]<p.thr),D]
  pval <- IC[which(IC[Xn:d,4]<p.thr)+(Xn-1),4]
  net.nom <- names(rho)
  id <- as.numeric(gsub("L.","", gsub("cc.","",gsub("k.","",net.nom))))
  roi <- LUT[match(id, LUT$id),"roi"]
  loadings <- rbind(id, roi, rho, pval); colnames(loadings) <- net.nom
  W <- list()
  W$L <- loadings[,grep(pattern = "L.", colnames(loadings))]
  W$C <- loadings[,grep(pattern = "cc.", colnames(loadings))]
  W$K <- loadings[,grep(pattern = "k.", colnames(loadings))]
  return(W)
}

# -------------------------------------------------------------------------------------------------------------------------- #
# rCCA.data - Obtains ordered loadings in a list from CCA
get.subcort <-function(Net.loadings,ModelName,LUT, roi.thr){
  # roi.thr is the value for subcortical ROIs usually:
  #     Destrieux+Volbrain: 152 & Schaefer 1001
  subCort <- data.frame(t(Net.loadings[,Net.loadings[1,]>=roi.thr]))
  colnames(subCort) <- rownames(Net.loadings)
  Indx <- match(subCort$id, LUT$id)
  LUT.sliced <- LUT[Indx,]
  subCort$roi.nom <- as.character(sapply(as.character(LUT.sliced$label.name), function(x) strsplit(x, "-")[[1]][2]))
  subCort$roi.side <- as.character(sapply(as.character(LUT.sliced$label.name), function(x) strsplit(x, "-")[[1]][1]))
  subCort$rCCA.mod <- rep(ModelName,dim(subCort)[1])
  subCort$measure <- sapply(colnames(Net.loadings)[Net.loadings[1,]>=roi.thr], function(x) strsplit(x, ".1")[[1]][1])
  subCort <- subCort[,c("rCCA.mod","id","roi","roi.nom","roi.side","rho","pval","measure")]
  return(subCort)
}

# -------------------------------------------------------------------------------------------------------------------------- #
# cca.loadings - Get the loadings of any model in a table
cca.loadings <- function(rCCA,IC,p.thr,D){
  d <- dim(IC)[1]
  Xn <- length(rCCA$cor)*2+1
  rho <- rCCA$scores$corr.X.xscores[which(IC[Xn:d,4]<p.thr),D]
  pval <- IC[which(IC[Xn:d,4]<p.thr)+(Xn-1),4]
  net.nom <- names(rho)
  loadings <- rbind(rho, pval)
  return(loadings)
}

# -------------------------------------------------------------------------------------------------------------------------- #
#### PLOTS for CCA ####
# -------------------------------------------------------------------------------------------------------------------------- #
plt.var.v2 <- function (res, d1, d2, int = 0.5, var.label = FALSE, Xnames = NULL,
                        Ynames = NULL) {
  require("scales")
  scaleBetween <- function(X,a,b) {
    Fx <- ((b-a)*(X-min(X))/(max(X)-min(X)))+a
    return(Fx)}
  back.lines <- function(int){
    abline(v = 0, h = 0, lwd=1.5, lty=2,col="gray65")
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)),lty=2,col="gray65")
    lines(int * cos(seq(0, 2 * pi, l = 100)), int * sin(seq(0,2 * pi, l = 100)),lty=2,col="gray65")}
  if (!var.label) {
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1),
         xlab = paste("Dimension ", d1), ylab = paste("Dimension ",
                                                      d2))
    abline(v = 0, h = 0, lwd=1.5, lty=2)
    points(res$scores$corr.X.xscores[, d1], res$scores$corr.X.xscores[,
                                                                      d2], pch = 20, cex = 1.2, col = rep(c("midnightblue","forestgreen","orange"),each=162))
    points(res$scores$corr.Y.xscores[, d1], res$scores$corr.Y.xscores[,
                                                                      d2], pch = 24, cex = 0.7, col = "red3")
  }
  else {
    if (is.null(Xnames))
      Xnames = res$names$Xnames
    if (is.null(Ynames))
      Ynames = res$names$Ynames
    plot(0, bty='n',type = "n", xlim = c(-1, 1), ylim = c(-1, 1),
         xlab = paste("CCA dimension ", d1), ylab = paste("CCA dimension ",
                                                          d2))
    back.lines(int)
    Alpha <- scaleBetween(abs(res$scores$corr.X.xscores[, d1]),0.15,0.5) + scaleBetween(abs(res$scores$corr.X.xscores[, d2]),0.15,0.5)
    Cex <- scaleBetween(abs(res$scores$corr.X.xscores[, d1]),0.25,0.875) + scaleBetween(abs(res$scores$corr.X.xscores[, d2]),0.25,0.875)
    text(res$scores$corr.X.xscores[, d1], res$scores$corr.X.xscores[,d2], Xnames, col = alpha(rep(c("purple3","chartreuse3","orange"),each=162),Alpha), font = 1,cex=Cex)
    text(res$scores$corr.Y.xscores[, d1], res$scores$corr.Y.xscores[,d2], Ynames, col = "cyan3", font = 3,cex=1.5)
  }
}

# -------------------------------------------------------------------------------------------------------------------------- #
#  plt.var.v3: Improvement to the color for canonical scores plotting
plt.var.v3 <- function (res, d1, d2, int = 0.5, var.label = FALSE, Xnames = NULL,
                        Ynames = NULL, Xcol=NULL,Ycol=NULL, NO.alpha=NULL) {
  require("scales")
  scaleBetween <- function(X,a,b) {
    Fx <- ((b-a)*(X-min(X))/(max(X)-min(X)))+a
    return(Fx)}
  back.lines <- function(int){
    abline(v = 0, h = 0, lwd=1.5, lty=2,col="gray65")
    lines(cos(seq(0, 2 * pi, l = 100)), sin(seq(0, 2 * pi, l = 100)),lty=2,col="gray65")
    lines(int * cos(seq(0, 2 * pi, l = 100)), int * sin(seq(0,2 * pi, l = 100)),lty=2,col="gray65")}
  if (!var.label) {
    plot(0, type = "n", xlim = c(-1, 1), ylim = c(-1, 1),
         xlab = paste("Dimension ", d1), ylab = paste("Dimension ",
                                                      d2))
    abline(v = 0, h = 0, lwd=1.5, lty=2)
    points(res$scores$corr.X.xscores[, d1], res$scores$corr.X.xscores[,d2], pch = 20, cex = 1.2, col = rep("gray65",each=length(res$scores$corr.X.xscores[, d1])))
    points(res$scores$corr.Y.xscores[, d1], res$scores$corr.Y.xscores[,d2], pch = 24, cex = 0.7, col = "red3")
  }
  else {
    if (is.null(Xnames))
      Xnames = res$names$Xnames
    if (is.null(Ynames))
      Ynames = res$names$Ynames
    if (is.null(Xcol))
      Xcol = rep("gray25",each=length(res$scores$corr.X.xscores[, d1]))
    if (is.null(Ycol))
      Ycol = "red3"
    
    plot(0, bty='n',type = "n", xlim = c(-1, 1), ylim = c(-1, 1),
         xlab = paste("CCA dimension ", d1), ylab = paste("CCA dimension ",d2))
    back.lines(int)
    # Length of X matrix
    N <- length(Xnames)
    
    # Size and Alpha dependt of First & Second Canonical Variates
    Cex <- scaleBetween(abs(res$scores$corr.X.xscores[, d1]),0.25,0.875) + scaleBetween(abs(res$scores$corr.X.xscores[, d2]),0.25,0.875)
    Alpha =  scaleBetween(abs(res$scores$corr.X.xscores[1:N, d1]),0.05,0.5) + scaleBetween(abs(res$scores$corr.X.xscores[1:N, d2]),0.05,0.5)
    
    # Size and Alpha dependt of First Canonical Variate
    # Cex <- scaleBetween(abs(res$scores$corr.X.xscores[, d1]),0.5,1.75)
    # Alpha =  scaleBetween(abs(res$scores$corr.X.xscores[1:N, d1]),0.05,0.95)
    
    # No alpha
    Alpha[NO.alpha] <- 1
    
    # PLot labels with color, font and size
    text(res$scores$corr.X.xscores[, d1], res$scores$corr.X.xscores[,d2], Xnames, col = alpha(Xcol,Alpha), font = 1,cex=Cex)
    text(res$scores$corr.Y.xscores[, d1], res$scores$corr.Y.xscores[,d2], Ynames, col = Ycol, font = 3,cex=1.5)
  }
}

# -------------------------------------------------------------------------------------------------------------------------- #
# optim.color:  Function to optimize the colormap 
optim.color <- function(Data,Colors) {
  mtx <- as.matrix(Data)
  # Following code limits the lowest and highest color to 5%, and 95% of your range, respectively
  quantile.range <- quantile(mtx, probs = seq(0, 1, 0.01))
  palette.breaks <- seq(quantile.range["5%"], quantile.range["95%"], 0.1)
  # Find optimal divergent color palette (or set own)
  color.function <- colorRampPalette(Colors)
  color.palette  <- color.function(length(palette.breaks) - 1)
  return(list(color.palette=as.vector(color.palette),palette.breaks=as.vector(palette.breaks)))
}

# -------------------------------------------------------------------------------------------------------------------------- #
img.matcor2 <- function (correl, type = 1, Col, Xmain, Ymain) {
  matcorX = correl$Xcor
  matcorY = correl$Ycor
  matcorXY = correl$XYcor
  lX = ncol(matcorX)
  lY = ncol(matcorY)
  def.par <- par(no.readonly = TRUE)
  if (type == 1) {
    par(mfrow = c(1, 1), pty = "s")
    image(1:(lX + lY), 1:(lX + lY), t(matcorXY[nrow(matcorXY):1,
                                               ]), zlim = c(-1, 1), main = "XY correlation", col = Col,
          axes = FALSE, xlab = "", ylab = "")
    box()
    abline(h = lY + 0.5, v = lX + 0.5, lwd = 2, lty = 2)
    image.plot(legend.only = TRUE, zlim = c(-1, 1), col = Col,
               horizontal = TRUE)
  }
  else {
    layout(matrix(c(1, 2, 3, 3, 0, 0), ncol = 2, nrow = 3,
                  byrow = TRUE), widths = 1, heights = c(0.8, 1, 0.06))
    par(pty = "s", mar = c(2, 2, 2, 1))
    image(1:lX, 1:lX, t(matcorX[lX:1, ]), zlim = c(-1, 1),
          main = Xmain, axes = FALSE, xlab = "",
          ylab = "", col = Col)
    box()
    image(1:lY, 1:lY, t(matcorY[lY:1, ]), zlim = c(-1, 1),
          main = Ymain, col = Col, axes = FALSE,
          xlab = "", ylab = "")
    box()
    partXY = matcorXY[lX:1, (lX + 1):(lX + lY)]
    if (lX > lY) {
      partXY = matcorXY[(lX + lY):(lX + 1), 1:lX]
      lX = ncol(matcorY)
      lY = ncol(matcorX)
    }
    par(pty = "m", mar = c(5, 4, 2, 3), mai = c(0.8, 0.5,
                                                0.3, 0.4))
    image(1:lY, 1:lX, t(partXY), zlim = c(-1, 1), main = "Cross-correlation",
          axes = FALSE, xlab = "", ylab = "", col = Col)
    box()
    image.plot(legend.only = TRUE, zlim = c(-1, 1), horizontal = TRUE,
               col = Col, legend.width = 2.5)
  }
  par(def.par)
}

# -------------------------------------------------------------------------------------------------------------------------- #
#  rCCA.plots - Barplot of canonical correlations
plot.cor <- function (IC,Title='Canonical correlations'){
  require(extrafont)
  par(family="Verdana")
  Col <- ifelse(IC[1:9,4]<=0.05,"gray35","gray")
  Bar <- barplot(IC[1:9,1], xlab = "Dimension", ylab = "",main=Title, names.arg = 1:9, ylim = c(0,1),col = Col)
  for (i in 1:9){
    lines(rep(Bar[i],2), IC[i,2:3],lwd=2)
    lines(c(Bar[i]-0.2,Bar[i]+0.2), rep(IC[i,2],2),lwd=2)
    lines(c(Bar[i]-0.2,Bar[i]+0.2), rep(IC[i,3],2),lwd=2)
    if (IC[i,4]<0.05) {text(x=Bar[i],y=0.8,labels = "*",cex=3)}}}

# -------------------------------------------------------------------------------------------------------------------------- #
#  rCCA.plots - Canonical Scores U1 vs V1
plot.scores <- function (rCCA,Cog,Title,Clust,sm){
  require(extrafont)
  NetwkMeasures <- rCCA$scores$xscores[ , 1]
  CogniMeasures <- rCCA$scores$yscores[ , 1]
  sdr <- sort(CogniMeasures)
  sdr <- sdr[c(1, length(sdr))]
  ext <- match(sdr, CogniMeasures)
  # first and last
  Col <- c("gray30","darkblue","forestgreen","red3")[sm$clust+1]
  Col <- Col[Clust]
  Cex <- ((as.numeric(cut(Cog[Clust],breaks = 3))-1)*2)+1
  par(family="Verdana")
  plot(0,0,col=NA,bty='n',main=Title,
       ylab = "Network canonical scores",
       xlab = "Cognitive canonical scores",xlim=c(-3,2.5),ylim=c(-3,2.5))
  abline(h=0,v=0,lty=3,col="gray45",lwd=2)
  #
  points(CogniMeasures, NetwkMeasures, pch=21,cex=Cex,bg=alpha(Col,0.3),col=Col)
}

# -------------------------------------------------------------------------------------------------------------------------- #
# rCCA.plots -  Y loadigs Cognitive Scores Stick & Ball
sort.Yloads <- function(rCCA,IC,Xlim=c(-1,1),Title='Y loadigs'){
  require(extrafont)
  require(RColorBrewer)
  par(family="Verdana")
  scaleBetween <- function(X,a,b) {
    Fx <- ((b-a)*(X-min(X))/(max(X)-min(X)))+a
    return(Fx)}
  Val <- rCCA$scores$corr.Y.xscores[,1]
  
  # Same Order of Cognitive scores
  Orden <- rev(c("PS","VCI","WMI","PR","VMI","IMI", "DMI", "AMI", "VWMI"))
  Val <- Val[match(Orden, names(Val))]
  rbPal <- colorRampPalette(rev(brewer.pal(10,"RdYlBu")))
  npsy.names <- rev(c("Processing Speed", "Verbal Comprehension", "Working Memory", "Perceptual Reasoning","Visual Memory", "Immediate memory", "Delayed memory", "Auditory memory", "Visual working memory"))
  Cex <- scaleBetween(abs(Val),1,1.75)
  # Color match
  S <- gsub(" ","",format(seq(-1,1,0.01),nsmall=2))
  Q <- gsub(" ","",format(round(Val,2),nsmall=2))
  Col <- rbPal(length(S))
  Col <- Col[match(Q,S)]
  # plot
  y <- seq(0.05,0.95,length.out = 9)
  plot(NA,xlim=c(Xlim[1]-0.5,Xlim[2]),ylim=c(0,1),xlab="",ylab="",bty="n",las=2,axes = F,main=paste0("Cognitive Loadings-CanVariate 1\n",Title))
  axis(side = 1,at=seq(Xlim[1],Xlim[2],0.5),las=1)
  axis(side = 2,at=seq(0.05,0.95,length.out = 9),las=1,labels = npsy.names,col.ticks = NA,lwd = 0,pos = Xlim[1], line = NA)
  abline(v=seq(Xlim[1],Xlim[2],0.5),h=seq(0,1,length.out = 4),lty=2,col="gray65")
  
  # Confidence Intervals
  Sign <- sign(Val)
  IC.sort <- IC[match(Orden,rownames(IC)),1:3]
  points(IC.sort[,1]*Sign,y,cex=1.5,col=Col,bg=Col,pch=21)
  
  for (i in 1:9){lines(IC.sort[i,2:3]*Sign[i],
                       rep(y[i],2),lwd=3,col=Col[i])}
}

# -------------------------------------------------------------------------------------------------------------------------- #
# rCCA.plots - Plots the loadings of subcortical areas 
plot.subW <- function(Loadings, Title, LUT, Indx, roi.thr) {
  # require(BrainSlice) my package in DEVELOMPMENT
  # Match ROI-id with plot.subcortical position
  Sliced <- get.subcort(Loadings, Title, LUT, roi.thr)
  Sliced <- Sliced[match(Indx,Sliced$roi),]
  Sliced[is.na(Sliced)] <- 0
  Ind <- which(Sliced$id!=0)
  plot.subcortical(ID=Ind, Rho=Sliced$rho[Ind],Title = Title,FrameCol = NA,Labels = FALSE)
}


