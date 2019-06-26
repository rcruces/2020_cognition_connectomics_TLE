#---------------------------------------------------------------#
####  Reads each Matrix and correct for ICV and Hippocampal Volume ####  
#---------------------------------------------------------------#
# Libraries

# Change this variable to your home path
koti="~/git_here/cognition_connectomics_TLE/"

# Sets the working directory
setwd(koti)

# Load libraries and functions
source("code/functions.R")

# Reads the files
docs<-list.files("./matrices/")

# Creates an emty matrix
N<-length(docs)
mtxs<-array(0,dim = c(162,162,N))
lab<-c()

# Reads and uploads all matrices
# NOTE: use COMMA separated file instead of space
for (i in 1:N){
  tmp<-mtx.load(paste0("matrices/",docs[i]))
  mtxs[,,i]<-tmp$matrix
  lab<-c(lab,tmp$id)
  print(paste0("Reading files from ",tmp$id))
  # QC png
  # png(filename = paste0("QC_images/",as.character(gsub("\\D", "", docs[i])),"_matrix.png"),width = 800,height = 800)
  # corrplot(mtx_log(mtxs[,,i]),is.corr = F,tl.col="black",method="color",cl.pos = "r", title = lab[i])
  # dev.off()
}
rows<-tmp$rows
rm(tmp, i)
