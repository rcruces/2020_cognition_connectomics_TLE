#-----------------------------------------#
####   LOADS matrix ####  
#-----------------------------------------#
mtx.load<- function(fichier) {
  # Nombre e Id del archivo
  id<-as.character(gsub("\\D", "", fichier))
  
  # Carga la matrix
  mtx<-as.matrix(read.csv(fichier,sep = ",",header = FALSE))
  
  # Agrega las etiqutas a las columnas
  rows<-colnames(mtx)<-rownames(mtx)<-1:dim(mtx)[1]
  
  # Remueve las filas vacias (118,76,42)
  #zero<-which(apply(mtx,2,mean)==0)
  zero<-c(118,76,42)
  mtx<-mtx[-zero,-zero]
  rows<-rows[-zero]
  
  # Counts the total number of stream lines
  diag(mtx)<-0
  
  # Mirror matrix Completa el triangulo inferior
  mtx[lower.tri(mtx)] <- t(mtx)[lower.tri(mtx)]
  
  tmp<-list("id"=id,"matrix"=mtx,"rows"=rows)
  return(tmp)
}