
#IMPORTACIÓN DE DATOS
workingDir <- "C:/Users/usuario/Desktop/Curso14-15/Aprendizaje_Computacional/Practica1"
setwd(workingDir)
datos <- read.delim2("spospread.txt", header=T)
datos.filt <- read.delim2("sporulation-filtered.txt", header=T)
datos.filt.num <- as.matrix(datos.filt[1:474,2:8])
dimnames(datos.filt.num)[[1]] <- datos.filt[1:474,1]
str(datos.filt.num)

#cLUSTERING POR RED AUTO-ORGANIZADA (SOM):
library(kohonen)
datos.filt.sc <- t(scale(t(datos.filt.num)))
datos.som <- som(data=datos.filt.sc, grid=somgrid(xdim=7, ydim=1))
plot(datos.som, main="SOM Plot")
#cuántos genes se han clasificado en cada cluster (lo guardamos en un vector de int):
datos.som$unit.classif
som.n.clust <- as.integer(summary(as.factor(datos.som$unit.classif)))
som.n.clust

#CLUSTERING JERÁRQUICO (HC):
misdatos.hc<-t(scale(t(datos.filt.num)))
hr<-hclust(as.dist(1-cor(t(misdatos.hc),method="pearson")),method="average")
plot(as.dendrogram(hr), main="Cluster dendrogram")


mycl7<-cutree(hr,k=7)
#cuantos genes se han clasificado en cada cluster (lo guardamos en un vector de int):
hc.n7.clust<-as.integer(summary(as.factor(mycl7)))
hc.n7.clust
#analizando el dendrograma y los datos de genes almacenados en cada cluster, observamos que la 
#agrupación es muy desequilibrada (hay 4 grupos con muy pocos genes), por lo que variamos K:
mycl4<-cutree(hr,k=4)
#cuantos genes se han clasificado en cada cluster (lo guardamos en un vector de int):
hc.n4.clust<-as.integer(summary(as.factor(mycl4)))
hc.n4.clust

#comprobamos que al disminuir el número de clusters, la agrupación sigue siendo muy desequilibrada,
#por lo que probamos con otro método diferente de aglomeración.
hc1<-hclust(as.dist(1-cor(t(misdatos.hc),method="pearson")),method="complete")
#para k=7:
mycl1.7<-cutree(hc1,k=7)
#cuantos genes se han clasificado en cada cluster (lo guardamos en un vector de int):
hc1.n7.clust<-as.integer(summary(as.factor(mycl1.7)))
hc1.n7.clust
#para k=6:
mycl1.6<-cutree(hc1,k=6)
#cuantos genes se han clasificado en cada cluster (lo guardamos en un vector de int):
hc1.n6.clust<-as.integer(summary(as.factor(mycl1.6)))
hc1.n6.clust
#para k=5:
mycl1.5<-cutree(hc1,k=5)
#cuantos genes se han clasificado en cada cluster (lo guardamos en un vector de int):
hc1.n5.clust<-as.integer(summary(as.factor(mycl1.5)))
hc1.n5.clust

#pintamos el dendrograma del clustering con el método "complete", marcando con una línea roja el
#"corte" del árbol para K=7 y con una línea azul, para K=6.
plot(as.dendrogram(hc1))
abline(h=1.0, col="red")
abline(h=1.2, col="blue")

#creamos una función para los colores del heatmap
my.colorFct <- function(n=50,low.col=0.45,high.col=1,saturation=1) {
  if(n<2)stop("n must be greater than 2")
  n1<-n%/%2
  n2<-n-n1
  c(hsv(low.col,saturation,seq(1,0,length=n1)),hsv(high.col,saturation,seq(0,1,length=n2)))
}
#pintamos el resultado del clustering jerárquico con método "complete" y K=6 en un heatmap
mycolhc<-sample(rainbow(256))
mycolhc<-mycolhc[as.vector(mycl1.6)]

hc<-hclust(as.dist(1-cor(misdatos.hc,method="spearman")),method="complete")

heatmap(datos.filt.num,Rowv=as.dendrogram(hc1),Colv=as.dendrogram(hc),col=my.colorFct(low.col=0.6,
        high.col=1,saturation=1),scale="row",RowSideColors=mycolhc,verbose=TRUE)


#VALIDACIÓN:
#creamos una función que realice SOM eliminando una columna de tiempo
som.func <- function(datos, t) {
  datos.val <- datos[,-t]
  datos.val.sc <- t(scale(t(datos.val)))  
  datos.val.som <- som(data=datos.val.sc, grid=somgrid(xdim=7, ydim=1))
  return(datos.val.som$unit.classif)
}
#llamamos a la función eliminando una columna en cada iteración y guardamos los resultado en una matriz
#(la primera columna contiene el clustering inicial sin eliminar instantes de t, para comparar)
som.val.mat<-matrix(c(datos.som$unit.classif),nrow=474,ncol=1)
for(i in 1:7) {
  som.val.mat<-cbind(som.val.mat,som.func(datos.filt.num, i))
}
som.val.mat

#creamos una función que realice HC eliminando una columna de tiempo
hc.func <- function(datos, t){
  datos.val <- datos[,-t]
  datos.val.sc<-t(scale(t(datos.val)))  
  hr.val<-hclust(as.dist(1-cor(t(datos.val.sc),method="pearson")),method="complete")
  mycl.val<-cutree(hr.val,k=6)
  return(mycl.val)
}
#llamamos a la función eliminando una columna en cada iteración y guardamos los resultado en una matriz
#(la primera columna contiene el clustering inicial sin eliminar instantes de t, para comparar)
hc.val.mat<-matrix(c(mycl1.6), nrow=474, ncol=1)
for(i in 1:7) {
  hc.val.mat<-cbind(hc.val.mat,hc.func(datos.filt.num, i))
}
hc.val.mat

#creamos una función para la validación APN
apn.func <- function (val.mat) {
  apn<-0  
  sum<-0
  for(j in 2:dim(val.mat)[2]) {
    for(i in 1:dim(val.mat)[1]) {
      if(val.mat[i,j]!=val.mat[i,1]){
        sum<-sum+1
      }
    }
  }
  apn <- sum/length(val.mat[,2:8])
  return(apn)
}

#llamamos a la función APN con el método hc:
hc.apn <- apn.func(hc.val.mat)
hc.apn

#llamamos a la función APN con el método som:
som.apn <- apn.func(som.val.mat)
som.apn

#creamos una función para la validación AD
ad.func <- function (val.mat) {
  ad<-0  
  sum<-0
  for(j in 2:dim(val.mat)[2]) {
    for(i in 1:dim(val.mat)[1]) {
      if(val.mat[i,j]!=val.mat[i,1]){
        sum<-sum+1
      }
    }
  }
  ad <- sum/length(val.mat[,2:8])
  return(ad)
}


#VALIDACIÓN UTILIZANDO EL PAQUETE clValid:
library(clValid)
#Estabilidad (APN, AD, ADM):
stab.valid <- clValid(datos.filt.num, 4:10, clMethods=c("hierarchical","som"), validation="stability", 
                    metric="correlation", method="complete")
summary(stab.valid)
plot(stab.valid, measure=c("APN","AD","ADM"),legend=FALSE)
plot(nClusters(stab.valid),measures(stab.valid,"APN")[,,1],type="n",axes=F, xlab="",ylab="")
legend("center", clusterMethods(stab.valid), col=1:9, lty=1:9, pch=paste(1:9))

#Interna (Silhuette):
intern.valid <- clValid(datos.filt.num, 4:10, clMethods=c("hierarchical","som"), validation="internal", 
                      metric="correlation", method="complete")
summary(intern.valid)
plot(intern.valid, measures=measNames(intern.valid)[3])
