require(cheddar)
require(igraph)
require(NetIndices)

setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
#setwd("S:/2018/pglemasle/Data/Script/")

#source("S:/2018/pglemasle/Data/Script/matrix_function.R")
source("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/matrix_function.R")
tab_trophic=read.delim("tab_trophic.txt")

#create_matrix(tab_trophic)

# Import manual matrix
latin=read.delim("latin.txt", header=F)
Matrix_Reg_Manual=read.delim2("Matrix_Manual_Chekck.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
row.names(Matrix_Reg_Manual)=latin$V3
colnames(Matrix_Reg_Manual)=latin$V3

g=graph_from_adjacency_matrix(Matrix_Reg_Manual, mode="undirected")
g1=graph_from_adjacency_matrix(Matrix_Reg_Manual[c('Detritus','Macrophyte','Phytoplankton'),
                                                 c('Detritus','Macrophyte','Phytoplankton')], mode="undirected")
troph=TrophInd(Matrix_Reg_Manual)
row.names(troph)=latin$V3



name_2=c("A.aquaticus","Bivalvia","Culicidae","Ephemeroptera",
         "Lymnaea","Physa","Potamogyrus")

for (i in 1:length(name_2)) {
  troph[row.names(troph)==name_2[i],1]=2
}

troph[row.names(troph)=='Detritus',1]=1
troph[row.names(troph)=='Macrophyte',1]=1
troph[row.names(troph)=='Phytoplankton',1]=1



quant=quantile(troph$TL, probs = seq(0, 1, 0.14))

layout.matrix.1<-matrix(
  nrow=length(V(g)),  # Rows equal to the number of vertices
  ncol=2)

layout.matrix.1[troph$TL==quant[1],1]=c(1.5,4.5,7.5)
layout.matrix.1[troph$TL>quant[1] & troph$TL<quant[2],1]=4.5
layout.matrix.1[troph$TL==quant[2],1]=seq(1,length(troph$TL[troph$TL==quant[2]]))
layout.matrix.1[troph$TL>quant[2] & troph$TL<quant[3],1]=seq(1,length(troph$TL[troph$TL>quant[2] & troph$TL<quant[3]]))
layout.matrix.1[troph$TL>=quant[3] & troph$TL<quant[4],1]=seq(1,length(troph$TL[troph$TL>=quant[3] & troph$TL<quant[4]]))
layout.matrix.1[troph$TL>=quant[4] & troph$TL<quant[5],1]=seq(1,length(troph$TL[troph$TL>=quant[4] & troph$TL<quant[5]]))
layout.matrix.1[troph$TL>=quant[5] & troph$TL<quant[6],1]=seq(1,length(troph$TL[troph$TL>=quant[5] & troph$TL<quant[6]]))
layout.matrix.1[troph$TL>=quant[6] & troph$TL<quant[7],1]=seq(1,length(troph$TL[troph$TL>=quant[6] & troph$TL<quant[7]]))
layout.matrix.1[troph$TL>=quant[7] & troph$TL<quant[8],1]=seq(1,length(troph$TL[troph$TL>=quant[7] & troph$TL<quant[8]]))
layout.matrix.1[troph$TL>=quant[8],1]=seq(1,length(troph$TL[troph$TL>=quant[8]]))

layout.matrix.1[levels(tab_trophic$Predator)=="Gasteropode spp.",1]=5.5
layout.matrix.1[levels(tab_trophic$Predator)=="Asellus aquaticus",1]=1.5
layout.matrix.1[levels(tab_trophic$Predator)=="Glossiphonie",1]=5.5
layout.matrix.1[levels(tab_trophic$Predator)=="Acilius spp.",1]=1.5






layout.matrix.1[,2]<-troph$TL # y-axis value based on trophic level

vertexLdist=rep(0.8,length(troph$TL))
names(vertexLdist)=latin$V3
#vertexLdist[levels(tab_trophic$Predator)=="Gasteropode spp."]=-0.5
#vertexLdist[levels(tab_trophic$Predator)=="Asellus aquaticus"]=-0.5
#vertexLdist[levels(tab_trophic$Predator)=="Breme spp."]=-0.1
vertexLdist["N.glauca"]=0.3
vertexLdist["Gammarus"]=0.3


par(mar=c(1,1,2,4),mfrow=c(1,1))

plot.igraph(g,
            vertex.label.cex=4,
            vertex.size=5,
            color="grey95",
            edge.arrow.size=.25,
            edge.width=1.7,
            layout=layout.matrix.1,
            canvas.width = 1000, canvas.height = 1000,
            vertex.label.dist=vertexLdist,
            edge.arrow.size=0.8,
            vertex.color = 'grey50',
            edge.color="grey90",
            vertex.frame.color = "grey10",
            vertex.label.color="black")


