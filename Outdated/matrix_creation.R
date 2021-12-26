tab_trophic=read.delim("S:/2018/pglemasle/Data/Script/tab_trophic.txt")
colnames(tab_trophic)[1]="Predator"
tab_releve=read.delim("S:/2018/pglemasle/Data/Script/tab_releve.txt")
envi_data=read.delim("S:/2018/pglemasle/Data/Script/envi_data.txt")




#############♦

tab_releve=tab_releve[tab_releve[,1]=="Mare",]


tab_releve=droplevels(tab_releve)
tab_releve[,2]=as.character(tab_releve[,2], trim = TRUE)
tab_releve[,2]=as.factor(tab_releve[,2])


#####

envi_data[,2]=gsub("M_", "", envi_data[,2], fixed = TRUE)
envi_data[,2]=sub("^0+", "", envi_data[,2])
envi_data[,2]=factor(envi_data[,2])

####

test=as.factor(tab_releve[!duplicated(tab_releve[,2]),2])
test=as.factor(tab_releve[,2])

####

test[(test %in% envi_data[,2])]

tab_releve=tab_releve[(test %in% envi_data[,2]),]
tab_releve=droplevels(tab_releve)

####

doublonPhyto=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonPhyto[,c(5,6,8,9)]=NA
doublonPhyto[,7]="Phytoplankton"
doublonZoo=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonZoo[,c(5,6,8,9)]=NA
doublonZoo[,7]="Zooplankton"
doublonDetri=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonDetri[,c(5,6,8,9)]=NA
doublonDetri[,7]="Detritus"


tab_releve=rbind(tab_releve,doublonPhyto,doublonZoo,doublonDetri)

rm(doublonZoo,doublonPhyto,doublonDetri)

####CrÃ©ation de la matrice rÃ©gionnal

create_matrix=function(x){
  
  nm=paste("Matrix_",deparse(substitute(tab_trophic)),sep="")
  
  temp=matrix(0,length(levels(x$Predator)),length(levels(x$Predator)))
  
  colnames(temp)=levels(x$Predator)
  row.names(temp)=levels(x$Predator)
  
  pb <- txtProgressBar(min = 0, max = length(x[,1]), style = 3)
  
  for (i in 1:length(x[,1])) {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
    for (j in 1:length(levels(x$Predator))) {
      for (k in 1:length(levels(x$Predator))) {
        if(x[i,1]==colnames(temp)[j] & x[i,2]==row.names(temp)[k]){
          temp[k,j]=1
        }
        
      }
      
    }
  }
  close(pb)
  assign(nm,temp, env = .GlobalEnv)
}

create_matrix(tab_trophic)

###Passage des relevÃ©es aux matrices locales


create_matrix_marre=function(x,y){
  start.time <- Sys.time()
  
  marres=levels(x$Id_Sites)
  
  pb <- txtProgressBar(min = 0, max = length(marres), style = 3)
  
  for (o in 1:length(marres)) {
    Sys.sleep(0.1)
    setTxtProgressBar(pb, o)
    
    temp_site=tab_releve[tab_releve[,2]==marres[o],]
    temp_site[,2]=droplevels(temp_site[,2])
    temp_site[,7]=droplevels(temp_site[,7])
    
    nm=paste("Matrix_site_",marres[o],sep="")
    
    
    temp=matrix(0,length(levels(temp_site$Taxon)),length(levels(temp_site$Taxon)))
    
    colnames(temp)=levels(temp_site$Taxon)
    row.names(temp)=levels(temp_site$Taxon)
    
    
    for (i in 1:length(y[,1])) {
      for (j in 1:length(levels(temp_site$Taxon))) {
        for (k in 1:length(levels(temp_site$Taxon))) {
          if(y[i,1]==colnames(temp)[j] & y[i,2]==row.names(temp)[k] ){
            temp[k,j]=1
          }
          
        }
        
      }
    }
    
    assign(nm,temp, env = .GlobalEnv)
  }
  close(pb)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
}

create_matrix_marre(tab_releve,tab_trophic)




####


library(enaR)

descriptors=as.data.frame(matrix(nrow=length(levels(tab_releve$Id_Sites)),ncol=15))
colnames(descriptors)=c("Site","n","L","C","LD","ppr","lam1A","mlam1A","rho","R",
                        "d","no.scc","no.scc.big","pscc", "Richesse_spe")
name_site=levels(tab_releve$Id_Sites)
name_matrix=1:length(levels(tab_releve$Id_Sites))

for (i in 1:length(levels(tab_releve$Id_Sites))) {
  descriptors[i,1]=name_site[i]
  name_matrix[i]=paste("Matrix_site_",name_site[i],sep="")
  descriptors[i,2:14]=enaStructure(get(name_matrix[i]))$ns
  descriptors[i,15]=length(get(name_matrix[i])[,1])

}

####

library(igraph)

g1=graph_from_adjacency_matrix(Matrix_tab_trophic)
plot(g1, rescale = FALSE, ylim=c(0,4),xlim=c(-4,4), asp = 0)
tkplot(g1, canvas.width = 500, canvas.height = 500,vertex.size=4,
       vertex.label.dist=1, vertex.color="red", edge.arrow.size=0.5,label.cex=24,label.color="black")
modularity(g1,membership(cluster_walktrap(g1)))


