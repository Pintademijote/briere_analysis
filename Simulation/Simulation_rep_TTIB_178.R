setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
#setwd("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/Data")

require("enaR")
require("igraph")
require("NetIndices")
require("nlme")
require("mgcv")
library("ggplot2")
library(scales)
source('utils_NAR.R') #This code is used to compute all the network metrics

# Import trophic link database
tab_trophic=read.delim("tab_trophic.txt")
colnames(tab_trophic)[1]="Predator"
# Import taxon occurrence data
tab_releve=read.delim("tab_releve.txt")
which(tab_releve$Id_Sites=="")
tab_releve=tab_releve[-which(tab_releve$Id_Sites==""),]
tab_releve=droplevels(tab_releve)

sites=levels(tab_releve$Id_Sites)
for (i in sites) {
  x=tab_releve[tab_releve$Id_Sites==i,]
  if(all(c("Aeshnidae","Anisoptere","Libellulidae") %in% x$Taxon)){
    tab_releve[tab_releve$Id_Sites==i & tab_releve$Taxon=="Anisoptere",7]="Libellulidae"
  }else{
    if(all(c("Aeshnidae","Anisoptere") %in% x$Taxon)){
      tab_releve[tab_releve$Id_Sites==i & tab_releve$Taxon=="Aeshnidae",7]="Anisoptere"
    }else{
      if(all(c("Libellulidae","Anisoptere") %in% x$Taxon)){
        tab_releve[tab_releve$Id_Sites==i & tab_releve$Taxon=="Libellulidae",7]="Anisoptere"
      }
    }
    
    
  }
}



tab_releve=tab_releve[tab_releve$Taxon!="Couleuvre a collier",]
tab_releve=droplevels(tab_releve)

# Import environmental variables
envi_data=read.table("envi_data.txt",header=T)
# Import distance from marsh
dist=read.delim("Dist_232.txt")
colnames(dist)[1]="Distance"
# Import pond surface
Surf=read.delim2("Surface.txt")
# Import manual matrix
Matrix_Reg_Manual=read.delim2("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/Matrix_Final.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
#data=read.delim("data.txt")

# Import function to create matrix
#source("S:/2018/pglemasle/Data/Script/matrix_function.R")

##########Marais/Mare

Surf[,1]=gsub("M_", "", Surf[,1], fixed = TRUE)
Surf[,1]=sub("^0+", "", Surf[,1])
Surf[,1]=factor(Surf[,1])


dist[,3]=gsub("M_", "", dist[,3], fixed = TRUE)
dist[,3]=sub("^0+", "", dist[,3])
dist[,3]=as.factor(dist[,3])
dist=droplevels(dist)

## Keep only ponds
#tab_releve=tab_releve[tab_releve[,1]=="Mare",]

## Reset factor levels
tab_releve=droplevels(tab_releve)
## Reset format
tab_releve[,2]=as.character(tab_releve[,2], trim = TRUE)
tab_releve[,2]=as.factor(tab_releve[,2])

test=as.factor(tab_releve[!duplicated(x[,2]),2])	# Keep only one line per site
test=as.factor(tab_releve[,2])


#### ADD non-observed taxa at all sites
doublonPhyto=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonPhyto[,c(5,6,8,9)]=NA				# Remove useless data but keep matrix size		
doublonPhyto[,7]="Phytoplankton"
doublonZoo=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonZoo[,c(5,6,8,9)]=NA
doublonZoo[,7]="Zooplankton"
doublonDetri=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonDetri[,c(5,6,8,9)]=NA
doublonDetri[,7]="Detritus"
doublonMacrop=tab_releve[!duplicated(tab_releve$Id_Sites),]
doublonMacrop[,c(5,6,8,9)]=NA
doublonMacrop[,7]="Macrophyte"
## Add these taxa to the dataset
tab_releve=rbind(tab_releve,doublonPhyto,doublonZoo,doublonDetri,doublonMacrop)

rm(doublonZoo,doublonPhyto,doublonDetri,doublonMacrop)

tab_releve$Id_Sites=as.factor(tab_releve$Id_Sites)
length(levels(tab_releve$Id_Sites))
length(levels(droplevels(tab_releve[tab_releve$Id_Sites %in% levels(Surf$Id_Mare),2])))


dist=dist[dist$Id_Mare %in% levels(droplevels(tab_releve[tab_releve$Id_Sites %in% levels(Surf$Id_Mare),2])),]
Surf=Surf[Surf$Id_Mare %in% levels(droplevels(tab_releve[tab_releve$Id_Sites %in% levels(Surf$Id_Mare),2])),]

envi=merge(Surf,dist,by.x = c("Id_Mare"),by.y = c("Id_Mare"))
envi=envi[,c("Distance","Surface")]

####


x=Matrix_Reg_Manual
y=envi_data
t=100
name_tab="data"
repeti=2
require("enaR")
require("igraph")
require("NetIndices")


setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
Matrix_Reg_Manual=read.delim2("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/Matrix_Final.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
envi_data=read.table("envi_data.txt",header=T)
tab_releve=read.delim("tab_releve.txt")
dist=read.delim2("Distance.txt")
colnames(dist)[2]="Distance"
####
source("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/matrix_function.R")
####
data_format(tab_releve,envi_data,dist)

####


x=Matrix_Reg_Manual
y=envi_data
t=100
name_tab="data"
repeti=2

simu_TTIB_Rep=function(x,y,t=200,name_tab="data",repeti=100){
  output <- NULL
  pb <- txtProgressBar(min = 0, max = repeti, style = 3)
  start.time = Sys.time()
  name_site=as.character(y$Mares)
  DisP=rescale(envi[,1], c(0.95,0.05))
  SurfP=rescale(envi[,2], c(0.95,0.05))
  P=length(x[,1])
  Base=c("Detritus","Macrophyte","Phytoplankton","Zooplankton")
  basal=c("Detritus","Macrophyte","Phytoplankton")
  
  for (repet in 1:repeti) {
    
    Sys.sleep(0.1)		
    setTxtProgressBar(pb, repet)
    
    for (mares in 1:length(DisP)) {
      
      
      temp=matrix(0,length(x[,1]),t)
      row.names(temp)=row.names(x)
      
      temp[Base,1]=1
      
      for (i in 2:t) {
        Colonization=0
        Extinction=0
        C=DisP[mares]*(1-(sum(temp[,i-1]))/P)
        E=SurfP[mares]*(sum(temp[,i-1])/P)
        
        if(i>1){
          temp[,i]=temp[,i-1]
          present=names(temp[temp[,i-1]==1,i-1])
          possible=colnames(x[present,colSums(x[present,])!=0])
          possible=possible[!(possible %in% present)]
          Colonization=sample(c(0,1),1,prob = c(1-C,C))
          Extinction=sample(c(0,1),1,prob = c(1-E,E))
          
          if(Colonization==1){
            add=sample(possible,1)
            temp[add,i]=1
          }
          if(Extinction==1 & length(present[!(present %in% Base)])!=0){
            rmv=sample(present[!(present %in% Base)],1)
            temp[rmv,i]=0
            cascade=TRUE
            while (cascade) {
              presenti=names(temp[temp[,i]==1,i])
              temp_matrix=x[presenti,presenti]
              legal=names(colSums(temp_matrix)>0)
              illegal=presenti[!(presenti %in% legal)]
              if(length(illegal)!=0){
                temp[illegal,i]=0
              }
              else{cascade=FALSE}
              
            }
          }
          
        }
        
      }
      present=names(temp[temp[,i]==1,i])
      
      temp_matrix=x[present,present]
      test=present[!(present %in% basal)]
      temp2_matrix=x[present,test]
      
      
      
      g=graph_from_adjacency_matrix(temp_matrix, mode="undirected")
      
      
      
      C <- sum(temp_matrix)/((dim(temp_matrix)[1])*(dim(temp_matrix)[1]-1)) #Connectance if we prevent cannibalism
      mfcl <- 0
      mfcls <- mfcl # Mean food chain length
      links <- sum(temp_matrix) # Number of links
      modularity_t=modularity(g,membership(cluster_walktrap(g)))
      Mean_Trophic_Level=mean(TrophInd(temp_matrix)[,1])
      Mean_Omnivory=mean(TrophInd(temp_matrix)[,2])
      links.species <- sum(temp_matrix)/length(present) # Number of links per species
      indegree <- MeanGenerality(temp_matrix)
      outdegree <- MeanVulnerability(temp_matrix)
      basal <- FractionOfBasal(temp_matrix) # Percentage of basal species
      top <- FractionOfTop(temp_matrix) # Percentage of top species
      intermediate <- FractionOfIntermediate(temp_matrix) # Percentage of intermediate species
      omnivory <- Omnivory(temp_matrix)
      S_top <- NumberOfTop(temp_matrix)# Number of top species
      S_intermediate <- NumberOfIntermediate(temp_matrix) # Number of intermediate species
      S_basal <- NumberOfBasal(temp_matrix) # Number of basal species
      ratio <- DisP[mares]/SurfP[mares] # Area
      Distance=envi[mares,1]
      Surface=envi[mares,2]
      
      colonization_rate <- DisP[mares]
      extinction_rate <- SurfP[mares]
      
      
      cur_out <- data.frame(replicate=repeti, S=length(present), C, links,modularity_t,
                            Mean_Trophic_Level,Mean_Omnivory,links.species, 
                            indegree, outdegree, basal, top, intermediate, omnivory, S_top, 
                            S_intermediate, S_basal, ratio, colonization_rate,
                            extinction_rate,Distance,Surface)
      
      
      if(is.null(output)){
        output <- cur_out
      }else{
        output <- rbind(output, cur_out)
      } 
      
    }
    
  }
  assign(paste(name_tab),
         output, env = .GlobalEnv) 
  close(pb)
  close(pb) # Close the progress bar
  end.time <- Sys.time() #Establish the time when ending the function
  time.taken <- end.time - start.time # Calculate the time the function worked
  time.taken 
}



simu_TTIB_Rep(Matrix_Reg_Manual,envi_data,t=1000,name_tab="data",repeti = 100)



write.table(data,"simulation_TTIB_178.txt", sep="\t")
