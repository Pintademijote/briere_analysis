require("enaR")
require("igraph")
require("NetIndices")


setwd("C:/Users/pglem/Documents/Master/Stage M1/Backup_INRA/Data/Script/")
Matrix_Reg_Manual=read.delim2("Matrix_Manual_Chekck.txt", row.names = 1)
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
normalized  <- function(x){(x-min(x))/(max(x)-min(x))}

x=Matrix_Reg_Manual
y=envi_data
t=100
name_tab="data"

simu_TTIB=function(x,y,t=200,name_tab="data"){
  
  start.time = Sys.time()
  descriptors=matrix(NA,nrow = length(y$Mares), ncol = 8)
  descriptors=as.data.frame(descriptors)
  colnames(descriptors)=c("Site","n","L","C","Modularity","Mean_trophic_level","Mean_Omnivory_Index", "Richesse_spe")
  name_site=as.character(y$Mares)
  DisP=normalized(y$Distance)
  SurfP=normalized(log(y$Surface))
  P=length(x[,1])
  Base=c("Detritus","Macrophyte","Phytoplankton","Zooplankton")
  pb <- txtProgressBar(min = 0, max = length(DisP), style = 3)
  
  for (mares in 1:length(DisP)) {
    
    Sys.sleep(0.1)		
    setTxtProgressBar(pb, mares)
    temp=matrix(0,length(x[,1]),t)
    row.names(temp)=row.names(x)
    
    temp[Base,1]=1
    
    for (i in 2:t) {
      Colonization=0
      Extinction=0
      C=(DisP[mares]*(P-sum(temp[,i-1])))/100
      E=(SurfP[mares]*sum(temp[,i-1]))/100
      
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
          presenti=names(temp[temp[,i]==1,i])
          temp_matrix=x[presenti,presenti]
          presenti=presenti[!(presenti %in% Base)]
          temp_matrix=temp_matrix[,presenti]
          if(class(temp_matrix)=="integer"){
            if(sum(temp_matrix)==0){
              temp[presenti,i]=0
            }
          }
          if(class(temp_matrix)=="matrix"){
            if(!all(colSums(temp_matrix)==0)){
              temp[presenti[colSums(temp_matrix)==0],i]=0
          }
          
          }
        }
        
      }
      
    }
    present=names(temp[temp[,i]==1,i])
    
    temp_matrix=x[present,present]
    
    g=graph_from_adjacency_matrix(temp_matrix, mode="undirected")
    descriptors[mares,1]=name_site[mares]
    descriptors[mares,2:4]=enaStructure(temp_matrix)$ns[,c(1,2,3)]
    descriptors[mares,5]=modularity(g,membership(cluster_walktrap(g)))
    descriptors[mares,6]=mean(TrophInd(temp_matrix)[,1])
    descriptors[mares,7]=mean(TrophInd(temp_matrix)[,2])
    descriptors[mares,8]=length(temp_matrix[,1])
    
  }
  descriptors[,1]=factor(descriptors[,1])
  assign(paste(name_tab),
         merge(descriptors,y,by.x = c("Site"),by.y = c("Mares")), env = .GlobalEnv) 
  close(pb)
  close(pb) # Close the progress bar
  end.time <- Sys.time() #Establish the time when ending the function
  time.taken <- end.time - start.time # Calculate the time the function worked
  time.taken 
}
