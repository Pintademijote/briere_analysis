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
normalized  <- function(x){(x-max(x))/(min(x)-max(x))}

x=Matrix_Reg_Manual
y=envi_data
t=100
name_tab="data"
repeti=2

simu_TTIB_Rep=function(x,y,t=200,name_tab="data",repeti=100){
  pb <- txtProgressBar(min = 0, max = repeti, style = 3)
  start.time = Sys.time()
  descriptors=matrix(NA,nrow = length(y$Mares), ncol = 1+2*6)
  descriptors=as.data.frame(descriptors)
  colnames(descriptors)=c("Site","n_mean","n_sd","L_mean","L_sd","C_mean","C_sd","Modularity_mean",
                          "Modularity_sd","Mean_trophic_level_mean","Mean_trophic_level_sd",
                          "Mean_Omnivory_Index_mean","Mean_Omnivory_Index_sd")
  name_site=as.character(y$Mares)
  DisP=normalized(y$Distance)
  SurfP=normalized(log(y$Surface))
  P=length(x[,1])
  Base=c("Detritus","Macrophyte","Phytoplankton","Zooplankton")
  Nodes=as.data.frame(matrix(NA,nrow = length(y$Mares), ncol = repeti))
  Links=as.data.frame(matrix(NA,nrow = length(y$Mares), ncol = repeti))
  Connectance=as.data.frame(matrix(NA,nrow = length(y$Mares), ncol = repeti))
  Modularity=as.data.frame(matrix(NA,nrow = length(y$Mares), ncol = repeti))
  Mean_Trophic_Level=as.data.frame(matrix(NA,nrow = length(y$Mares), ncol = repeti))
  Mean_Omnivory=as.data.frame(matrix(NA,nrow = length(y$Mares), ncol = repeti))
  countValid=as.data.frame(matrix(NA,nrow = length(y$Mares), ncol = repeti))
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
      test=present[!(present %in% basal)]
      temp2_matrix=x[present,test]
      
      if(class(temp2_matrix)!="integer"){
        if(all(colSums(temp2_matrix)>0)){
          countValid[mares,repet]=1
        }else{countValid[mares,repet]=0}}else{if(sum(temp2_matrix)>0){
          countValid[mares,repet]=1}else{countValid[mares,repet]=0}}
      
      g=graph_from_adjacency_matrix(temp_matrix, mode="undirected")
      
      Nodes[mares,repet]=enaStructure(temp_matrix)$ns[,1]
      Links[mares,repet]=enaStructure(temp_matrix)$ns[,2]
      Connectance[mares,repet]=enaStructure(temp_matrix)$ns[,3]
      Modularity[mares,repet]=modularity(g,membership(cluster_walktrap(g)))
      Mean_Trophic_Level[mares,repet]=mean(TrophInd(temp_matrix)[,1])
      Mean_Omnivory[mares,repet]=mean(TrophInd(temp_matrix)[,2])
      
    }
    
  }
    descriptors[,1]=name_site
    descriptors[,2]=rowMeans(Nodes)
    descriptors[,3]=apply(Nodes,1,sd) 
    descriptors[,4]=rowMeans(Links)
    descriptors[,5]=apply(Links,1,sd) 
    descriptors[,6]=rowMeans(Connectance)
    descriptors[,7]=apply(Connectance,1,sd) 
    descriptors[,8]=rowMeans(Modularity)
    descriptors[,9]=apply(Modularity,1,sd)
    descriptors[,10]=rowMeans(Mean_Trophic_Level)
    descriptors[,11]=apply(Mean_Trophic_Level,1,sd) 
    descriptors[,12]=rowMeans(Mean_Omnivory)
    descriptors[,13]=apply(Mean_Omnivory,1,sd) 
    descriptors[,1]=factor(descriptors[,1])
  assign(paste(name_tab),
         merge(descriptors,y,by.x = c("Site"),by.y = c("Mares")), env = .GlobalEnv) 
  assign("countValid",countValid, env = .GlobalEnv)
  close(pb)
  close(pb) # Close the progress bar
  end.time <- Sys.time() #Establish the time when ending the function
  time.taken <- end.time - start.time # Calculate the time the function worked
  time.taken 
}

simu_TTIB_Rep(Matrix_Reg_Manual,envi_data,t=500,name_tab="data",repeti = 1000)

hist(data$n_mean, breaks = 50)

write.table(data,"simulation_TTIB.txt", sep="\t")
