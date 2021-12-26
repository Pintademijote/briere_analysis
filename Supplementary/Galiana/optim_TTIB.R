setwd('C:/Users/pglemasle/Downloads/Code_Galiana_EtAl.2018')
source('utils_NAR.R') #This code is used to compute all the network metrics

Matrix_Reg_Manual=read.delim2("S:/2018/pglemasle/Data/Script/Matrix_Manual_Chekck.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
library(ade4)
require("enaR")
require("igraph")
require("NetIndices")


setwd("S:/2018/pglemasle/Data/Script")
Matrix_Reg_Manual=read.delim2("Matrix_Manual_Chekck.txt", row.names = 1)
colnames(Matrix_Reg_Manual)=row.names(Matrix_Reg_Manual)
Matrix_Reg_Manual=as.matrix(Matrix_Reg_Manual)
envi_data=read.table("envi_data.txt",header=T)
tab_releve=read.delim("tab_releve.txt")
dist=read.delim2("Distance.txt")
colnames(dist)[2]="Distance"
####
source("S:/2018/pglemasle/Data/Script/matrix_function.R")
####
data_format(tab_releve,envi_data,dist)

create_matrix_marre(tab_releve,Matrix_Reg_Manual,envi_data,Count = T, year=F,name_tab = "data")




TTIB=function(par){
  
timesteps <- 1000
output <- NULL
c=par[1]
e=par[2]
square=rep(NA,100)
for (re in 1:100) {
  
      food_web=Matrix_Reg_Manual # This loads the networks previously generated with the niche model that will be used as the metaweb
      # For this to work you need to save the networks generated in a folder called networks within your working directory.
      
      S <- length(Matrix_Reg_Manual[,1])
      presences <- array(0, S)
      
      for(i in 1:timesteps){ # This loop runs the colonisation-extinction dynamics with the trophic constraint of the TTIB
        #basal species
        basal_sps <- which(colSums(food_web) == 0)
        cur_com <- presences
        
        #this is what species are present
        present <- which(cur_com != 0)
        
        if(sum(cur_com) != 0){
          
          #this is the set of species in the network that have prey in the current community
          with_prey <- which(colSums(as.matrix(food_web[present, ])) > 0)
          #this are the potential colonisers
          potential <- append(basal_sps , with_prey)
          #to which we remove the ones that are already present
          potential <- setdiff(potential, present)
          
        }else{
          potential <- basal_sps
        }
        
        if(length(potential) != 0){
          colonisers <- potential[which(runif(length(potential), 0, 1) < c)]
          presences[colonisers] <- 1
        }
        
        extinctions <- FALSE
        if(length(present) != 0){
          extinct <- present[which(runif(length(present), 0, 1) < e)]
          if(length(extinct) != 0){
            presences[extinct] <- 0
            extinctions <- TRUE
          }
        }
        
        if(extinctions){
          cur_com <- presences
          #this is what species are present
          present <- which(cur_com != 0)
          #this is the set of species in the network that have prey in the current community
          with_prey <- which(colSums(as.matrix(food_web[present, ])) > 0)
          legal <- union(with_prey, basal_sps)
          illegal <- present[!(present %in% legal)]
          presences[illegal] <- 0  
        }
        
      }
      
      results <- presences
      net <- food_web
      
      current_species <- which(results > 0) # Number of species
      
      if(length(current_species) < 2){
        next
      }
      
      M <- food_web[current_species, current_species] # This is the assembled community after the colonization-extinction dynamics
      
      #we can check if the network is connected
      graf <- igraph::graph.adjacency(M); # Making sure the network is connected
      y <- igraph::no.clusters(graf, mode='weak');
      attempts <- 0
      if(y > 1 ){
        print('network not found')
        square[re]=NA
        next
      }

      
 
      square[re]=(datacount[1,2]-length(current_species))^2
}
return(sum(square,na.rm = T))
}
  

TTIB(c(0.1,0.2))

optim(par = c(0.1, 0.1), TTIB)

TTIB(c(0.09412248,0.07003506))
