setwd("S:/2018/pglemasle/Data/Script/")
#setwd("~/Documents/Stages_theses/2018/M1_reseaux_trophiques/Data")

list.of.packages <- c("enaR", "igraph","NetIndices") 
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages) 
lapply(new.packages, require, character.only=T)
rm(list.of.packages,new.packages)

# Import trophic link database
tab_trophic=read.delim("tab_trophic.txt")
colnames(tab_trophic)[1]="Predator"
# Import taxon occurrence data
tab_releve=read.delim("tab_releve.txt")
# Import environmental variables
envi_data=read.delim2("envi_data.txt")
# Import distance from marsh
dist=read.delim2("Distance.txt")



#############
## Keep only ponds
tab_releve=tab_releve[tab_releve[,1]=="Mare",]

## Reset factor levels
tab_releve=droplevels(tab_releve)
## Reset format
tab_releve[,2]=as.character(tab_releve[,2], trim = TRUE)
tab_releve[,2]=as.factor(tab_releve[,2])


#####
## Standardize site IDs to fit with with tab_releve
envi_data[,2]=gsub("M_", "", envi_data[,2], fixed = TRUE)
envi_data[,2]=sub("^0+", "", envi_data[,2])
envi_data[,2]=factor(envi_data[,2])

#####


#### Remove sites not present in envi-data to shorten the loop
test=as.factor(tab_releve[!duplicated(tab_releve[,2]),2])	# Keep only one line per site
test=as.factor(tab_releve[,2])
####
tab_releve=tab_releve[(test %in% envi_data[,2]),]
tab_releve=droplevels(tab_releve)
rm(test)

####
##

dist[,1]=gsub("Mare_", "", dist[,1], fixed = TRUE)
dist[,1]=sub("^0+", "", dist[,1])
dist[,1]=factor(dist[,1])

####
##



#####
##
envi_data=merge(envi_data,tab_releve[!duplicated(tab_releve$Id_Sites),c(2,15)],by.x = c("Mares"),by.y = c("Id_Sites"))
envi_data=merge(envi_data,dist,by.x = c("Mares"),by.y = c("Id_Mare"))
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

####Création de la matrice régionnale

create_matrix=function(x){
  # Name of the matrix
  nm="Matrix_tab_trophic"
  # 0-full matrix of the desired size
  temp=matrix(0,length(levels(x$Predator)),length(levels(x$Predator)))
  # Assign names for rows and columns of the 0-full matrix
  colnames(temp)=levels(x$Predator)
  row.names(temp)=levels(x$Predator)
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = length(levels(x$Predator)), style = 3)
  # Start building the matrix

    for (j in 1:length(levels(x$Predator))) {	# Loop through each 
      Sys.sleep(0.1)		# Refreshing rate of progress bar
      setTxtProgressBar(pb, j)	# Update progress bar at each iteration
      for (k in 1:length(levels(x$Predator))) {
        if(row.names(temp)[k] %in% x[x[,1]==colnames(temp)[j],2]){
          temp[k,j]=1
        }
        
      }
      
    }
  close(pb)
  assign(nm,temp, env = .GlobalEnv)
  }



create_matrix(tab_trophic)

###Passage des relevÃ©es aux matrices locales

x=tab_releve
y=tab_trophic

create_matrix_marre=function(x,y,z,Matrix=F,descriptor=T,Count=F, year=F, name_tab = "Descriptors"){
  #x is taxon occurrence data
  #y is trophic link database
  #if true arg matrix specified if it will return the local matrix
  #if true descriptor specified if it will return the table with descriptors of local matrix
  
  start.time = Sys.time() #Establish the time when running the function
  
  
  if(year==F){
  marres=levels(x$Id_Sites) #Assign in an object with the name of the site
  
  
  pb <- txtProgressBar(min = 0, max = length(marres), style = 3) #Create the progress bar
  
  temp_count=matrix(nrow = length(marres),ncol = length(levels(x$Taxon))+1)
  colnames(temp_count)=c(as.character(as.data.frame(table(tab_releve[tab_releve[,2],7]))[,1]),"Site")
  
  descriptors=as.data.frame(matrix(nrow=length(levels(x$Id_Sites)),ncol=8)) #Create an empty data frame for the descriptors
  colnames(descriptors)=c("Site","n","L","C","Modularity","Mean_trophic_level","Mean_Omnivory_Index", "Richesse_spe") #Specified the name of the column for the descriptos object
  name_site=levels(x$Id_Sites) #Assign in an object with the name of the site
  name_matrix=1:length(levels(x$Id_Sites)) #Assign a numeric level for the name of matrix
  
  for (o in 1:length(marres)) {#Loop through each ponds
    Sys.sleep(0.1) # Define refresh rate for progress bar
    setTxtProgressBar(pb, o) #Set the current status of the progress bar
    
    temp_count[o,1:62]=as.data.frame(table(x[x[,2]==marres[o],7]))[,2]
    temp_count[o,63]=name_site[o]
    
    temp_site=x[x[,2]==marres[o],] #Subset the taxon occurrence data for a define pond
    temp_site[,2]=droplevels(temp_site[,2]) #Forget old levels
    temp_site[,7]=droplevels(temp_site[,7]) #Forget old levels
    
    nm=paste("Matrix_site_",marres[o],sep="") #Define the name of the matrix returned for the current pond
    
    
    temp=matrix(0,length(levels(temp_site$Taxon)),length(levels(temp_site$Taxon))) # 0-full matrix of the desired size
    
    colnames(temp)=levels(temp_site$Taxon) # Define species name for columns of the temporary matrix
    row.names(temp)=levels(temp_site$Taxon) # Define species name for rows of the temporary matrix
    
    
      for (j in 1:length(levels(temp_site$Taxon))) { # Loop through columns of the temp matrix
        for (k in 1:length(levels(temp_site$Taxon))) { # Loop through rows of the temp matrix
          if(row.names(temp)[k] %in% y[y[,1]==colnames(temp)[j],2]){ # If for the current position in the matrix,
            #the specie in current row is present in the subset taxon occurrence data for the predator in the current column: 
            temp[k,j]=1 # Add 1 in the current cell check in the temporay matrix
          }
          
        }
        
      }
    if(descriptor==T){ # If arg "descriptor" is true, calculate the descriptor
    g=graph_from_adjacency_matrix(temp, mode="undirected")
    descriptors[o,1]=name_site[o]
    name_matrix[o]=paste("Matrix_site_",name_site[o],sep="")
    descriptors[o,2:4]=enaStructure(temp)$ns[,c(1,2,3)]
    descriptors[o,5]=modularity(g,membership(cluster_walktrap(g)))
    descriptors[o,6]=mean(TrophInd(temp)[,1])
    descriptors[o,7]=mean(TrophInd(temp)[,2])
    descriptors[o,8]=length(temp[,1])
    }
    if(Matrix==T){assign(nm,temp, env = .GlobalEnv)} # If arg "Matrix" is true, put in global environnement adjancy matrix of the current pond
    }
    
  descriptors$Site=as.factor(descriptors$Site)
  descriptors$Anne=as.factor(descriptors$Anne)
    
  if(descriptor==T & Count==F){assign("Descriptors",
                           merge(descriptors,z,by.x = c("Site"),by.y = c("Mares")), env = .GlobalEnv)} # If "descriptor" is true, put in the global environnement the data frame of descriptors
  col_temp=temp_count[,length(levels(x$Taxon))+1]
  temp_count[temp_count>0]=1
  temp_count[,length(levels(x$Taxon))+1]=col_temp
  if(Count==T){assign("Descriptors_Count",
                     merge(merge(descriptors,z,by.x = c("Site"),by.y = c("Mares")),temp_count,by.x = c("Site"),by.y = c("Site")), env = .GlobalEnv)}
  }
  
  if(year==T){
    marres=levels(x$Id_Sites) #Assign in an object with the name of the site
    years=levels(x$Annee)
    
    pb <- txtProgressBar(min = 0, max = length(marres), style = 3) #Create the progress bar
    
    temp_count=matrix(nrow = length(marres),ncol = length(levels(x$Taxon))+1)
    colnames(temp_count)=c(as.character(as.data.frame(table(tab_releve[tab_releve[,2],7]))[,1]),"Site")
    
    count_year=as.matrix(table(x[,c(2,3)]))[,1:6]
    count_year[count_year>0]=1
    
    
    
    descriptors=as.data.frame(matrix(nrow=sum(count_year),ncol=9)) #Create an empty data frame for the descriptors
    colnames(descriptors)=c("Site","n","L","C","Modularity","Mean_trophic_level","Mean_Omnivory_Index", "Richesse_spe","Annee") #Specified the name of the column for the descriptos object
    name_site=levels(x$Id_Sites) #Assign in an object with the name of the site
    name_matrix=1:length(levels(x$Id_Sites)) #Assign a numeric level for the name of matrix
    year_matrix=1:length(levels(x$Annee))
    
    for (o in 1:length(marres)) {#Loop through each ponds
      Sys.sleep(0.1) # Define refresh rate for progress bar
      setTxtProgressBar(pb, o) #Set the current status of the progress bar
      
      for(h in 1:length(levels(x$Annee))){
        if(years[h] %in% x[x[,2]==marres[o],3]){
          
          temp_count[o,1:62]=as.data.frame(table(x[x[,2]==marres[o] & x[,3]==years[h],7]))[,2]
          temp_count[o,63]=name_site[o]
          
          temp_site=x[x[,2]==marres[o] & x[,3]==years[h],] #Subset the taxon occurrence data for a define pond
          temp_site[,2]=droplevels(temp_site[,2]) #Forget old levels
          temp_site[,7]=droplevels(temp_site[,7]) #Forget old levels
          
          nm=paste("Matrix_site_",marres[o],"_",years[h],sep="") #Define the name of the matrix returned for the current pond
          
          
          temp=matrix(0,length(levels(temp_site$Taxon)),length(levels(temp_site$Taxon))) # 0-full matrix of the desired size
          
          colnames(temp)=levels(temp_site$Taxon) # Define species name for columns of the temporary matrix
          row.names(temp)=levels(temp_site$Taxon) # Define species name for rows of the temporary matrix
          
          
          for (j in 1:length(levels(temp_site$Taxon))) { # Loop through columns of the temp matrix
            for (k in 1:length(levels(temp_site$Taxon))) { # Loop through rows of the temp matrix
              if(row.names(temp)[k] %in% y[y[,1]==colnames(temp)[j],2]){ # If for the current position in the matrix,
                #the specie in current row is present in the subset taxon occurrence data for the predator in the current column: 
                temp[k,j]=1 # Add 1 in the current cell check in the temporay matrix
              }
              
            }
            
          }
          if(descriptor==T){ # If arg "descriptor" is true, calculate the descriptor
            g=graph_from_adjacency_matrix(temp, mode="undirected")
            descriptors[o,1]=name_site[o]
            name_matrix[o]=paste("Matrix_site_",name_site[o],sep="")
            descriptors[o,2:4]=enaStructure(temp)$ns[,c(1,2,3)]
            descriptors[o,5]=modularity(g,membership(cluster_walktrap(g)))
            descriptors[o,6]=mean(TrophInd(temp)[,1])
            descriptors[o,7]=mean(TrophInd(temp)[,2])
            descriptors[o,8]=length(temp[,1])
            descriptors[o,9]=years[h]
          }
          if(Matrix==T){assign(nm,temp, env = .GlobalEnv)} # If arg "Matrix" is true, put in global environnement adjancy matrix of the current pond
        }
        
        descriptors$Site=factor(descriptors$Site)
        descriptors$Anne=factor(descriptors$Anne)
        
        if(descriptor==T & Count==F){assign(paste(name_tab),
                                            merge(descriptors,envi_data,by.x = c("Site"),by.y = c("Mares"),all.x = T), env = .GlobalEnv)} # If "descriptor" is true, put in the global environnement the data frame of descriptors
        col_temp=temp_count[,length(levels(x$Taxon))+1]
        temp_count[temp_count>0]=1
        temp_count[,length(levels(x$Taxon))+1]=col_temp
        if(Count==T){assign(paste(name_tab,"count",sep=""),
                            merge(merge(descriptors,envi_data,by.x = c("Site"),by.y = c("Mares")),temp_count,by.x = c("Site"),by.y = c("Site")), env = .GlobalEnv)}
        
      }
    }
  }
  
  
  
  close(pb) # Close the progress bar
  end.time <- Sys.time() #Establish the time when ending the function
  time.taken <- end.time - start.time # Calculate the time the function worked
  time.taken # Return the time taking for the function to work
}

create_matrix_marre(tab_releve,tab_trophic,envi_data,Count = F, year=T,name_tab = "data")

####


g1=graph_from_adjacency_matrix(Matrix_tab_trophic)
tkplot(g1, canvas.width = 500, canvas.height = 500,vertex.size=4,
       vertex.label.dist=1, vertex.color="red", edge.arrow.size=0.8,label.cex=24,label.color="black")


