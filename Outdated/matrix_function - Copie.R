
####CrÃ©ation de la matrice rÃ©gionnale

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



###Passage des relevÃÂ©es aux matrices locales


create_matrix_marre=function(x,y,z,Matrix=F,descriptor=T,Count=F, year=F, name_tab = "data"){
  
  start.time = Sys.time() #Establish the time when running the function
  
  if(year==F){
    marres=levels(x$Id_Sites) #Assign in an object with the name of the site
    
    
    pb <- txtProgressBar(min = 0, max = length(marres), style = 3) #Create the progress bar
    
    temp_count=matrix(nrow = length(marres),ncol = length(levels(x$Taxon))+1)
    colnames(temp_count)=c(as.character(as.data.frame(table(x[x[,2],7]))[,1]),"Site")
    
    descriptors=as.data.frame(matrix(nrow=length(levels(x$Id_Sites)),ncol=8)) #Create an empty data frame for the descriptors
    colnames(descriptors)=c("Site","n","L","C","Modularity","Mean_trophic_level","Mean_Omnivory_Index", "Richesse_spe") #Specified the name of the column for the descriptos object
    name_site=levels(x$Id_Sites) #Assign in an object with the name of the site
    name_matrix=1:length(levels(x$Id_Sites)) #Assign a numeric level for the name of matrix
    
    for (o in 1:length(marres)) {#Loop through each ponds
      Sys.sleep(0.1) # Define refresh rate for progress bar
      setTxtProgressBar(pb, o) #Set the current status of the progress bar
      
      temp_count[o,1:length(levels(x$Taxon))]=as.data.frame(table(x[x[,2]==marres[o],7]))[,2]
      temp_count[o,length(levels(x$Taxon))+1]=name_site[o]
      
      temp_site=x[x[,2]==marres[o],] #Subset the taxon occurrence data for a define pond
      temp_site=droplevels(temp_site) #Forget old levels
      
      
      nm=paste("Matrix_site_",marres[o],sep="") #Define the name of the matrix returned for the current pond
      
      

      level=as.character(c(levels(temp_site$Taxon))) #Take species species column name of species present in ponds
      
      temp=y[level,level] #subset the regionnal network with only present species
      
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
    
    if(descriptor==T & Count==F){assign(paste(name_tab),
                                        merge(descriptors,z,by.x = c("Site"),by.y = c("Mares")), env = .GlobalEnv)} # If "descriptor" is true, put in the global environnement the data frame of descriptors
    col_temp=temp_count[,length(levels(x$Taxon))+1]
    temp_count[temp_count>0]=1
    temp_count[,length(levels(x$Taxon))+1]=col_temp
    if(Count==T){assign(paste(name_tab,"count",sep=""),
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
          
          temp_count[o,1:length(levels(x$Taxon))]=as.data.frame(table(x[x[,2]==marres[o] & x[,3]==years[h],7]))[,2]
          temp_count[o,length(levels(x$Taxon))+1]=name_site[o]
          
          temp_site=x[x[,2]==marres[o] & x[,3]==years[h],] #Subset the taxon occurrence data for a define pond
          temp_site[,2]=droplevels(temp_site[,2]) #Forget old levels
          temp_site[,7]=droplevels(temp_site[,7]) #Forget old levels
          
          nm=paste("Matrix_site_",marres[o],"_",years[h],sep="") #Define the name of the matrix returned for the current pond
          
          
          level=as.character(c(levels(temp_site$Taxon)))
          
          
          if(length(level)==1){level=as.character(c(levels(temp_site$Taxon),"Zooplankton","Phytoplankton","Detritus","Macrophyte"))}
          
          temp=y[level,level]
          
          
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
        
        
        if(descriptor==T & Count==F){assign(paste(name_tab),
                                            merge(descriptors,envi_data,by.x = c("Site"),by.y = c("Mares"),all.x = T), env = .GlobalEnv)} # If "descriptor" is true, put in the global environnement the data frame of descriptors
        col_temp=temp_count[,length(levels(x$Taxon))+1]
        temp_count[temp_count>0]=1
        temp_count[,length(levels(x$Taxon))+1]=col_temp
        if(Count==T){assign(paste(name_tab,"count",sep = ""),
                            merge(merge(descriptors,envi_data,by.x = c("Site"),by.y = c("Mares")),temp_count,by.x = c("Site"),by.y = c("Site")), env = .GlobalEnv)}
        
      }
    }
  }
  
  
  
  close(pb) # Close the progress bar
  end.time <- Sys.time() #Establish the time when ending the function
  time.taken <- end.time - start.time # Calculate the time the function worked
  time.taken # Return the time taking for the function to work
  
}


#####


data_format=function(x,y,z){
  #############
  ## Keep only ponds
  x=x[x[,1]=="Mare",]
  
  ## Reset factor levels
  x=droplevels(x)
  ## Reset format
  x[,2]=as.character(x[,2], trim = TRUE)
  x[,2]=as.factor(x[,2])
  
  
  #####
  ## Standardize site IDs to fit with with x
  y[,2]=gsub("M_", "", y[,2], fixed = TRUE)
  y[,2]=sub("^0+", "", y[,2])
  y[,2]=factor(y[,2])
  
  #####
  
  
  #### Remove sites not present in envi-data to shorten the loop
  test=as.factor(x[!duplicated(x[,2]),2])	# Keep only one line per site
  test=as.factor(x[,2])
  ####
  x=x[(test %in% y[,2]),]
  x=droplevels(x)
  rm(test)
  
  ####
  ##
  
  z[,1]=gsub("Mare_", "", z[,1], fixed = TRUE)
  z[,1]=sub("^0+", "", z[,1])
  z[,1]=as.factor(z[,1])
  z=droplevels(z)
  
  
  #####
  ##
  #y=merge(y,x[!duplicated(x$Id_Sites),c(2,15)],by.x = c("Mares"),by.y = c("Id_Sites"))
  y=merge(y,z,by.x = c("Mares"),by.y = c("Id_Mare"))
  #### ADD non-observed taxa at all sites
  doublonPhyto=x[!duplicated(x$Id_Sites),]
  doublonPhyto[,c(5,6,8,9)]=NA				# Remove useless data but keep matrix size		
  doublonPhyto[,7]="Phytoplankton"
  doublonZoo=x[!duplicated(x$Id_Sites),]
  doublonZoo[,c(5,6,8,9)]=NA
  doublonZoo[,7]="Zooplankton"
  doublonDetri=x[!duplicated(x$Id_Sites),]
  doublonDetri[,c(5,6,8,9)]=NA
  doublonDetri[,7]="Detritus"
  doublonMacrop=x[!duplicated(x$Id_Sites),]
  doublonMacrop[,c(5,6,8,9)]=NA
  doublonMacrop[,7]="Macrophyte"
  ## Add these taxa to the dataset
  x=rbind(x,doublonPhyto,doublonZoo,doublonDetri,doublonMacrop)
  
  rm(doublonZoo,doublonPhyto,doublonDetri,doublonMacrop)
  
  assign("tab_releve",x, env = .GlobalEnv)
  assign("envi_data",y, env = .GlobalEnv)
}




