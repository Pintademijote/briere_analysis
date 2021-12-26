setwd('./Code_Galiana_EtAl.2018')

source('utils_NAR.R') #This code is used to compute all the network metrics
source('space.r') # This code is used to define the spatial structure of the meta-community. It defines the connections between the patches based on the distance between them and the threshold for connectivity defined.

#Defining the structure of the meta-community.
islands <- 10 #Number of patches/islands you want
#S <- 50
#C <- connectance(S)
landscape <- space(islands, 0.3) #the function 'space' takes as inputs the number of patches and the threshold for connectivity between them desired. 
landscape[landscape > 0] <- 1

#parameters
c <- 0.1    #colonisation rate
e <- 0.4    #extinction rate
dispersals <- c(0, .01, .05, .1, .2, .3, .4, .5) #dispersal rates


timesteps <- 1000 # Number of time steps used in the colonisation-extinction dynamics.
replicates <- 1 # Number of replicates desired for the simulations. In the paper we used 100.
max_attempts <- 100

output <- NULL
for(d in dispersals){
  print(d)
  for(j in 1:replicates){
    load(paste('./networks/network-',j,sep=''))# This loads the networks previously generated with the niche model that will be used as the metaweb
    # For this to work you need to save the networks generated in a folder called networks within your working directory.
    S <- food_web$S
    presences <- array(0, c(islands,S))
    
    print(j)
    for(i in 1:timesteps){
      #basal species
      basal_sps <- which(colSums(food_web$M) == 0)
      for(isl in 1:islands){
        #this is the current community
        cur_com <- presences[isl,]
        #this is what species are present
        present <- which(cur_com != 0)
        #this is the set of species in the network that have prey in the current community
        with_prey <- which(colSums(as.matrix(food_web$M[present, ])) > 0)
        #this are the potential colonisers
        potential <- append(basal_sps , with_prey)
        #to which we remove the ones that are already present
        potential <- setdiff(potential, present)
        
        if(length(potential) != 0){
          colonisers <- potential[which(runif(length(potential), 0, 1) < c)]
          presences[isl, colonisers] <- 1
        }
        
        extinctions <- FALSE
        if(length(present) != 0){
          extinct <- present[which(runif(length(present), 0, 1) < e)]
          if(length(extinct) != 0){
            presences[isl, extinct] <- 0
            extinctions <- TRUE
          }
        }
        
        if(extinctions){
          cur_com <- presences[isl,]
          #this is what species are present
          present <- which(cur_com != 0)
          #this is the set of species in the network that have prey in the current community
          with_prey <- which(colSums(as.matrix(food_web$M[present, ])) > 0)
          legal <- union(with_prey, basal_sps)
          illegal <- present[!(present %in% legal)]
          presences[isl, illegal] <- 0  
        }
        
      }
      
      
      ### dispersal after colonisation and extinction processes
      
      randCol = matrix(runif(islands*S,0,1),nr=islands,nc=S)
      ColMat = matrix(0,nr=islands,nc=S)
      
      ConPop = landscape%*%presences
      LocalCol = 1-(1-d)^ConPop
      
      LocalPrey = presences %*% food_web$M
      LocalPrey[LocalPrey>0] = 1
      LocalPrey[,basal_sps] = 1
      
      ColMat[presences == 0 & LocalPrey == 1 & randCol<LocalCol] = 1
      
      presences = presences + ColMat
      
    }
    
    
    
    results <- presences # This is the community matrix presences
    net <- food_web
    
    #### after obtaining the results we aggregate local communities and obtain network attributes 
    species <- c()
    area <-c()
    C <- c()
    links <- c()
    links.species <- c()
    mfcls <- c()
    indegree <- c()
    outdegree <- c()
    basal <- c()
    top <- c()
    intermediate <- c()
    omnivory <- c()
    S_top <- c()
    S_intermediate <- c()
    S_basal <- c()
    
    for(i in 1:islands){
      areas <- i
      # This is the aggregation process to generate the increase in area. 
      # In this model area correspond to the number of patches that have been aggregated. The order of aggreggation in this case is random.
      patches_selected <- as.matrix(results[sample(nrow(results), i),])
      
      if(i > 1){
        current_species_ran <- which(colSums(patches_selected[,]) > 0)
        
      }else{
        current_species_ran <- which(patches_selected[,] > 0)
        
      }
      if(length(unique(current_species_ran)) <= 1){
        next
      }
      
      M <- net$M[current_species_ran, current_species_ran]
      
      #we can check if the network is connected
      graf <- igraph::graph.adjacency(M); # Making sure the network is connected
      y <- igraph::no.clusters(graf, mode='weak');
      attempts <- 0
      while(y > 1 & attempts < max_attempts){
        patches_selected <- as.matrix(results[sample(nrow(results), i),])
        
        if(i > 1){
          current_species_ran <- which(colSums(patches_selected[,]) > 0)
          
        }else{
          current_species_ran <- which(patches_selected[,] > 0)
          print(i)
        }    
        M <- net$M[current_species_ran, current_species_ran]
        graf <- igraph::graph.adjacency(M); # Making sure the network is connected
        y <- igraph::no.clusters(graf, mode='weak');
        
        attempts <- attempts + 1
      }
      
      if(length(unique(current_species_ran)) < 2 | attempts >= max_attempts){
        print('network not found')
        next
      }
      
      species <- append(species, length(current_species_ran))
      C <- append(C, sum(M)/((dim(M)[1])*(dim(M)[1]-1))) #Connectance if we prevent cannibalism
      
      mfcl <- tryCatch({
        MeanFoodChainLength(M)
      }, warning = function(w) {
        'NA'
      }, error = function(e) {
        'NA'
      }, finally = {
        
      })
      
      mfcls <- append(mfcls, mfcl) # Mean food chain length
      area <- append(area, areas)
      links <- append(links, sum(M)) # Number of Links
      links.species <- append(links.species, sum(M)/dim(M)[1]) # Number of links per species
      indegree <- append(indegree, MeanGenerality(M)) 
      outdegree <- append(outdegree, MeanVulnerability(M))
      basal <- append(basal, FractionOfBasal(M)) # Percentage of basal species
      top <- append(top, FractionOfTop(M)) # Percentage of top species
      intermediate <- append(intermediate, FractionOfIntermediate(M)) # Percentage of intermediate species
      omnivory <- append(omnivory, Omnivory(M))
      S_top <- append(S_top, NumberOfTop(M)) # Number of top species
      S_intermediate <- append(S_intermediate, NumberOfIntermediate(M)) # Number of intermediate species
      S_basal <- append(S_basal, NumberOfBasal(M)) # Number of basal species
      
      cur_out <- data.frame(rep(j,length(species)),rep(d,length(species)), species, C, links, links.species, area, mfcls, indegree, outdegree, basal, top, intermediate, omnivory, S_top, S_intermediate, S_basal)
      
    }
      
    if(is.null(output)){
      output <- cur_out
    }else{
      output <- rbind(output, cur_out)
    } 
  
  }
  
}

colnames(output)[2]<-"dispersal"
write.csv(output, file='output_metacom_test.csv') # Save the output


### Plotting Network-Area relationships; In this model area correspond to the number of patches aggreggated.

ggplot(data=output, aes(area, species, colour=as.factor(output$dispersal))) + geom_smooth(linetype=1, size=1,fullrange=FALSE)


