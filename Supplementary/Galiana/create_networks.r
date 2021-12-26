


#### create niche networks
source('NicheNetwork.R') # code for the niche model (Williams and Martinez, 2000)

n_networks <- 1 # Number of networks that you want to generate using the niche model.

# In the niche model you need to define the number of species of the network and its connectance

S <- 200 # Number of species
C <- connectance(S) # Connectance of the network; Here is defined based on a scaling relationship 
                    # with the number of species defined in the code NicheNetwork.R that you can modify.


for(n in 1:n_networks){
  food_web <- NicheNetwork(S,C)
  save(food_web, file=paste('network-',n, sep=''))
}
