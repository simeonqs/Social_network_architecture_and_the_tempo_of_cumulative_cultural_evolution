# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Project: Cantor et al. 2020
# Date started: 04-07-2020
# Date last modified: 28-09-2020
# Author: Simeon Quirinus Smeele
# Description: This is the R version of the first agent based model in Cantor et al. 2020. It loads the 
# networks and saves the results from the paths specified in the DATA section. It will install and load 
# required packages. 
# Step 1: Change the paths in the DATA sections to fit your folder location.
# Step 2 (optional): Change the settings in the DATA sections.
# Step 3: Run the code selecting all and clicking 'Run' or line by line by clicking 'Run' multiple times. 
# WARNING: Last line saves the data needed to create the diversity plot. Only save this for a few iterations
# otherwise the file will be very large. 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Loading libraries
libraries = c('data.table', 'tidyverse', 'parallel', 'gridExtra')
for(i in libraries){
  if(! i %in% installed.packages()) lapply(i, install.packages)
  lapply(libraries, require, character.only = TRUE)
}

# Clean R
rm(list=ls()) 
dev.off() # ignore the error message
cat("\014")  

# Paths
path_networks = '/Users/ssmeele/Library/Mobile Documents/com~apple~CloudDocs/Downloads/Social_network_architecture_and_the_tempo_of_cumulative_cultural_evolution-master/code/3_R_agent_based_models/networks'
path_out = '/Users/ssmeele/Library/Mobile Documents/com~apple~CloudDocs/Downloads/Social_network_architecture_and_the_tempo_of_cumulative_cultural_evolution-master/code/3_R_agent_based_models/output/model_1/'

# Settings
n_iter = 3 # number of iterations
n_cores = 4 # number of cores to use (needs to be 1 on Windows)
min_epoch = 50 # how many epochs to run for diversity

# Read in networks
list_networks = list.files(path_networks, full.names = T, pattern = '*pajek', recursive = T) %>% 
  lapply(read.table, skip = 2) %>%
  lapply(function(x){
    l = x %>% unlist %>% max
    mat = matrix(0, l, l)
    mat[as.matrix(x)] = 1
    mat[as.matrix(x[,2:1])] = 1
    mat
  })
names(list_networks) = list.files(path_networks, pattern = '*pajek', recursive = T) %>% 
  str_remove('.pajek') %>% str_split('/') %>% sapply(`[`, 2)

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# ANALYSIS ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Set seed
set.seed(1)

# Small network to debug
network = matrix(c(1,   0,   0.2, 0.4,
                   0,   1,   0.1, 0.1,
                   0.2, 0.1, 1  , 0.3,
                   0.4, 0.1, 0.3, 1), nrow = 4, ncol = 4)
name_network = 'test_d00'

# Plants, lower case are plants, upper case are medicine
plants = list(a1 = 6, a2 = 8, a3 = 10, 
              b1 = 6, b2 = 8, b3 = 10)
combinations = list(A1 = 48, A2 = 109, A3 = 188, AC = 358, 
                    B1 = 48, B2 = 109, B3 = 188, BC = 358)
ingredients = list(A1 = c('a1', 'a2', 'a3'), 
                   A2 = c('A1', 'a1', 'b2'), 
                   A3 = c('A2', 'b2', 'b3'), 
                   AC = c('A3', 'B3', 'A2'), 
                   B1 = c('b1', 'b2', 'b3'), 
                   B2 = c('B1', 'a2', 'a3'), 
                   B3 = c('B2', 'b1', 'a2'), 
                   BC = c('B3', 'A3', 'B2'))

# Function
simCE = function(network, name_network, n_iter){
  
  # Make sure agents names are numeric
  rownames(network) = colnames(network) = 1:nrow(network)
  
  # Include edge to self in network
  diag(network) = 1
  
  # Run through iterations
  out = mclapply(1:n_iter, function(i){
    
    # Equip agents
    agents_pockets = lapply(1:nrow(network), function(x) plants)
    
    # Dataframe to save innovation timings
    timings = data.frame(network_type = character(),
                         degree = numeric(),
                         pop_size = numeric(),
                         iteration = numeric(),
                         epoch = numeric(),
                         timestep = numeric(),
                         innovation_score = numeric(),
                         lineage = character())
    
    # Dataframe to save the proportions
    props = data.frame(network_type = character(),
                       degree = numeric(),
                       pop_size = numeric(),
                       iteration = numeric(),
                       epoch = numeric(),
                       timestep = numeric(),
                       medicin = character(),
                       proportion = numeric())
    
    # Get basics
    network_type = str_split(name_network, '_')[[1]][1]
    degree = str_split(name_network, '_')[[1]][2] %>% str_remove('d') %>% as.integer
    pop_size = nrow(network)
    
    # Running through epochs
    cont = T
    t = 0 # epochs
    tt = 0 # time-steps/interactions
    while(cont | t < min_epoch){ # run until crossover and reached at least 400
      
      # Update t
      t = t+1
      
      # Start with all agents
      agents_left = 1:nrow(network)
      
      # Keep running untill all have been focal
      while(length(agents_left) > 1){
        
        # Update tt
        tt = tt + 1
        
        # Pick random agent from agents left to be focal and remove for next round
        agent_i = sample(as.character(agents_left), 1) %>% as.numeric # otherwise it doesn't work when one left
        agents_left = agents_left[agents_left != agent_i]
        
        # Pick neighbour, make sure it's not the agent itself, weigh by edge-weight
        options = which(network[agent_i,] > 0)
        options = options[options != agent_i]
        weights = network[agent_i, options] / sum(network[agent_i, options])
        agent_j = sample(as.character(options), 1, prob = weights) %>% as.numeric
        
        # Pick ingredients
        pi = as.numeric(agents_pockets[[agent_i]])/sum(as.numeric(agents_pockets[[agent_i]]))
        pj = as.numeric(agents_pockets[[agent_j]])/sum(as.numeric(agents_pockets[[agent_j]]))
        how_many = sample(2)
        ing_choosen = c(sample(names(agents_pockets[[agent_i]]), how_many[1], prob = pi),
                        sample(names(agents_pockets[[agent_j]]), how_many[2], prob = pj))
        
        # Check if invention
        invention = NA
        for(row in 1:length(ingredients)){ # run trough all potential inventions
          if(ingredients[[row]][1] %in% ing_choosen &
             ingredients[[row]][2] %in% ing_choosen &
             ingredients[[row]][3] %in% ing_choosen){
            invention = list(combinations[[row]])       # if all ingredients match, extract value
            names(invention) = names(combinations)[row] # and name
            break
          } # End if-loop
        } # End row loop
        
        # Save if invention
        if(!is.na(invention)){
          
          # Save innovation level
          timings = rbind(timings,
                          data.frame(network_type = network_type,
                                     degree = degree,
                                     pop_size = pop_size,
                                     iteration = i,
                                     epoch = t,
                                     timestep = tt,
                                     innovation_score = as.numeric(invention),
                                     lineage = names(invention)))
          
          # Check if cross-over is found
          if(invention == 358) cont = F

          # Get neighbours to share, include i and j
          n_i = which(network[agent_i,] > 0)
          n_j = which(network[agent_j,] > 0)
          to_receive = unique(c(n_i, n_j))

          # Run through them
          for(agent in to_receive){
            
            if(!names(invention) %in% names(agents_pockets[[agent]])){ # prevent duplicates
              agents_pockets[[agent]] = append(agents_pockets[[agent]], invention)
            }
            
          } # End agent loop
          
        } # End save loop
        
      } # End while loop
      
      # Save proportions
      p = sapply(c('A1', 'A2', 'A3', 'AC',
                   'B1', 'B2', 'B3', 'BC'), function(ing){
                     pockets = sapply(agents_pockets, function(pocket) ing %in% names(pocket))
                     sum(as.numeric(pockets))/length(agents_pockets)
                   })
      props = rbind(props,
                    data.frame(network_type = network_type,
                               degree = degree,
                               pop_size = pop_size,
                               iteration = i,
                               epoch = t,
                               timestep = tt,
                               medicin = c('A1', 'A2', 'A3', 'AC',
                                           'B1', 'B2', 'B3', 'BC'),
                               proportion = p))
      
    } # End t loop
    
    # Return timings dataframe
    return(list(timings, props))
    
  }, mc.cores = n_cores) # End iteration loop
  
  # Return all values
  timings_binded = lapply(out, function(o) o[[1]]) %>% bind_rows()
  props_binded = lapply(out, function(o) o[[2]]) %>% bind_rows()
  out_binded = list(timings_binded, props_binded)
  return(out_binded)
  
} # End simCE

# Apply function to each network and bind the two dataframes
all_out = lapply(1:length(list_networks), function(x){
  if(length(list_networks) == 0) stop('List networks empty.You probably made a mistake in the path to the network folder.')
  simCE(list_networks[[x]], names(list_networks)[x], n_iter)
})
timings_all = lapply(all_out, function(o) o[[1]]) %>% bind_rows()
props_all = lapply(all_out, function(o) o[[2]]) %>% bind_rows()

# Get time to cross-over and others
timings_all$network_type = as.character(timings_all$network_type)
dat = timings_all %>% 
  group_by(network_type, degree, pop_size, iteration, innovation_score) %>%
  summarise(epoch = min(epoch))
innovations = c('A1/B1', 'A2/B2', 'A3/B3', 'AC/BC')
scores = c(48, 109, 188, 358)
dat$innovation = innovations[sapply(dat$innovation_score, function(x) which(scores == x))]
dat$combined = paste(dat$pop_size, dat$network_type, dat$degree, sep = '_')
summary_dat = dat %>% group_by(combined, innovation) %>% summarise(med_epoch = median(log(epoch)))
score_to_crossover = summary_dat[summary_dat$innovation == 'AC/BC',]
score_to_crossover = score_to_crossover[order(score_to_crossover$med_epoch, decreasing = T),]
dat$combined = factor(dat$combined, c(score_to_crossover$combined))

# Summarise diversity
props_all$combined = paste(props_all$pop_size, props_all$network_type, props_all$degree, sep = '_')
div = props_all %>% 
  group_by(combined, timestep, medicin) %>% 
  summarise(proportion_mean = mean(proportion))
div$lineage = ifelse(div$medicin %>% str_detect('A'), 'A', 'B')
div$progress = ifelse(div$medicin %>% str_detect('2'), 2, 1)
div$progress = ifelse(div$medicin %>% str_detect('3'), 3, div$progress)
div$progress = ifelse(div$medicin %>% str_detect('C'), 4, div$progress)

# Save the data
timings_all = crossdat
save(timings_all,
     file = paste0(path_out,'results_ABM1_', str_replace(Sys.time(), ' ', '_'), '.RData'))
save(div,
     file = paste0(path_out,'results_ABM1_div_', str_replace(Sys.time(), ' ', '_'), '.RData'))

# Save data for diversity plots ----
# !!!Only save data for few runs!!!, otherwise the file will be very large.
# This data is needed to replicate the diversity plots.

save(props_all, file = paste0(path_out, 'props_all.RData'))
