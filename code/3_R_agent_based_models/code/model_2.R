# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Project: Cantor et al. 2020
# Date started: 04-07-2020
# Date last modified: 28-09-2020
# Author: Simeon Smeele
# Description: This is the R version of the second agent based model in Cantor et al. 2020. It loads the 
# networks and saves the results from the paths specified in the DATA section. It will install and load 
# required packages. Make sure to update the paths to fit your folder location. 
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# DATA ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Loading libraries
libraries = c('data.table', 'tidyverse', 'parallel')
for(i in libraries){
  if(! i %in% installed.packages()) lapply(i, install.packages)
  lapply(libraries, require, character.only = TRUE)
}

# Clean R
rm(list=ls()) 
dev.off()
cat("\014")  

# Paths
path_networks = 'Social_network_architecture_and_the_rate_of_cumulative_cultural_evolution/code/3_R_agent_based_models/networks'
path_out = 'Social_network_architecture_and_the_rate_of_cumulative_cultural_evolution/code/3_R_agent_based_models/output/model 2/'

# Settings
n_iter = 10000 # number of iterations
n_cores = 40 # number of cores to use (needs to be 1 on Windows)

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

# Load latest results if present
results = list.files(path_out, '*RData', full.names = T)
if(length(results) > 0) load(results[length(results)])

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

# Plants, lower case are plants, upper case are medicin
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
  
  # Make sure agents names are numerics
  rownames(network) = colnames(network) = 1:nrow(network)
  
  # Include edge to self in network
  diag(network) = 1
  
  # Run through iterations
  timings = mclapply(1:n_iter, function(i){
    
    # Equip agents
    agents_pockets = lapply(1:nrow(network), function(x) plants)
    
    # Get basics
    network_type = str_split(name_network, '_')[[1]][1]
    degree = str_split(name_network, '_')[[1]][2] %>% str_remove('d') %>% as.integer
    pop_size = nrow(network)
    
    # Running through epochs
    cont = T
    cross_known = F
    t = 0 # epochs
    tt = 0 # time-steps/interactions
    while(cont){ # run until 50% has crossover
      
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
        
        # Agents comparing notes, is a bit of a workaround, first find the unique names, then take the items
        # from plants and combinations, then remove the ones where e.g. a plant name was not found in the combination
        # list
        all_known = unique(c(names(agents_pockets[[agent_i]]), names(agents_pockets[[agent_j]])))
        new_pocket = append(plants[all_known], combinations[all_known]) 
        new_pocket = new_pocket[!is.na(names(new_pocket))]
        agents_pockets[[agent_i]] = new_pocket
        agents_pockets[[agent_j]] = new_pocket
        
        # Pick ingredients
        pi = as.numeric(agents_pockets[[agent_i]])/sum(as.numeric(agents_pockets[[agent_i]]))
        pj = as.numeric(agents_pockets[[agent_j]])/sum(as.numeric(agents_pockets[[agent_j]]))
        how_many = sample(2) # randomise which agent picks one vs two items
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
        
        # Check if the crossover was an ingredient
        if('AC' %in% ing_choosen | 'BC' %in% ing_choosen){
          choose_crossover = ing_choosen[ing_choosen %in% c('AC', 'BC')][1]
          invention = list(358)
          names(invention) = choose_crossover
        }
        
        # Save if invention, make sure they didn't already have it
        if(! is.na(invention)){
          if(! names(invention) %in% names(new_pocket)){
            agents_pockets[[agent_i]] = append(agents_pockets[[agent_i]], invention)
            agents_pockets[[agent_j]] = append(agents_pockets[[agent_j]], invention)
            
            # Save if the first time crossover
            if(!cross_known & invention == 358){
              timings = data.frame(network_type = network_type,
                                   degree = degree,
                                   pop_size = pop_size,
                                   iteration = i,
                                   epoch = t,
                                   first_or_50 = 'first')
              cross_known = T
            }
            
          }
        } # End if invention
        
      } # End while loop for agents
      
      # Check if 50% found crossover
      found_crossover = sapply(agents_pockets, function(x) max(unlist(x)) == 358)
      fraction_crossover = sum(as.numeric(found_crossover))/length(found_crossover)
      if(fraction_crossover > 0.5){
        # Save innovation level
        timings = rbind(timings, 
                        data.frame(network_type = network_type,
                                   degree = degree,
                                   pop_size = pop_size,
                                   iteration = i,
                                   epoch = t,
                                   first_or_50 = '50'))
        # Stop while loop (t)
        cont = F
      } 
      
    } # End t loop
    
    # Return timings dataframe
    return(timings)
    
  }, mc.cores = n_cores) # End iteration loop
  
  # Return all values
  timings = bind_rows(timings)
  return(timings)
  
} # End simCE

# Apply function to each network and bind the two dataframes
timings_all = lapply(1:length(list_networks), function(x) 
  simCE(list_networks[[x]], names(list_networks)[x], n_iter)) %>% bind_rows

# Save the data
save(timings_all,
     file = paste0(path_out,'results_ABM2_', str_replace(Sys.time(), ' ', '_'), '.RData'))

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# PLOTTING ----
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Combine network_type, degree and pop_size, do colour and order popsize
timings_all$combined = paste(timings_all$pop_size, timings_all$network_type, timings_all$degree, sep = '_')
colours = c('#922B21', '#633974', '#21618C', '#0E6655', '#9A7D0A', '#1C2833', '#EC407A')
timings_all$pop_size = factor(timings_all$pop_size, c('324', '144', '64'))

# Plot time to crossover
crossdat = timings_all[timings_all$first_or_50 == 'first',]
summary_dat = crossdat %>% group_by(combined) %>% summarise(med_epoch = median(log(epoch)))
summary_dat = summary_dat[order(summary_dat$med_epoch, decreasing = T),]
crossdat$combined = factor(crossdat$combined, c(summary_dat$combined))
pdf(paste0(path_out, 'time_to_crossover.pdf'), 10, 20)
ggplot(crossdat, aes(combined, epoch)) + 
  geom_jitter(height = 0, colour = alpha(colours[crossdat$network_type %>% as.factor %>% as.integer], 0.4)) + 
  geom_boxplot(outlier.colour = alpha('white', 0.0), 
               fill = alpha('white', 0.7)) +
  scale_y_continuous(trans = 'log2') +
  coord_flip() + 
  facet_grid(pop_size ~., scales = 'free') +
  labs(y = 'time to crossover', x = 'network type') +
  theme_light()
dev.off()

# Plot time to 50% knows crossover
dat_50 = timings_all[timings_all$first_or_50 == '50',]
summary_dat = dat_50 %>% group_by(combined) %>% summarise(med_epoch = median(log(epoch)))
summary_dat = summary_dat[order(summary_dat$med_epoch, decreasing = T),]
dat_50$combined = factor(dat_50$combined, c(summary_dat$combined))
pdf(paste0(path_out, 'time_to_50.pdf'), 10, 20)
ggplot(dat_50, aes(combined, epoch)) + 
  geom_jitter(height = 0, colour = alpha(colours[dat_50$network_type %>% as.factor %>% as.integer], 0.4)) + 
  geom_boxplot(outlier.colour = alpha('white', 0.0), 
               fill = alpha('white', 0.7)) +
  scale_y_continuous(trans = 'log2') +
  coord_flip() + 
  facet_grid(pop_size ~., scales = 'free') +
  labs(y = 'time to 50% knows crossover', x = 'network type') +
  theme_light()
dev.off()






