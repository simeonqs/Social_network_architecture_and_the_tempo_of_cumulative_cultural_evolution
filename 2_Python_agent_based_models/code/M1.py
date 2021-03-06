#!/usr/bin/env python
# coding: utf-8

import numpy as np
import networkx as nx
from networkx.algorithms import community
import random
import warnings
warnings.filterwarnings('ignore')
import os
from os import listdir
from os.path import isfile, join
from sys import argv

#this parses the command line argument, which sets the graph architecture to be simulated. Can take string values: clustered, degree, full, modular_and_clustered, modular, multilevel, small_world
graph_type = str(argv[1]) if len(argv) > 1 else 0

#parameters
masterID = 0 # an integer used to keep track of individual agents
master_sim_no = 0 # an integer used to keep track of individual simulations
#unused parameters to test pop turnover
turnover = False
num_turnover = 10

#loads graphs of the network type selected by cmd line argument
def get_graphs(graph_type):
    graph_files = [join("../networks",graph_type,file) for file in listdir("../networks/{}".format(graph_type)) if not file.startswith('.')]
    return graph_files

#function that generates the underlying network from an adjacency list and populates it with instances of the agent class
def generate_network(graph_file):
    G = nx.read_adjlist(graph_file,nodetype=int)
    for x in list(G.nodes()):
        G.nodes[x]['data'] = agent()
    return G

#this is the agent class, an instance of which is generated for each node in the graph at the beginning of simulations.
#it keeps track of the agents knowledge state, as well as contains methods for learning and producing behaviors
class agent:

    '''
    instances of agents are created and attached to nodes in the network
    agent in node x can be referenced with G.nodes[x]["data"]
    '''
    def __init__(self):
        global masterID
        #this variable keeps agents identifiable across simulations
        self.id = masterID
        masterID += 1

        '''
        this is the inventory of the agent
        each element has 3 entries: name, medicinal value, knowledge (0:naive,1:knowledgable)
        x[:,0] #output list of ingredient names
        x[:,1] #output list of ingredient values
        x[:,2] #output list of ingredient knowledge
        labels: i1a,i2a,i3a,A1,A2,A3,CA
        payoffs = [6,8,10,48,109,188,358]
        '''
        self.inventory = np.array([
              [1,6,1], #ia1
              [2,8,1], #ia2
              [3,10,1], #ia3
              [4,48,0], #A1
              [5,109,0], #A2
              [6,188,0], #A3
              [7,358,0], #CrossoverA
              [8,6,1], #ib1
              [9,8,1], #ib2
              [10,10,1], #ib3
              [11,48,0], #B1
              [12,109,0], #B2
              [13,188,0], #B3
              [14,358,0] #CrossoverB
             ])

    '''
    calculates the probability weights for choosing different items in inventory
    returns list of items
    num_items indicates whether they return 1 or 2 items, determined at chance during interaction
    '''
    def produce_items(self,num_items):
        #generate a list of medicinal values for items which are in inventory
        current_inventory = self.inventory[:,1] * self.inventory[:,2]
        #generate a list of weights that sum to 1 for use when choosing production
        probs = np.divide(current_inventory, np.sum(current_inventory))
        #print(probs)
        #randomly choose items without replacement from inventory
        production = np.random.choice(self.inventory[:,0],size=num_items, replace=False, p=probs)

        return list(production)



def choose_dyad(G,focal):
    '''
    chooses random agent and partner
    chooses parter in proportion to edgeweights if present
    '''
    agent_i = focal
    neighbors = [node for node in G.neighbors(agent_i)]
    if nx.is_weighted(G):
        weights = [G[agent_i][node]['weight'] for node in neighbors]
        probs = np.divide(weights, sum(weights))
        agent_j = np.random.choice(neighbors, p = probs)
    else:
        #randomly selects connected neighbor of agent_i
        agent_j = np.random.choice(neighbors)

    dyad = [agent_i,agent_j]

    return dyad

def discovery_check(combination):
    '''
    checks whether the combination of inventory items is a valid triad
    returns associated discovery with valid combination
    '''
    combination_set = set(combination)
    possible_discoveries = (4,11,5,12,6,13,7,14)
    names = ("A1","B1","A2","B2","A3","B3","XA","XB")
    successful_triads = [(1,2,3),(8,9,10),(4,1,9),(11,2,3),(5,9,10),(12,8,2),(6,13,5),(13,6,12)]

    if 7 in combination:
        discovery = 7
        name = "XA"

    elif 14 in combination:
        discovery = 14
        name = "XB"

    else:
        for index, triad in enumerate(successful_triads):
            if combination_set == set(triad):
                discovery = possible_discoveries[index]
                name = names[index]
                break
            else:
                discovery = False
                name = "not valid"
    return discovery,name

def interaction(G, dyad):
    '''
    dyad interacts, each agents produces 1 and 2 items each, with order randomized
    returns list of 3 chosen items
    '''
    num_items = [1,2]
    np.random.shuffle(num_items)
    #produce items without replacement from both agents
    items_i = G.nodes[dyad[0]]["data"].produce_items(num_items[0])
    items_j = G.nodes[dyad[1]]["data"].produce_items(num_items[1])
    #concatenate items
    combination = items_i + items_j
    return combination

def ind_learning(G,dyad,discovery):
    '''
    updates personal information of interaction dyad by flipping memory from 0 to 1
    '''
    for agent in dyad:
        #print("Memory before learning {}".format(G.nodes[agent]["data"].inventory[discovery-1,2]))
        G.nodes[agent]["data"].inventory[discovery-1,2] = 1
        #print("Memory after learning {}".format(G.nodes[agent]["data"].inventory[discovery-1,2]))

def diffusion(G, dyad, discovery):
    '''
    diffuses innovations to all other agents connected to dyad
    this occurs any time a product is successfully produced, but need to check with original
    '''
    neighbors_i = list(G.neighbors(dyad[0]))
    #print(neighbors_i)
    neighbors_j = list(G.neighbors(dyad[1]))
    #print(neighbors_j)
    neighbors = neighbors_i + list(set(neighbors_j) - set(neighbors_i))
    #print(neighbors)
    for agent in neighbors:
        G.nodes[agent]["data"].inventory[discovery-1,2] = 1
        #print(G.nodes[agent]["data"].inventory[discovery-1,2])
    return neighbors

def turnover_event(G):
    '''
    bonus turnover function, replaces num_turnover agents randomly with naive agents
    '''
    turnover_list = np.random.choice(exposure_list,replace = False, size = num_turnover)
    for agent in turnover_list:
        G.nodes[agent]["data"] = agent()

def create_csv():
    #writes header for main data
    with open("../results/m1/TTC_"+str(graph_type)+".csv","w") as f:
        f.write("sim,timestep,epoch,graph_type,pop_size,agent_i,agent_j,discovery,innov_level,a_track,b_track\n")

def create_csv_proportions():
    #writes header for proportions data
    with open("../results/m1/proportions_"+str(graph_type)+".csv","w") as f:
        f.write("sim,epoch,graph_type,pop_size,count_A1,count_A2,count_A3,count_XA,count_B1,count_B2,count_B3,count_XB\n")

def write_csv(sim_num,timestep,epoch,graph_condition,pop_size,agent_i,agent_j,discovery,innov_level,a_track,b_track):
    #writes row of data
    with open("../results/m1/TTC_"+str(graph_type)+".csv","a") as f:
        f.write("{},{},{},{},{},{},{},{},{},{},{}\n".format(sim_num,timestep,epoch,graph_condition,pop_size,agent_i,agent_j,discovery,innov_level,a_track,b_track))

def count_pockets(G,sim_num,epoch,graph_condition,pop_size):
    count_A1 = 0
    count_A2 = 0
    count_A3 = 0
    count_XA = 0
    count_B1 = 0
    count_B2 = 0
    count_B3 = 0
    count_XB = 0
    crossover_spread = False

    for n in G.nodes():
        count_A1 += G.nodes[n]["data"].inventory[3,2]
        count_A2 += G.nodes[n]["data"].inventory[4,2]
        count_A3 += G.nodes[n]["data"].inventory[5,2]
        count_XA += G.nodes[n]["data"].inventory[6,2]
        count_B1 += G.nodes[n]["data"].inventory[10,2]
        count_B2 += G.nodes[n]["data"].inventory[11,2]
        count_B3 += G.nodes[n]["data"].inventory[12,2]
        count_XB += G.nodes[n]["data"].inventory[13,2]

        if count_XA >= (pop_size/2)+1 or count_XB >= (pop_size/2)+1:
            crossover_spread = True
    with open("../results/m1/proportions_"+str(graph_type)+".csv","a") as f:
        f.write("{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(sim_num,epoch,graph_condition,pop_size,count_A1,count_A2,count_A3,count_XA,count_B1,count_B2,count_B3,count_XB,crossover_spread))


def simulation(num_replicates):
    global master_sim_no
    graph_files = get_graphs(graph_type)
    print(graph_files)
    for graph_file in graph_files:
        #grab graph filename as graph_condition
        graph_condition = os.path.basename(graph_file).rpartition('.')[0]
        for sim_num in range(num_replicates):
            print("simulation {}".format(sim_num))
            G = generate_network(graph_file)
            N = len(G.nodes())
            epoch = 0
            timestep = 0
            end = 0
            innovations = []
            a_track = 0
            b_track = 0
            #loops through interactions and exits once the final product is innovated
            while not end: #force sims to continue running
                #generate list of agents and randomly shuffle them
                focal_list = np.arange(1,N+1)
                np.random.shuffle(focal_list)

                #loops through agents, now in random order
                for focal in focal_list:
                    #picks a neighboring agent of focal agent
                    dyad = choose_dyad(G,focal)
                    #print("Agents {} and {} are interacting".format(dyad[0],dyad[1]))

                    #they produce some combination of items
                    combination = interaction(G, dyad)
                    #print("Agents produced combination {}".format(combination))

                    #check to see if there has been a discovery
                    discovery, discovery_name = discovery_check(combination)

                    #if so, and if it is not an innovation, add it to innovations
                    if discovery:
                        if discovery_name not in innovations:
                            innovations.append(discovery_name)
                            if discovery <= 7:
                                a_track += 1
                            else:
                                b_track += 1
                            #print("Combination led to discovery {} ({})".format(discovery_name,discovery))
                            write_csv(sim_num, timestep,epoch,graph_condition,N,dyad[0],dyad[1],discovery,len(innovations),a_track,b_track)

                        #both agents in the dyad learn the discovery
                        ind_learning(G,dyad,discovery)

                        #diffuse to all neighbors of the dyad (one to many "broadcast" diffusion)
                        neighbors = diffusion(G,dyad,discovery)
                        #print("information diffused to {} agents".format(len(neighbors)))

                        #If it is an innovation, and it is the end product, end the simulation
                        if discovery in [7,14]:
                            print("Crossover found. Simulation ending at timestep {}, epoch {}.".format(timestep, epoch))
                            end = 1
                            break
                        #write data here

                    #if the timestep is divisible by 100, write out the data
                    elif timestep%100==0 and not end:
                        write_csv(sim_num, timestep,epoch,graph_condition,N,dyad[0],dyad[1],discovery,len(innovations),a_track,b_track)

                    timestep +=1 #increment timestep

                #count_pockets(G,sim_num,epoch,graph_condition,N)
                epoch+=1 #increment epoch
            else:
                #code here will run when while loop condition evaluates False
                master_sim_no += 1

create_csv()
create_csv_proportions()
simulation(5000)
