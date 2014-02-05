#!/usr/bin/python


"""
flu_simulator.py is a simulator for a flu-like infection spreading across a 
network between airports (nodeS) via air travel routes (edges). The goal of this
simulation is to test a genetic algorithim to find the optimal vaccination 
strategy for the given network. Data is loaded with command-line arguments such
as:

    flu_simulator.py <airport database> <route database>

"""

# Title:  air_travel_fl.py
# Author: Nicholas A. Yager
# Date:   2013-01-12

import copy
import getopt
import math
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from multiprocessing import Process
import networkx as nx
import os
import random
import sys 
import time


def main(argv):
    """
    Primary function that initiates network creation and handles execution of
    infection simulations.

    Args:
        argv: A list of command-line arguments passed to the application.

    Returns:
        Void

    """
    
    target = 3682

    random.seed(100)
    chromosome_ids = 1

    # Create the network using the command arguments.
    network = create_network(argv[0], argv[1])

    # Make a directory for the data, and change into that directory.
    currenttime = time.strftime("%Y-%m-%dT%H%M%S", time.gmtime())
    os.makedirs(currenttime)
    os.chdir(currenttime)

    # Simulate the infection for a baseline.
    baseline = infection(network, None, target)

    # Create 50 vaccination strategies.
    n_strat = 50
    max_generations = 15
    nodes = network.nodes()
    vaccinations = list()
    for i in range(0,n_strat):
        number_of_airports = random.randint(1,3000)
        airports = random.sample(nodes,number_of_airports)
        vaccinations.append(Vaccination(airports, chromosome_ids))
        chromosome_ids += 1

    
        
    for generation in range(0,max_generations):

        print(generation)

        # Test the strategies

        test_vaccinations(network, vaccinations, target)
        calculate_fitnesses(vaccinations)

        vaccinations = copy.deepcopy(sorted(vaccinations, key=lambda k: k.shared_fitness,reverse=True) )

        infected = list()
        closures = list()
        colors = list()
        highest_fitness = vaccinations[0].shared_fitness
        lowest_fitness = vaccinations[-1].shared_fitness
        print("High: ",highest_fitness," Low: ",lowest_fitness)
        
        if generation is max_generations - 1:
            print(vaccinations[0].airports)
            continue

        for strat in vaccinations:
            infected.append(strat.infected)
            closures.append(int(strat.closures))
            relative_fitness = (strat.shared_fitness - lowest_fitness) / (highest_fitness - lowest_fitness)
            colors.append(str(1-relative_fitness))
            plt.text(strat.closures,strat.infected,strat.ID)

        for index in range(int(n_strat/2), n_strat):
            parent = copy.deepcopy(random.choice(vaccinations[0:int(n_strat/5)]))
            vaccinations[index] = Vaccination(parent.airports,
                                              chromosome_ids)
            option = random.randint(1,10)
            if option == 1:
                vaccinations[index].mutate(network.nodes())
            elif option == 2:
                other_strat = copy.deepcopy(random.choice(vaccinations[0:int(n_strat/5)]))
                vaccinations[index].recombine(other_strat)
            else:
                print(vaccinations[index].ID, "- No genetic change")
            chromosome_ids += 1


    plt.axis([0,3300,0,3300])
    plt.scatter(closures, infected, color=colors, cmap=plt.cm.spectral)
    plt.show()

class Vaccination:

    def __init__(self, airports, ID):
        self.closures = len(airports)
        self.airports = copy.deepcopy(airports)
        self.fitness = 0
        self.shared_fitness= 0
        self.ID = ID

    def cleanDuplicates(self):

        self.airports = list(set(self.airports))
        self.closures = len(self.airports)

    def recombine(self, other_strategy):
        print(self.ID, "- Recombining with another strategy.")
        number_to_take = random.randint(1, other_strategy.closures)
        self.airports.extend(copy.deepcopy(random.sample(other_strategy.airports, number_to_take)))
        self.closures = len(self.airports)

        self.cleanDuplicates()

    def mutate(self, nodes):
        operation = random.randint(1,2)

        if self.closures < 3:
            operation = 1

        if operation is 1:
            # Add a new airport.
            print(self.ID, "- Adding an airport")
            self.airports.append(random.choice(nodes))
            self.closures += 1

        if operation is 2:
            print(self.ID, "- Removing an airport")
            num_to_remove = random.randint(1,int(self.closures/2))
            to_remove = sorted(random.sample(range(0,self.closures-1), 
                                             num_to_remove),
                               reverse=True)
            for item in to_remove:
                self.airports.pop(item)
            self.closures -= 1
        self.cleanDuplicates()

def test_vaccinations(network, vaccinations, target):

    for vaccination in vaccinations:
        print(vaccination.ID)
        results = infection(network, vaccination.airports, target, False)
        d = math.sqrt( math.pow(results["Suscceptable"],2) + math.pow(vaccination.closures,2))
        vaccination.fitness = (results["Suscceptable"] - vaccination.closures) +\
                              d
        vaccination.infected = results["Infected"] + results["Recovered"]


def calculate_fitnesses(vaccinations):
    """
    Calculate the shared Pareto Niche fitness for a list of vaccination 
    strategies.

    Args:
        vaccinations: A list of vaccination strategies.
    """
    # Calculate the fitness with a pareto niche algorithim
    for i in range(0,len(vaccinations)-1):
        vaccination = vaccinations[i]
        infected = vaccination.infected
        closures = vaccination.closures

        Mi = 0
        # Compare with the other vaccinations
        for j in range(0,len(vaccinations)-1):
            # Dont compare yourself.
            if j is i:
                continue

            other_vaccination = vaccinations[j]
            other_infected = other_vaccination.infected
            other_closures = other_vaccination.closures

            # Distance calculation
            delta_infected = math.pow(other_infected - infected, 2)
            delta_closures = math.pow(other_closures - closures, 2)
            distance = math.sqrt( delta_infected + delta_closures)

            # Shared function
            Sigma_share = 500
            if distance <= Sigma_share:
                Mi += 1 - (distance/Sigma_share)
            elif distance == Sigma_share:
                Mi += 1

        # Shared fitness:
        if Mi > 0:
            vaccination.shared_fitness = vaccination.fitness / Mi    
        else:
            vaccination.shared_fitness = vaccination.fitness

    return vaccinations

def create_network(nodes, edges):
    """
    Create a NetworkX graph object using the airport and route databases.

    Args:
        nodes: The file path to the nodes .csv file.
        edeges: The file path to the edges .csv file.

    Returns:
        G: A NetworkX DiGraph object populated with the nodes and edges assigned
           by the data files from the arguments.
           
    """

    print("Creating network.")
    G = nx.DiGraph()

    print("\tLoading airports", end="")
    sys.stdout.flush()
    # Populate the graph with nodes.
    with open(nodes) as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")

           # if entries[3] == "United States":
           #     G.add_node(entries[0], 
           #                name=entries[1], 
           #                lat=entries[6],
           #                lon=entries[7])
            G.add_node(int(entries[0]), 
                       name=entries[1], 
                       lat=entries[6],
                       lon=entries[7])


    print("\t\t\t\t\t[Done]")
    
    print("\tLoading routes",end="")
    sys.stdout.flush()
    # Populate the graph with edges.
    edge_count = 0
    error_count = 0
    duplicate_count = 0
    line_num = 1
    with open(edges) as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")
            try:
                if G.has_edge(int(entries[3]),int(entries[5])):
                    duplicate_count += 1
                else:
                    if line_num > 1:
                        G.add_edge(int(entries[3]), int(entries[5]) )
                        edge_count += 1
            except ValueError:
                # The value doesn't exist
                error_count += 1
                pass
            line_num += 1
    
    print("\t\t\t\t\t\t[Done]")


    
    # Remove nodes without inbound edges
    indeg = G.in_degree()
    outdeg = G.out_degree()
    to_remove = [n for n in indeg if (indeg[n] + outdeg[n] < 1)]
    
    G.remove_nodes_from(to_remove)

    return G

def infection(input_network, vaccination, start, visualize = False):
    """
    Simulate an infection within network, generated using seed, and with the
    givin vaccination strategy. This function will write data from each timestep
    to "data.csv".

    Args:
        network: A NetworkX DiGraph object.
        vaccination: A list of node indices to label as recovered from the 
                     begining.

    Returns:
        state: A dictionary of the total suscceptable, infected, and recovered.

    """

    print("Simulating infection.")

    network = input_network.copy()


    # Open the data file
    f = open("data.csv", "w")
    f.write("time, s, i, r\n")

    # Set the default to susceptable
    print("\tInitializing network for infection.", end="")
    sys.stdout.flush()
    for node in network.nodes():
        network.node[node]["status"] =  "s"
        network.node[node]["color"] = "b"
        network.node[node]["age"] = 0
    
    # Add in the recovering
    if vaccination is not None:
        vaccination = sorted(vaccination)
        for node in vaccination:
            if node != start:
                # Get the node's predecessors and sucessors
                remove_predecessors = [ (node, suc) for suc in network.predecessors(node)]
                remove_successors = [ (node, suc) for suc in network.successors(node)]
                network.remove_edges_from(remove_predecessors)
                network.remove_edges_from(remove_successors)

    # Assign the infected
    infected = start
    network.node[infected]["status"] = "i"
    network.node[infected]["color"]  = "orange"

    print("\t\t\t[Done]")

    print("\tInitial vector: "+network.node[infected]["name"])
    if vaccination is not None:
        print("\tVaccinated: ", len(vaccination) )
    else: 
        print("\tVaccinated: None")

    if visualize:
        # Calculate the layout of the network to ensure even plotting.
        pos = nx.spring_layout(network)

    # Iterate through the evolution of the disease.
    for step in range(0,99):

        # Create variables to hold the outcomes as they happen
        S, I, R = 0, 0, 0

        for node in network.nodes():
            status = network.node[node]["status"]
            age = network.node[node]["age"]
            color = network.node[node]["color"]

            if status is "i" and age >= 5:
                # The infected has reached its recovery time
                network.node[node]["status"] = "r"
                network.node[node]["color"] = "g"

            elif status is "i":
                # Propogate the infection.
                if age > 0:
                    possible_victims = network.successors(node)

                    #n = int(len(possible_victims)/)
                    
                    #n = 3

                    #if n > len(possible_victims):
                    #    n = len(possible_victims) 

                    #victims = random.sample(possible_victims, n)
                    victims = possible_victims

                    for infectees in victims:
                        infect_status = network.node[infectees]["status"]
                        if infect_status != "r" and infect_status != "i":
                            network.node[infectees]["status"] = "i"
                            network.node[infectees]["age"] = 0
                            network.node[infectees]["color"] = "orange"
                network.node[node]["age"] += 1

        # Loop twice to prevent bias.
        for node in network.nodes():
            status = network.node[node]["status"]
            age = network.node[node]["age"]
            color = network.node[node]["color"]

            if status is "s":
                # Count those susceptable
                S += 1

            elif status is "r":

                R += 1

            elif status is "i":
                
                I += 1

        printline = "{0}, {1}, {2}, {3}".format(step, S, I, R)
        f.write(printline + "\n")

       # print("\t"+printline)

        if I is 0:
            break

        if visualize:
            visualize(network, pos)
        
    print("\t----------\n\tS: {0}, I: {1}, R: {2}".format(S,I,R))

    return {"Suscceptable":S,"Infected":I, "Recovered":R}

       
def visualize(network, pos):
    """
    Visualize the network given an array of posisitons.
    """
    colors = []
    for node in network.nodes():
        colors.append(network.node[node]["color"])

    plt.figure(figsize=(8,8))

    nx.draw_networkx_nodes(network,
            pos,
            node_size=25,
            with_labels=False,
            node_color = colors)

    nx.draw_networkx_edges(network,pos,
                           width=1,
                           edge_color='black',
                           arrows=False)

    plt.axis('off')

    number_files = str(len(os.listdir()))
    while len(number_files) < 3:
        number_files = "0" + number_files

    
    plt.savefig("infection-{0}.png".format(number_files),
                bbox_inches='tight', dpi=100 
            )

if __name__ == "__main__":
       main(sys.argv[1:])
