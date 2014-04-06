#!/usr/bin/python3
"""
flu_simulator.py is a simulator for a flu-like infection spreading across a 
network between airports (nodes) via air travel routes (edges). The goal of this
simulation is to test a genetic algorithim to find the optimal vaccination 
strategy for the given network. Data is loaded with command-line arguments such
as:

    flu_simulator.py -gbdrus <airport database> <route database>

    Flags:
        -g: Run a genetic algorithm quarantine simulation.
        -b: Run a betweenness-based quarantine simulation.
        -d: Add a delay to quarantine simulations.
        -w: Run a simulation based on edge weights.
        -r: Run a random quarantine simulation.
	    -c: Run a simulation based on making clusters.
        -u: Convert the network to an undirected network.
        -s: Run a naive simulation and output the SIR data.
        -i: Filter to only quarantine international flights.
        -q: Filter to only quarantine domestic flights.
        -y: Disable recalculation of network edges.

"""

# Title:  flu_simulator.py
# Author: Nicholas A. Yager
# Date:   2013-01-12

import copy
import getopt
import math
import networkx as nx
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap as Basemap
import operator
import os
import random
import sys
from scipy import stats
import time

global VISUALIZE

def main():
    """
    Primary function that initiates network creation and handles execution of
    infection simulations.

    Args:
        argv: A list of command-line arguments passed to the application.

    Returns:
        Void

    """

    # Flag defaults
    GENETIC = False 
    BETWEENNESS = False 
    WEIGHT = False 
    RANDOM = False 
    SIR = False
    UNDIRECT = False
    VISUALIZE = False
    EDGES = False
    INTERNATIONAL = False
    DOMESTIC = False
    CLUSTER = False
    RECALCULATE = True
    DELAY = 0

    # Determine the parameters of the current simulation.
    opts, args = getopt.getopt(sys.argv[1:], "wbgd:iqcrusvey", [ "Weight",
                                                            "Betweenness",
                                                            "Genetic",
                                                            "Delay",
                                                            "Cluster",
                                                            "International",
                                                            "Domestic",
                                                            "Random",
                                                            "Undirect",
                                                            "SIR",
                                                            "visualize",
                                                            "Recalculate",
                                                            "ByEdge" ]
                                                            )

    # Check if the data arguments are available
    if len(args) < 2:
        print("\nflu_simulator.py -gbdrus <airport database> <route database>\n")
        print("Flags:\n\t-g: Run a genetic algorithm quarantine simulation.")
        print("\t-y: Disable network edge recalculation.")
        print("\t-b: Run a betweenness-based quarantine simulation.")
        print("\t-d: Run a degree-based quarantine simulation.")
        print("\t-r: Run a random quarantine simulation.")
        print("\t-u: Convert the network to an undirected network.")
        print("\t-s: Run a naive simulation and output the SIR data.")
        print("\t-c: Run a simulation based on local clustering coefficients.")
        print("\t-q: Run a simulation based on domestic quarantine.")
        print("\t-i: Run a simulation based on international quarantine.\n")
        exit()


    AIRPORT_DATA = args[0]
    ROUTE_DATA = args[1]



    for o, a in opts:
        if o == "-g":
            GENETIC = True
            BETWEENNESS = False
            DEGREE = False
        elif o == "-b":
            GENETIC = False
            BETWEENNESS = True
        elif o == "-w":
            GENETIC = False
            WEIGHT = True
        elif o == "-d":
            DELAY = a
        elif o == "-r":
            GENETIC = False
            RANDOM = True
        elif o == "-s":
            SIR = True
        elif o == "-u":
            UNDIRECT = True
        elif o == "-v":
            VISUALIZE = True
        elif o == "-e":
            EDGES = True
        elif o == "-i":
            INTERNATIONAL = True
        elif o == "-q":
            DOMESTIC = True
        elif o == "-c":
            CLUSTER = True
        elif o == "-y":
            RECALCULATE = False
            
    #NUM_SIMULATIONS = 344
    NUM_SIMULATIONS = 100

    seed = 100
    random.seed(seed)

    # Identify the script.
    print("Flu Simulator 0.10.6")
    print("Created by Nicholas A. Yager and Matthew Taylor\n\n")

    # Create the network using the command arguments.
    network = create_network(AIRPORT_DATA, ROUTE_DATA)
  
    degrees = network.degree() 
    weights = dict()
    for airport, degree in degrees.items():
        weights[airport] = network.out_degree(airport) +\
                           network.in_degree(airport)
    
    target = list()
    for ind in range(0,NUM_SIMULATIONS):
        target_round = list()
        while len(target_round) < 10:
             chosen_airport = weighted_random(weights)
             if chosen_airport not in target_round:
                 target_round.append(chosen_airport)
        target.append(target_round)

    if UNDIRECT:
        network = network.to_undirected()

    # Make a directory for the data, and change into that directory.
    currenttime = time.strftime("%Y-%m-%dT%H%M%S", time.gmtime())
    os.makedirs(currenttime)
    os.chdir(currenttime)

    # Record relevent data about the simulation.
    # simulation_data(network, currenttime, target, seed)

    if GENETIC:
        genetic_simulations(network, 50, 50, target)

    if BETWEENNESS:
        betweenness_simulations(network, target, VISUALIZE, EDGES, DELAY,
                                INTERNATIONAL, DOMESTIC, RECALCULATE)
    
    if WEIGHT:
        degree_simulations(network, target, VISUALIZE, EDGES, DELAY,
                           INTERNATIONAL, DOMESTIC, RECALCULATE)

    if RANDOM:
        random_simulations(network, target, VISUALIZE, EDGES, DELAY,
                           INTERNATIONAL, DOMESTIC, RECALCULATE)

    if SIR:
        sir_simulations(network, target, VISUALIZE, EDGES, DELAY, RECALCULATE)

    if CLUSTER:
        cluster_simulations(network, target, VISUALIZE, EDGES, DELAY,
                            INTERNATIONAL, DOMESTIC, RECALCULATE)

def weighted_random(weights):
    number = random.random() * sum(weights.values())
    for k,v in weights.items():
        if number <= v:
            break
        number -= v
    return k

def pad_string(integer, n):
    """
    Add "0" to the front of an interger so that the resulting string in n 
    characters long.

    Args:
        integer: The number to pad.
        n: The desired length of the string

    Returns
        string: The padded string representation of the integer.
        
    """

    string = str(integer)

    while len(string) < n:
        string = "0" + string

    return string

def sir_simulations(network, targets, VISUALIZE, EDGES, DELAY, RECALCULATE):
    """
    Run an infection simulation across the network for each of the given
    targets, and determine the median number of infected per day.

    Args:
        network: A NetworkX graph object.
        targets: A list of initial infection targets.

    Returns
        VOID.

    IO:
        sir.csv: The proportion of a population that is infected by time step.
    """

    print("SIR Mode")

    # Make a new folder for the degree data.
    os.makedirs("sir")

    iteration = 0
    
    for target in targets:
        print(target)
        sir_file = "sir/sir_{0}.csv".format(pad_string(iteration,4))

        results = infection(network, None, target, vis=VISUALIZE, file_name=sir_file )
        N = results["Suscceptable"] + results["Infected"] + results["Recovered"]
        iteration += 1


def simulation_data(network, time, targets, seed):
    """
    Output various statistics of the nature of the network to a file, including
    the diameter, the number of cycles, number of verticies and edges, and the
    average in and out degrees.

    Args:
        network: A NetworkX network graph.

    Returns:
        VOID

    IO:
        network.dat: A file with all of the relevent netowkr information.

    """

    print("\tWriting graph edgelist")

    print("\tCalculating edge betweenness centrality")
    edge_betweenness = nx.edge_betweenness_centrality(network)
    print(edge_betweenness)

    edgelist = open("edgelist.csv", "w")
    edgelist.write('"IATA_From","IATA_To","Edge_Betweenness"\n')

    for edge in network.edges():
        from_vertex = edge[0]
        to_vertex = edge[1]

        # Collect data
        IATAFrom = network[from_vertex][to_vertex]['IATAFrom']
        IATATo = network[from_vertex][to_vertex]['IATATo']
        betweenness = edge_betweenness[(from_vertex, to_vertex)]
        edgelist.write('"{0}","{1}",{2}\n'.format(IATAFrom, 
                                                  IATATo, 
                                                  betweenness))
        

    edgelist.close()

        

    print("\tDetermining network type.")
    # Determine if the graph is directed or undirected
    if isinstance(network,nx.DiGraph):
        network_type = "Directed"
    else:
        network_type = "Undirected"

    print("\tCalculaing edges and verticies.")
    # Number of verticies and edges
    edges = network.number_of_edges()
    verticies = network.number_of_nodes()

    
    # Not every vertex can lead to every other vertex.
    # Create a subgraph that can.
    print("\tTemporarily converting to undirected.")
    undirected = network.to_undirected()
    print("\tFinding subgraphs.")
    subgraphs = nx.connected_component_subgraphs(undirected)


    # Find the number of vertices in the diameter of the network

    print("\tFinding network diameter.")
    diameter = nx.diameter(subgraphs[0])

    # Find the number of simple cycles in the network.
    #print("\tCalcuating network cycles.")
    #simple_cycles = nx.simple_cycles(network)


    # Find the largest and smallest cycle in the network.
    smallest_cycle = 9999999999
    largest_cycle = 0
   
    number_cycles = 0

    #for cycle in simple_cycles:
    #    print(cycle)
    #    length = len(cycle)
    #    if length < smallest_cycle:
    #        smallest_cycle = length
    #    if length > largest_cycle:
    #        largest_cycle = length


    print("\tStoring network parameters")

    data_file = open("network.dat", "w")
    data_file.write("Simulation name: {0}\n\n".format(time))
    data_file.write("Network properties\n===============\n")
    data_file.write("Network type: {0}\n".format(network_type))
    data_file.write("Number of verticies: {0}\n".format(verticies))
    data_file.write("Number of edges: {0}\n".format(edges))
    data_file.write("Diameter: {0}\n".format(diameter))
    data_file.write("Number of cycles: {0}\n".format(number_cycles))
    data_file.write("Largest cycle: {0}\n".format(largest_cycle))
    data_file.write("Smallest cycle: {0}\n".format(smallest_cycle))

    data_file.close()

def random_simulations(network, targets, VISUALIZE, EDGES, DELAY, I, Q, RECALCULATE):
    """
    Simulate the spread of infection for increasing vaccination efforts by
    quarantining airports randomly.

    Args:
        network: A NetworkX graph object.

    Returns:
        VOID

    IO:
        random.csv: A gile with the number of total ingected people in the
                        network for each quarantine effort.
    """

    print("Random Mode")

    # Make a new folder for the degree data.
    os.makedirs("random")

    # Filter international flights if necessary
    randoms = random.sample(network.edges(), len(network.edges()))

    if I:
        for i,j in randoms:
            if network[i][j]["international"] == False:
                randoms.remove((i,j))
    if Q:
        for i,j in randoms:
            if network[i][j]["domestic"] == False:
                randoms.remove((i,j))

    iteration = 0
    for target in targets:

        # Open a file for this targ'ets dataset
        random_file = open("random/random_{0}.csv".format(pad_string(iteration,4)),"w")
        random_file.write('"effort","total_infected"\n')


        # Generate a baseline
        results = infection(network, None, target,DELAY=DELAY, vis=VISUALIZE, 
                            title="Random - 0%")
        total_infected = results["Infected"] + results["Recovered"]
        random_file.write("{0},{1}\n".format(0,total_infected))
        
        # Perform a check for every strategy
        for effort in range(1,101,5):
            max_index = int(len(randoms) * (effort/100))-1
            strategy = randoms[0:max_index]

            title = "random - {0}%".format(effort/100)
            results = infection(network, strategy, target, vis=VISUALIZE,
                                DELAY=DELAY,
                                title=title, inf_type=EDGES)
            total_infected = results["Infected"] + results["Recovered"]
            random_file.write("{0},{1}\n".format(effort/100,total_infected))
     
            if total_infected == 1:
                for remaining_effort in range(effort+5,101,5):
                    random_file.write("{0},{1}\n".format(remaining_effort/100,
                                                         total_infected))
                break
   
        iteration += 1
        random_file.close()



def degree_simulations(network, targets, VISUALIZE, EDGES, DELAY, I, Q, RECALCULATE):
    """
    Simulate the spread of infection for increasing vaccination efforts by 
    closing routes of decreasing degree-based weight.

    Args:
        network: A NetworkX graph object.
        targets: A list of the initial infection vertices.

    Returns:
        Void

    IO:
        degree.csv: A file with the number of total infected people in the
                         network for each vaccination effort.
    """
    print("Degree Mode.")
    print("\tCalculating degrees", end="")
    sys.stdout.flush()
    degrees = network.edges(data = True)

    index = 0
    if I:
        for i,j,data in degrees:
            if data["international"] == False:
                degrees.remove((i,j,data))
            index += 1
    if Q:
        for i,j,data in degrees:
            if data["domestic"] == False:
                degrees.remove((i,j,data))
            index += 1

    sorted_degree = sorted(degrees, key=lambda k: k[2]['weight'], reverse=True)
    degree = list()
    for degree_item in sorted_degree:
        degree.append((degree_item[0], degree_item[1]))
    print("\t\t\t\t\t[Done]")

    # Make a new folder for the degree data.
    os.makedirs("weight")

    iteration = 0
    for target in targets:

        # Open a file for this targ'ets dataset
        degree_file = open("weight/weight_{0}.csv".format(pad_string(iteration,4)),"w")
        degree_file.write('"effort","total_infected, edges_closed"\n')


        # Generate a baseline
        results = infection(network, None, target, vis=VISUALIZE,
                            DELAY=DELAY, title="weight - 0%")
        total_infected = results["Infected"] + results["Recovered"]
        degree_file.write("{0},{1}\n".format(0,total_infected))

        # Perform a check for every strategy
        for effort in range(1,101,5):
            max_index = int(len(degree) * (effort/100))-1
            strategy = degree[0:max_index]

            edges_closed = len(strategy)

            title = "weight - {0}%".format(effort/100)

            results = infection(network, strategy, target, DELAY=DELAY,
                                vis=VISUALIZE, 
                                title=title, inf_type=EDGES)
            total_infected = results["Infected"] + results["Recovered"]
            degree_file.write("{0},{1}\n".format(effort/100,
                                                 total_infected,
                                                 edges_closed))

            if total_infected == 1:
                for remaining_effort in range(effort+5,101,5):
                    degree_file.write("{0},{1}\n".format(remaining_effort/100,
                                                         total_infected))
                break

        
        iteration += 1
        degree_file.close()


def betweenness_simulations(network,targets, VISUALIZE, EDGES, DELAY, I, Q, RECALCULATE):
    """
    Simulate the spread of infection for increasing vaccination efforts by 
    quarantining airports of decreasing betweenness.

    Args:
        network: A NetworkX graph object.
        targets: A list of the initial infection vertices.

    Returns:
        Void

    IO:
        betweenness.csv: A file with the number of total infected people in the
                         network for each vaccination effort.
    """

    print("Betweenness Centrality Mode.")
    print("\tCalculating betweenness centrality", end="")
    sys.stdout.flush()

    betweennesses = nx.edge_betweenness_centrality(network,weight="weight")
    betweenness = sorted(betweennesses.keys(), 
                    key=lambda k: betweennesses[k], reverse=True)

    if I:
        for i,j in betweenness:
            if not network[i][j]["international"]:
                betweenness.remove((i,j))
    if Q:
        for i,j in betweenness:
            if not network[i][j]["domestic"]:
                betweenness.remove((i,j))

    print("\t\t\t\t[Done]")


    os.makedirs("betweenness")

    iteration = 0
    for target in targets:
    

        # Write the betweenness data to a folder.
        betweenness_file = open(
                            "betweenness/betweenness_{0}.csv".format( 
                                            pad_string(iteration,4)),
                            "w")
                           
        betweenness_file.write('"effort","total_infected"\n')

        # Generate a baseline
        results = infection(network, None, target, vis=VISUALIZE,
                            title="Betweenness - 0%",DELAY=DELAY)
        total_infected = results["Infected"] + results["Recovered"]
        betweenness_file.write("{0},{1}\n".format(0,total_infected))

        # Perform a check for every strategy
        for effort in range(1,101,5):
            max_index = int(len(betweenness) * (effort/100))-1
            strategy = betweenness[0:max_index]

            title = "betweenness - {0}%".format(effort/100)
            results = infection(network, strategy, target, vis=VISUALIZE,
                                title=title, inf_type=EDGES, DELAY=DELAY)
            total_infected = results["Infected"] + results["Recovered"]
            betweenness_file.write("{0},{1}\n".format(effort/100,total_infected))
            
            if total_infected == 1:
                for remaining_effort in range(effort+5,101,5):
                    betweenness_file.write("{0},{1}\n".format(remaining_effort/100,
                                                              total_infected))
                break

        iteration += 1
        betweenness_file.close()

def genetic_simulations(network, num_strategies, max_generations, target):
    """
    A procedure to stochastically simulate quarantine strategies using a
    niched Pareto genetic algorithm to determine the most fit strategies.

    Args:
        network: A NetworkX graph object.
        num_strategies: The number of strategies to test per simulation.
        max_generations: The number of generations to simulate.
        target: The initial infection vertex.
    
    Returns:
        Void.

    IO:
        pareto_*.csv: Files containing data on the fitness values for
            each strategy for each generation.
    """
    # Create 50 vaccination strategies.
    os.makedirs("pareto")
    nodes = network.nodes()
    vaccinations = list()
    chromosome_ids = 1

    for i in range(0,num_strategies):
        number_of_airports = random.randint(1,3000)
        airports = random.sample(nodes,number_of_airports)
        vaccinations.append(Vaccination(airports, chromosome_ids))
        chromosome_ids += 1

    
        
    for generation in range(0,max_generations):

        print(generation)

        # Test the strategies

        genetic_test(network, vaccinations, target)
        calculate_fitnesses(vaccinations)

        vaccinations = copy.deepcopy(sorted(vaccinations, 
                                            key=lambda k: k.shared_fitness,
                                            reverse=True) 
                                    )

        infected = list()
        closures = list()
        colors = list()
        highest_fitness = vaccinations[0].shared_fitness
        lowest_fitness = vaccinations[-1].shared_fitness
        print("High: ",highest_fitness," Low: ",lowest_fitness)
        

        # Open a file for the pareto data.
        pareto_file = open("pareto/pareto_{0}.csv".format(generation), 'w')
        pareto_file.write('"ID","infected","closures","fitness"\n')

        for strat in vaccinations:
            infected.append(strat.infected)
            closures.append(int(strat.closures))
            relative_fitness = (strat.shared_fitness - lowest_fitness) / \
                               (highest_fitness - lowest_fitness)
            pareto_file
            colors.append(str(1-relative_fitness))
            #plt.text(strat.closures,strat.infected,strat.ID)
            pareto_file.write('{0},{1},{2},{3}\n'.format(strat.ID,
                                                         infected[-1], 
                                                         closures[-1], 
                                                         strat.shared_fitness))


        pareto_file.close()

        if generation is max_generations - 1:
            print(vaccinations[0].airports)
            continue


        for index in range(int(num_strategies/2), num_strategies):
            limit = int(num_strategies/5)
            parent = copy.deepcopy(random.choice(vaccinations[0:limit]))
            vaccinations[index] = Vaccination(parent.airports,
                                              chromosome_ids)
            option = random.randint(1,2)
            if option == 1:
                vaccinations[index].mutate(network.nodes())
            elif option == 2:
                other_strat = copy.deepcopy(random.choice(vaccinations[0:limit]))
                vaccinations[index].recombine(other_strat)
            else:
                print(vaccinations[index].ID, "- No genetic change")
            chromosome_ids += 1
    return 


class Vaccination:

    def __init__(self, airports, ID):
        self.closures = len(airports)
        self.airports = copy.deepcopy(airports)
        self.fitness = 0
        self.shared_fitness= 0
        self.ID = ID

    def cleanDuplicates(self):

        #self.airports = list(set(self.airports))
        #self.closures = len(self.airports)
        pass

    def recombine(self, other_strategy):
        print(self.ID, "- Recombining with another strategy.")
        other_strategy = copy.deepcopy(other_strategy)
        number_to_take = random.randint(1, other_strategy.closures-1)
        print(number_to_take, len(other_strategy.airports))
        self.airports.extend(copy.deepcopy(random.sample(other_strategy.airports,
                                                         number_to_take)))
        self.closures = len(self.airports)

        self.cleanDuplicates()

    def mutate(self, nodes):
        operation = random.randint(1,2)

        if self.closures < 3:
            operation = 1

        if operation is 1:
            # Add a new airport.
            print(self.ID, "- Adding an airport")
            number_to_add = random.randint(1,int(self.closures/2))
            self.airports.extend(copy.deepcopy(random.sample(nodes,
                                                             number_to_add)))

        if operation is 2:
            print(self.ID, "- Removing an airport")
            num_to_remove = random.randint(1,int(self.closures/2))
            to_remove = sorted(random.sample(range(0,self.closures-1), 
                                             num_to_remove),
                               reverse=True)
            for item in to_remove:
                self.airports.pop(item)
        self.closures = len(self.airports)
        self.cleanDuplicates()

def cluster_simulations(network,targets, VISUALIZE, EDGES, DELAY, I, Q, RECALCULATE):
    """
    Simulate the spread of infection for increasing vaccination efforts by 
    quarantining airports by decreasing local clustering coefficients.

    Args:
        network: A NetworkX graph object.
        targets: A list of the initial infection vertices.

    Returns:
        Void

    IO:
        cluster.csv: A file with the number of total infected people in the
                         network for each vaccination effort.
    """

    print("Local Clustering Coefficient Mode.")
    sys.stdout.flush()
    
    clusters = network.edges(data = True)

    if I:
        for i,j,data in clusters:
            if not network[i][j]['international']:
                clusters.remove((i,j,data))
    if Q:
        for i,j,data in clusters:
            if not network[i][j]['domestic']:
                clusters.remove((i,j,data))
    
    sorted_cluster = sorted(clusters, key=lambda k: k[2]['cluster'], reverse=True)
    cluster = list()
    for cluster_item in sorted_cluster:
        cluster.append((cluster_item[0], cluster_item[1]))

    os.makedirs("cluster")

    iteration = 0
    for target in targets:
    

        # Write the cluster data to a folder.
        cluster_file = open(
                            "cluster/cluster_{0}.csv".format( 
                                            pad_string(iteration,4)),
                            "w")
                           
        cluster_file.write('"effort","total_infected"\n')

        # Generate a baseline
        results = infection(network, None, target, vis=VISUALIZE,
                            title="Cluster - 0%", DELAY=DELAY)
        total_infected = results["Infected"] + results["Recovered"]
        cluster_file.write("{0},{1}\n".format(0,total_infected))

        # Perform a check for every strategy
        for effort in range(1,101,5):
            max_index = int(len(cluster) * (effort/100))-1
            strategy = cluster[0:max_index]

            title = "cluster - {0}%".format(effort/100)
            results = infection(network, strategy, target, vis=VISUALIZE,
                                title=title, inf_type=EDGES, DELAY=DELAY)
            total_infected = results["Infected"] + results["Recovered"]
            cluster_file.write("{0},{1}\n".format(effort/100,total_infected))
            
            if total_infected == 1:
                for remaining_effort in range(effort+5,101,5):
                    cluster_file.write("{0},{1}\n".format(remaining_effort/100,
                                                              total_infected))
                break

        iteration += 1
        cluster_file.close()


def genetic_test(network, vaccinations, target):
    """
    genetic_test
    """
    for vaccination in vaccinations:
        print(vaccination.ID)
        results = infection(network, vaccination.airports, target, False)
        d = math.sqrt( math.pow(results["Infected"]+results["Recovered"],2) +\
                        math.pow(vaccination.closures,2))
        vaccination.fitness = -1 * (d)^2
        vaccination.infected = results["Infected"] + results["Recovered"]


def calculate_fitnesses(vaccinations):
    """
    Calculate the shared Pareto Niche fitness for a list of vaccination 
    strategies.

    Args:
        vaccinations: A list of vaccination strategies.
    """
    # Calculate the fitness with a pareto niche algorithim
    for strategy in vaccinations:

        Mi = 0
        # Compare with the other vaccinations
        for other_strategy in vaccinations:
            # Dont compare yourself.
            if strategy == other_strategy:
                continue

            # Distance calculation
            delta_infected = math.pow(other_strategy.infected -\
                                      strategy.infected,2)
            delta_closures = math.pow(other_strategy.closures -\
                                      strategy.closures, 2)
            distance = math.sqrt( delta_infected + delta_closures)

            # Shared function
            Sigma_share = 100
            if distance <= Sigma_share:
                Mi += 1 - (distance/Sigma_share)
            elif distance <= 1:
                Mi += 1

        # Shared fitness:
        if Mi > 0:
            strategy.shared_fitness = strategy.fitness / Mi    
        else:
            strategy.shared_fitness = strategy.fitness

    return
   
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
    with open(nodes, 'r', encoding='utf-8') as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")

           # if entries[3] == "United States":
           #     G.add_node(entries[0], 
           #                name=entries[1], 
           #                lat=entries[6],
           #                lon=entries[7])
            G.add_node(int(entries[0]),
                       country=entries[3],
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
    with open(edges, 'r', encoding="utf-8") as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")
            try:
                if G.has_edge(int(entries[3]),int(entries[5])):
                    duplicate_count += 1
                else:
                    if line_num > 1:
                        from_vertex = int(entries[3])
                        to_vertex = int(entries[5])
                        G.add_edge(from_vertex, to_vertex )
                        G.edge[from_vertex][to_vertex]['IATAFrom'] = entries[2]
                        G.edge[from_vertex][to_vertex]['IATATo'] = entries[4]
                        edge_count += 1
            except ValueError:
                # The value doesn't exist
                error_count += 1
                pass
            line_num += 1
    
    print("\t\t\t\t\t\t[Done]")

    # Limit to the first subgraph
    print("\tFinding largest subgraph",end="")
    undirected = G.to_undirected()
    subgraphs = nx.connected_component_subgraphs(undirected)
    subgraph_nodes = subgraphs[0].nodes()
    to_remove = list()
    for node in G.nodes():
        if node not in subgraph_nodes:
            to_remove.append(node)
    G.remove_nodes_from(to_remove)
    print("\t\t\t\t[Done]")

    
    print("\tRemoving isolated vertices",end="")
    # Remove nodes without inbound edges
    indeg = G.in_degree()
    outdeg = G.out_degree()
    to_remove = [n for n in indeg if (indeg[n] + outdeg[n] < 1)] 
    G.remove_nodes_from(to_remove)
    print("\t\t\t\t[Done]")



    # Calculate the edge weights
    print("\tCalculating edge weights",end="")
    G = calculate_weights(G)
    print("\t\t\t\t[Done]")
    
    # Add clustering data
    print("\tCalculating clustering coefficents",end="")
    cluster_network = nx.Graph(G)
    lcluster = nx.clustering(cluster_network)
    for i,j in G.edges():
        cluster_sum = lcluster[i] + lcluster[j]
        G[i][j]['cluster'] = cluster_sum
    print("\t\t\t[Done]")

    # Flag flights as domestic or international.
    print("\tCategorizing international and domestic flights",end="")
    for i,j in G.edges():
        if G.node[i]["country"] == G.node[j]['country']:
            G[i][j]['international'] = False
        else:
            G[i][j]['international'] = True
        if G.node[i]["country"] == G.node[j]['country']:
            G[i][j]['domestic'] = True
        else:
            G[i][j]['domestic'] = False
    print("\t\t[Done]")

    return G

def calculate_weights(input_network):
    """
    Add weights to the edges of a network based on the degrees of the connecting
    verticies, and return the network.

    Args:
        input_network: A NetworkX graph object
    Returns:
        G: A weighted NetworkX graph object.
    """
    
    G = input_network.copy()

    # Add weights to edges
    for node in G.nodes():
        successors = G.successors(node)
        weights = dict()

        # Calculate the total out degree of all succs
        total_degree = 0
        for successor in successors:
  
            try:
                total_degree += G.out_degree(successor)
            except TypeError:
                # Don't add anything
                pass

        # Find the weight for all possible successors
        for successor in successors:
            successor_degree = G.out_degree(successor)

            try:
                int(successor_degree)
            except TypeError:
                successor_degree = 0

            if total_degree > 0:
                probability_of_infection = successor_degree / \
                                           total_degree
            else:
                probability_of_infection = 0

            weights[successor] = probability_of_infection
        
        largest_weight = 0
        smallest_weight = 2
        for successor, weight in weights.items():
            if weight > largest_weight:
                largest_weight = weight
            elif weight < smallest_weight:
                smallest_weight = weight
        #(strat.shared_fitness - lowest_fitness) / \
        #                       (highest_fitness - lowest_fitness)

        for successor in successors:
            if largest_weight != smallest_weight:
                relative_weight = (weights[successor] - smallest_weight) /\
                                  (largest_weight - smallest_weight)
            else:
                relative_weight = 0
            G[node][successor]['weight'] = relative_weight

    return G

def infection(input_network, vaccination, starts,DELAY=0, vis = False, 
              file_name = "sir.csv", title="", inf_type=False, RECALCULATE = True):
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
    
    # Recalculate the weights of the network as per necessary

    # Open the data file
    f = open(file_name, "w")
    f.write("time, s, e, i, r\n")

    # Set the default to susceptable
    sys.stdout.flush()
    for node in network.nodes():
        network.node[node]["status"] =  "s"
        network.node[node]["color"] = "#A0C8F0"
        network.node[node]["age"] = 0
    
    # Assign the infected
    for start in starts:
        infected = start
        network.node[infected]["status"] = "i"
        network.node[infected]["color"]  = "green"
        
        if isinstance(network,nx.DiGraph):
            in_degree = network.in_degree()[infected] 
            out_degree = network.out_degree()[infected]
            degree = in_degree + out_degree
        else:
            degree = network.degree()[infected]

        print("\t",network.node[infected]["name"],"[",degree,"]")


    if vaccination is not None:
        print("\tVaccinated: ", len(vaccination) )
    else: 
        print("\tVaccinated: None")

    if vis:
        pos = nx.spring_layout(network, scale=2)

    # Iterate through the evolution of the disease.
    for step in range(0,99):
        # If the delay is over, vaccinate.
        # Convert the STRING! 
        if int(step) == int(DELAY):
            if vaccination is not None:
                network.remove_edges_from(vaccination)
                # Recalculate the weights of the network as per necessary
                if RECALCULATE == True:
                    network = calculate_weights(network)


        # Create variables to hold the outcomes as they happen
        S,E,I,R = 0,0,0,0

        for node in network.nodes():
            status = network.node[node]["status"]
            age = network.node[node]["age"]
            color = network.node[node]["color"]

            if status is "i" and age >= 11:
                # The infected has reached its recovery time
                network.node[node]["status"] = "r"
                network.node[node]["color"] = "purple"
                
            if status is "e" and age >= 3 and age < 11:
                # The infected has reached its recovery time
                network.node[node]["status"] = "i"
                network.node[node]["color"] = "green"

            elif status is "e":
                network.node[node]["age"] += 1

            elif status is "i":
                # Propogate the infection.
                if age > 0:
                    victims = network.successors(node)
                    number_infections = 0
                    for victim in victims:
                        infect_status = network.node[victim]["status"]
                        infect = False # Set this flag to False to start 
                                       # weighting.


                        if random.uniform(0,1) <= network[node][victim]['weight']:

                            infect = True
                            number_infections+=1

                        if infect_status == "s" and infect == True:
                            network.node[victim]["status"] = "e"
                            network.node[victim]["age"] = 0
                            network.node[victim]["color"] = "#FF6F00"
                network.node[node]["age"] += 1


        # Loop twice to prevent bias.
        for node in network.nodes():
            status = network.node[node]["status"]
            age = network.node[node]["age"]
            color = network.node[node]["color"]

            if status is "s":
                # Count those susceptable
                S += 1

            if status is "e":
                E += 1

            if status is "v":
                S += 1

            elif status is "r":

                R += 1

            elif status is "i":
                
                I += 1
        print("{0}, {1}, {2}, {3}, {4}".format(step, S, E, I, R))

        printline = "{0}, {1}, {2}, {3}, {4}".format(step, S, E, I, R)
        f.write(printline + "\n")

       # print("\t"+printline)

        if I is 0:
            break

        if vis:
            #write_dot(network, title+".dot")
            visualize(network, title, pos)
        
    print("\t----------\n\tS: {0}, I: {1}, R: {2}".format(S,I,R))

    return {"Suscceptable":S,"Infected":I, "Recovered":R}

       
def visualize(network, title,pos):
    """
    Visualize the network given an array of posisitons.
    """
    print("-- Starting to Visualize --")
    MAP = False

    if MAP:
        m = Basemap(
            projection='cea',
            llcrnrlat=-90, urcrnrlat=90,
            llcrnrlon=-180, urcrnrlon=180,
            resolution=None
            )

        pos = dict()

        for pos_node in network.nodes():
            # Normalize the lat and lon values
            x,y = m(float(network.node[pos_node]['lon']),
                    float(network.node[pos_node]['lat']))
        
            pos[pos_node] = [x,y]


    colors = []
    i_edge_colors = []
    d_edge_colors = []
    default = []
    infected = []
    for node in network.nodes():
        colors.append(network.node[node]["color"])
    for i,j in network.edges():
        color = network.node[i]["color"]
        alpha = 0.75
        if color == "#A0C8F0" or color == "#FF6F00" or color == "purple":
            color = "#A6A6A6"
            default.append((i,j))
            d_edge_colors.append(color)
        else:
            color = "#29A229"
            infected.append((i,j))
            i_edge_colors.append(color)

    plt.figure(figsize=(7,7))

    # Fist pass - Gray lines
    nx.draw_networkx_edges(network,pos,edgelist=default,
            width=0.5,
            edge_color=d_edge_colors,
            alpha=0.5,
            arrows=False)
   
    # Second Pass - Colored lines
    nx.draw_networkx_edges(network,pos,edgelist=infected,
            width=0.5,
            edge_color=i_edge_colors,
            alpha=0.75,
            arrows=False)

    nx.draw_networkx_nodes(network,
            pos,
            linewidths=0.5,
            node_size=10,
            with_labels=False,
            node_color = colors)
    
    # Adjust the plot limits
    cut = 1.05
    xmax = cut * max(xx for xx,yy in pos.values())
    xmin =  min(xx for xx,yy in pos.values())
    xmin = xmin - (cut * xmin)


    ymax = cut * max(yy for xx,yy in pos.values())
    ymin = (cut) * min(yy for xx,yy in pos.values())
    ymin = ymin - (cut * ymin)

    plt.xlim(xmin,xmax)
    plt.ylim(ymin,ymax)

    if MAP:
        # Draw the map
        m.bluemarble()
    plt.title=title

    plt.axis('off')

    number_files = str(len(os.listdir()))
    while len(number_files) < 3:
        number_files = "0" + number_files

    plt.savefig("infection-{0}.png".format(number_files),
                bbox_inches='tight', dpi=600 
            )
    plt.close()

def write_dot(network, path):
    """
    Write a DOT file for external network rendering. This function is necessary
    due to networkx incompentence.

    Args:
        network: A NetworkX graph object.
        path: A local path for the output file.
    """

    dot_file = open(path, "w")

    # Write header info
    dot_file.write("digraph network {\n")
    
    for node in network.nodes():
        neighbors = network.neighbors(node)
        for neighbor in neighbors:
            dot_file.write("\t{0} -> {1};\n".format(node, neighbor))

    dot_file.write("}")

if __name__ == "__main__":
    main()
