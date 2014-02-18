#!/usr/bin/python


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
        -d: Run a degree-based quarantine simulation.
        -r: Run a random quarantine simulation.
        -u: Convert the network to an undirected network.
        -s: Run a naive simulation and output the SIR data.

"""

# Title:  air_travel_fl.py
# Author: Nicholas A. Yager
# Date:   2013-01-12

import copy
import getopt
import math
import networkx as nx
import operator
import os
import random
import sys 
import time


def main():
    """
    Primary function that initiates network creation and handles execution of
    infection simulations.

    Args:
        argv: A list of command-line arguments passed to the application.

    Returns:
        Void

    """
   
    # Determine the parameters of the current simulation.
    opts, args = getopt.getopt(sys.argv[1:], "bgdrus", ["Betweenness",
                                                            "Genetic",
                                                            "Degree",
                                                            "Random",
                                                            "Undirect",
                                                            "SIR"])

    GENETIC = False 
    BETWEENNESS = False 
    DEGREE = False 
    RANDOM = False 
    SIR = False
    UNDIRECT = False

    for o, a in opts:
        if o == "-g":
            GENETIC = True
            BETWEENNESS = False
            DEGREE = False
        elif o == "-b":
            GENETIC = False
            BETWEENNESS = True
        elif o == "-d":
            GENETIC = False
            DEGREE = True
        elif o == "-r":
            GENETIC = False
            RANDOM = True
        elif o == "-s":
            SIR = True
        elif o == "-u":
            UNDIRECT = True
            
    AIRPORT_DATA = args[0]
    ROUTE_DATA = args[1]

    NUM_SIMULATIONS = 50

    seed = 100

    random.seed(seed)

    # Create the network using the command arguments.
    network = create_network(AIRPORT_DATA, ROUTE_DATA)
    target = random.sample(network.nodes(), NUM_SIMULATIONS)

    if UNDIRECT:
        network = network.to_undirected()

    # Make a directory for the data, and change into that directory.
    currenttime = time.strftime("%Y-%m-%dT%H%M%S", time.gmtime())
    os.makedirs(currenttime)
    os.chdir(currenttime)

    # Record relevent data about the simulation.
    simulation_data(network, currenttime, target, seed)

    if GENETIC:
        genetic_simulations(network, 50, 50, target)

    if BETWEENNESS:
        betweenness_simulations(network, target)
    
    if DEGREE:
        degree_simulations(network, target)

    if RANDOM:
        random_simulations(network, target)

    if SIR:
        sir_simulations(network, target)
        

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

def sir_simulations(network, targets):
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

        sir_file = "sir/sir_{0}.csv".format(pad_string(iteration,4))

        results = infection(network, None, target, False, sir_file)
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

def random_simulations(network, targets):
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

    iteration = 0
    for target in targets:

        # Open a file for this targ'ets dataset
        random_file = open("random/random_{0}.csv".format(pad_string(iteration,4)),"w")
        random_file.write('"effort","total_infected"\n')


        # Generate a baseline
        results = infection(network, None, target, False)
        total_infected = results["Infected"] + results["Recovered"]
        random_file.write("{0},{1}\n".format(0,total_infected))

        randoms = random.sample(network.nodes(), len(network.nodes()))

        # Perform a check for every strategy
        for effort in range(1,101,5):
            max_index = int(len(randoms) * (effort/100))-1
            strategy = [x for x in randoms[0:max_index]]

            results = infection(network, strategy, target, False)
            total_infected = results["Infected"] + results["Recovered"]
            random_file.write("{0},{1}\n".format(effort/100,total_infected))
        
        iteration += 1
        random_file.close()



def degree_simulations(network, targets):
    """
    Simulate the spread of infection for increasing vaccination efforts by 
    quarantining airports of decreasing degree.

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
    degrees = network.degree()
    degree = sorted(degrees.keys(), key=lambda k: degrees[k], reverse=True)
    print("\t\t\t\t\t[Done]")

    print("\tHighest Degree:{0}".format(degrees[degree[0]]))

    # Make a new folder for the degree data.
    os.makedirs("degree")

    iteration = 0
    for target in targets:

        # Open a file for this targ'ets dataset
        degree_file = open("degree/degree_{0}.csv".format(pad_string(iteration,4)),"w")
        degree_file.write('"effort","total_infected"\n')


        # Generate a baseline
        results = infection(network, None, target, False)
        total_infected = results["Infected"] + results["Recovered"]
        degree_file.write("{0},{1}\n".format(0,total_infected))

        # Perform a check for every strategy
        for effort in range(1,101,5):
            max_index = int(len(degree) * (effort/100))-1
            strategy = [x for x in degree[0:max_index]]

            results = infection(network, strategy, target, False)
            total_infected = results["Infected"] + results["Recovered"]
            degree_file.write("{0},{1}\n".format(effort/100,total_infected))
        
        iteration += 1
        degree_file.close()


def betweenness_simulations(network,targets):
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
    betweennesses = nx.betweenness_centrality(network)
    betweenness = sorted(betweennesses.keys(), 
                    key=lambda k: betweennesses[k], reverse=True)

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
        results = infection(network, None, target, False)
        total_infected = results["Infected"] + results["Recovered"]
        betweenness_file.write("{0},{1}\n".format(0,total_infected))

        # Perform a check for every strategy
        for effort in range(1,101,5):
            max_index = int(len(betweenness) * (effort/100))-1
            strategy = [x for x in betweenness[0:max_index]]

            results = infection(network, strategy, target, False)
            total_infected = results["Infected"] + results["Recovered"]
            betweenness_file.write("{0},{1}\n".format(effort/100,total_infected))

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


    # Limit to the first subgraph
    undirected = G.to_undirected()
    subgraphs = nx.connected_component_subgraphs(undirected)

    subgraph_nodes = subgraphs[0].nodes()

    to_remove = list()

    for node in G.nodes():
        if node not in subgraph_nodes:
            to_remove.append(node)
            
    G.remove_nodes_from(to_remove)

    return G

def infection(input_network, vaccination, start, visualize = False, file_name = "sir.csv"):
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
    f = open(file_name, "w")
    f.write("time, s, i, r\n")

    # Set the default to susceptable
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
                remove_successors = [ (suc, node) for suc in network.successors(node)]
                network.remove_edges_from(remove_predecessors)
                network.remove_edges_from(remove_successors)
                network.node[node]["status"] = 'v'

    # Assign the infected
    infected = start
    network.node[infected]["status"] = "i"
    network.node[infected]["color"]  = "orange"

    in_degree = network.in_degree()[infected] 
    out_degree = network.out_degree()[infected]
    print("\tInitial vector: "+network.node[infected]["name"])
    print("\tIn degree: ",in_degree)
    print("\tOut degree:",out_degree)
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
                    possible_victims = network[node]

                    victims = possible_victims
                    
                    # Calculate the total out degree of all infectees
                    total_degree = 0
                    for victim in victims:
                        total_degree += network.out_degree(victim)

                    # Infect possible victims
                    for infectees in victims:
                        infect_status = network.node[infectees]["status"]
                        victim_degree = network.out_degree(infectees)

                        infect = True # Set this flag to False to start 
                                      # weighting.

                        if total_degree > 0:
                            probability_of_infection = victim_degree / \
                                                       total_degree
                        else:
                            probability_of_infection = 0

                        if random.uniform(0,1) <= probability_of_infection:

                            infect = True

                        if infect_status == "s" and infect == True:
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

            if status is "v":
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
    main()
