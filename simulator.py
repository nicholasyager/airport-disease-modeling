#!/usr/bin/python3
"""
simulator.py is a simulator for an example infection spreading across a 
network between airports (nodes) via air travel routes (edges). The goal of this
simulation is to test edge-based quarantine strategies for the given network. 

Usage:

    simulator.py -brcsidv [--delay=<days>] [--nsim=<n>]
        <airport database> <route database>

Flags:
    -b: Run a betweenness-based quarantine simulation.
    -r: Run a random quarantine simulation.
    -c: Run a simulation based on vertex clustering coefficent.
    -s: Run a naive simulation and output the SIR data.
    -i: Filter to only quarantine international flights.
    -d: Filter to only quarantine domestic flights.
    -v: Visualize the network by plotting each time step.

Option:
    --delay=<days>  The number of days to delay a cancellation strategy.
    --nsim=<n>      The number of simulations to perform per strategy.
"""

# Title:  simulator.py
# Authors: Nicholas A. Yager and Matthew Taylor
# Date:   2013-01-12

import copy
import getopt
import math
import networkx as nx
import matplotlib.pyplot as plt
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
    VISUALIZE = False
    INTERNATIONAL = False
    DOMESTIC = False
    DELAY = 0
    NUM_SIMULATIONS = 100

    # Determine the parameters of the current simulation.
    opts, args = getopt.getopt(sys.argv[1:], "brcsidv", ["delay=",
                                                         "nsim="]
                                                            )

    # Check if the data arguments are available
    if len(args) < 2:
        print(__doc__)
        exit()


    AIRPORT_DATA = args[0]
    ROUTE_DATA = args[1]

    simulations = list()

    for o, a in opts:
        if o == "-b":
            simulations.append("betweenness")
        elif o == "-r":
            simulations.append("random")
        elif o == "-s":
            simulations.append("sir")
        elif o == "-c":
            simulations.append("clustering")
        elif o == "-v":
            VISUALIZE = True
        elif o == "-i":
            INTERNATIONAL = True
        elif o == "-q":
            DOMESTIC = True
        elif o == "-y":
            RECALCULATE = False
        elif o == "--delay":
            DELAY = int(a)
        elif o == "--nsim":
            NUM_SIMULATIONS = int(a)


            

    seed = 100
    random.seed(seed)

    # Identify the script.
    print("Air Disease Simulator 2.0.0")
    print("Created by Nicholas A. Yager and Matthew Taylor\n\n")

    # Create the network using the command arguments.
    network = create_network(AIRPORT_DATA, ROUTE_DATA)
  
    # Generate target-selection weights, and choose target vertices to infect.
    degrees = network.degree()
    weights = dict()
    for airport, degree in degrees.items():
        weights[airport] = network.out_degree(airport) +\
                           network.in_degree(airport)
    targets = list()
    for ind in range(0,NUM_SIMULATIONS):
        target_round = list()
        while len(target_round) < 10:
             chosen_airport = weighted_random(weights)
             if chosen_airport not in target_round:
                 target_round.append(chosen_airport)
        targets.append(target_round)


    # Make a directory for the data, and change into that directory.
    currenttime = time.strftime("%Y-%m-%dT%H%M%S", time.gmtime())
    os.makedirs(currenttime)
    os.chdir(currenttime)

    # Record relevent data about the simulation.
    # TMP simulation_data(network, currenttime, target, seed)


    # Prepare simulations

    # Remove some edges as per necessary from the pool of edges that can be
    # used as cancelled edges.

    edgepool = network.edges(data=True)
    if INTERNATIONAL:
        for i,j,data in edgepool:
            if data["international"] == False:
                edgepool.remove((i,j,data))
            index += 1
    elif DOMESTIC:
        for i,j,data in edgepool:
            if data["domestic"] == False:
                degrees.edgepool((i,j,data))
            index += 1


    for strategy in simulations:
    
        print("{0} Mode.".format(strategy) )

        index = 0

        # Generate a list a sorted list of flights to cancel based on the
        # strategy.
        
        cancellist = list()
        if strategy == "random":
            # Sort the edges randomly

            cancellist = random.sample(edgepool, len(edgepool))

        elif strategy == "clustering":
            # Sort the edges based on the sum of the clustering coefficent.

            sorted_cluster = sorted(edgepool, key=lambda k: k[2]['cluster'],
                            reverse=True)
            for cluster_item in sorted_cluster:
                if network[cluster_item[0]][cluster_item[1]]['cluster'] < 2:
                    if network[cluster_item[0]][cluster_item[1]]['cluster'] > 0:
                        cancellist.append((cluster_item[0], cluster_item[1]))

        elif strategy == "betweenness":
            # Sort the edges based on weighted edge-betweenness.

            betweennesses = nx.edge_betweenness_centrality(network,
                                                           weight="weight")
            cancellist = sorted(betweennesses.keys(), 
                                key=lambda k: betweennesses[k], reverse=True)


        print(cancellist[:20])
        # Make a new folder for the data.
        os.makedirs(strategy)
        
        iteration = 0
        efforts = [0]
        efforts.extend(range(1,101,5))
        for target in targets:

            # Open a file for this targets dataset
            output_file = open("{0}/{0}_{1}.csv".format(strategy,
                                                        pad_string(iteration,4)
                                                        ),"w")
            output_file.write('"effort","total_infected, edges_closed"\n')


            for effort in efforts:
                if effort != 0:
                    max_index = int(len(cancellist) * (effort/100))-1
                    cancelled = cancellist[0:max_index]
                else:
                    cancelled = None

                title = "{0} - {1}%".format(strategy, effort/100)
                results = infection(network, cancelled, target, vis=VISUALIZE,
                                    title=title, DELAY=DELAY)
                total_infected = results["Infected"] + results["Recovered"]
                output_file.write("{0},{1}\n".format(effort/100,total_infected))
                
                if total_infected == 1:
                    for remaining_effort in range(effort+5,101,5):
                        output_file.write("{0},{1}\n".format(remaining_effort/100,
                                                              total_infected))
                    break

            iteration += 1
            output_file.close()



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

def simulation_data(network, time, targets, seed):
    """
    Output various statistics of the nature of the network to a file, including
    the diameter, the number of verticies and edges, and the
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


    print("\tStoring network parameters")

    data_file = open("network.dat", "w")
    data_file.write("Simulation name: {0}\n\n".format(time))
    data_file.write("Network properties\n===============\n")
    data_file.write("Network type: {0}\n".format(network_type))
    data_file.write("Number of verticies: {0}\n".format(verticies))
    data_file.write("Number of edges: {0}\n".format(edges))
    data_file.write("Diameter: {0}\n".format(diameter))

    data_file.close()

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
              file_name = "sir.csv", title="",  RECALCULATE = True):
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
                print(DELAY,"on step",step)
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

if __name__ == "__main__":
    main()
