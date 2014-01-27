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

import getopt
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
import random
import sys 

def main(argv):
    """
    Primary function that initiates network creation and handles execution of
    infection simulations.

    Args:
        argv: A list of command-line arguments passed to the application.

    Returns:
        Void

    """

    # Create the network using the command arguments.
    network = create_network(argv[0], argv[1])

    # Simulate the infection
    infection(network, None, 3682)

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

    print("\tLoading airports.", end="")
    # Populate the graph with nodes.
    with open(nodes) as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")

           # if entries[3] == "United States":
           #     G.add_node(entries[0], 
           #                name=entries[1], 
           #                lat=entries[6],
           #                lon=entries[7])
            G.add_node(entries[0], 
                       name=entries[1], 
                       lat=entries[6],
                       lon=entries[7])


    print("\t\t\t\t\t[Done]")
    
    
    print("\tLoading routes.", end="")
    # Populate the graph with edges.
    with open(edges) as f:

        for line in f.readlines():
            entries = line.replace('"',"").rstrip().split(",")
            
            if entries[3] in G.nodes() and entries[5] in G.nodes():
                G.add_edge(entries[3], entries[5] )

    
    print("\t\t\t\t\t\t[Done]")

    # Remove nodes without inbound edges
    deg = G.in_degree()
    to_remove = [n for n in deg if deg[n] < 1]
    G.remove_nodes_from(to_remove)

    return G

def infection(network, vaccination, start):
    """
    Simulate an infection within network, generated using seed, and with the
    givin vaccination strategy. This function will write data from each timestep
    to "data.csv".

    Args:
        network: A NetworkX DiGraph object.
        vaccination: A list of node indices to label as recovered from the 
                     begining.

    Returns:
        void

    """

    print("Simulating infection.")

    # Open the data file
    f = open("data.csv", "w")
    f.write("time, s, i, r\n")

    # Set the default to susceptable
    for node in network.nodes():
        network.node[node]["status"] =  "s"
        network.node[node]["color"] = "b"
        network.node[node]["age"] = 0
    
    # Add in the recovering
    if vaccination is not None:
        for node in vaccination:
            network.node[node]["status"] = "r"
            network.node[node]["color"] = "g"

    # Assign the infected
    infected = str(start)
    network.node[infected]["status"] = "i"
    network.node[infected]["color"]  = "orange"

    print("\tInitial vector: "+network.node[infected]["name"])
    #print("\t\t Degree: "+str(network.degree(infected)))
    #print("\t\t Betweenness: "+str(nx.betweenness_centrality(network)[infected]))

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
                            network.node[infectees]["color"] = "r"
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

        print("\t"+printline)

        if I is 0:
            break

        colors = []
        for node in network.nodes():
            colors.append(network.node[node]["color"])

        pos = nx.spring_layout(network)

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
        plt.show()


if __name__ == "__main__":
       main(sys.argv[1:])
