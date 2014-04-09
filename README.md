# Controlling Disease Propagation Through an International Air Traffic Network

## Abstract

Airplanes serve as vital transportation in a highly globalized and
interconnected economy. Unfortunately, airplanes make it extremely easy for
specific pathogens to propagate quickly and effectively to distant parts of
the world. To model the spread of a pathogen through air routes, a directed
network was generated using existing airports as vertices and their associated
routes as edges. Airports were infected for a certain amount of time before
they recover. Multiple strategies to decrease the number of infected airports
were studied at a variety of efforts, or the percentage of total flights to be
canceled, to determine the most appropriate flights to cancel. We found that
by canceling flights based on edge betweenness centrality and clustering
coefficient resulted in a significant reduction in the number of airports
harboring infectious individuals. We also examined the effect of delaying our
cancellation strategy on the number of airports harboring infectious
individuals. We found that the greatest increase in the number of infections
occurs after one week of delay, although any delay causes a statistically
significant increase in the number of infections. Based on these results we
suggest that the World Health Organization and the IATA develop an adequate
response plan to cancel flights based on betweenness centrality and clustering
coefficient in the event of a serious epidemic.

## Introduction

## Methods

Using python and NetworkX, a network consisting of airports as nodes and
airline routes as edges was constructed with data from
[openflights.org](http://openflights.org/data) as shown in Figure 1. Nodes
without inbound or outbound edges were removed. In addition, redundant edges
were removed from the network.
[Documentation](./documentation/flu_simulator.html) on this simulator is
provided in pydoc form.

<center>
<img align=center width=600 src=./documentation/images/network_plot.png />
<div style="text-align:justify;width:550px;font-size:0.8em;"><b>Figure 1: Plot 
of the modified network.</b> The light blue vertices represent susceptible 
airports. Simulations start with 10 airports with infectious individuals. The 
edges these infividuals can travel on are colored green.</div>
</center>

Examination of the properties of our network yielded 39,467 edges, 3,308 
vertices and has a diameter of three. We also examined the degree distribution
for the airports in the network. As shown in Figure 2, the distribution follows
a power law, where the vast majority of airports worldwide are smaller regional
airports with 2 to 20 inbound and outbound flights. Additionally, there are few
airports that have between 300 and 400 inbound and outbound flights.

<center>
<img align=center width=600 src=./documentation/images/degree_disributions/network-dist.png />
<div style="text-align:justify;width:550px;font-size:0.8em;"><b>Figure 2: Degree
distribution for the airport network.</b> The network exhibits a power law, with
over 700 airports with a degree of 3 or less. Additionally, there are very few 
airports with a degree higher than 400. This implies that the topology of the 
network favors using a small number of international hubs to carry traffic 
between a large number of smaller regional airports.</div>
</center>

We also examined the degree distribution of both a high and low degree airport's
neighbors to better understand the sort of airports hub airports and regional 
airports are connected to. As seen in Figure 3, the distributions for large 
airports and small airports are highly distinct. In the case of the degree 
distribution for Hartfield-Jackson Atlanta International Airport, ATL is 
connected to a wide variety of airports, including low-degree airports such as
regional airports, as well as other high-degree hub airports. In contrast, 
Greater Rochester International Airport, a low-degree regional airport 
regardless of name, is connected almost entirely to airports with a degrees 
between 200 and 500.

<center>
<img align=center width=600 src=./documentation/images/degree_disributions/atl-dist.png />
<img align=center width=600 src=./documentation/images/degree_disributions/roc-dist.png />
<div style="text-align:justify;width:550px;font-size:0.8em;"><b>Figure 2: Degree
distribution for the airport network.</b> The network exhibits a power law, with
over 700 airports with a degree of 3 or less. Additionally, there are very few 
airports with a degree higher than 400. This implies that the topology of the 
network favors using a small number of international hubs to carry traffic 
between a large number of smaller regional airports.</div>
</center>

We implemented a disease based upon influenza A. As such, there were four states,
susceptible, exposed, infectious, and recovered. The infection took three days to
incubate and 7 days to recover from, and as such was observed to have a basic 
reproductive rate of 2.32. Using the aforementioned example disease, 
the model propagated the disease by following four rules:

  1. One plane leaves on every outbound route every tick and arrives at its 
        destination.
  2. Airports are used as proxies for individuals harbored within the airport.
        This assumption is held witht he understanding that transient passengers
        are able to infect perminant employees of the airport. These employees
        are then able to pass the infection on to other individuals in the 
        airport.
  3. Routes are directional, and each edge is weighted to represent the 
        probability of carrying infectious individuals based upon the degree of 
        the source and destination airports.
  4. Edge weights are recalculated after each cancellation strategy to properly
        model the flow of individuals around our cancellations. 

Propagation of the infection through the network can be seen in Figure 4. 
Starting with 10 randomly chosen airports of high degree, infectious 
individuals, as depicted in green, propagate along the green-colored edges. As
time continues, the majority of airports infected are the highly connected 
airports in the center of the graph. This can be compared to the less-connected
regional or public airports seen on the periphery of the graph. The simulation
seen in Figure 4 suggests that infections are most likely to occur in the larger
international airports, making them valid targets for cancellations and 
quarantine.

<center>
<img align=center width=600 src=./documentation/images/infection.gif />
<div style="text-align:justify;width:550px;font-size:0.8em;"><b>Figure 4: Degree
distribution for the airport network.</b> The network exhibits a power law, with
over 700 airports with a degree of 3 or less. Additionally, there are very few 
airports with a degree higher than 400. This implies that the topology of the 
network favors using a small number of international hubs to carry traffic 
between a large number of smaller regional airports.</div>

</center>

## Results

## Discussion
