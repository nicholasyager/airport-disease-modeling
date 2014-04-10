# Controlling Disease Propagation Through an International Air Traffic Network
Nicholas A. Yager <sup>1</sup> and Matthew Taylor <sup>2</sub>

<sup>1</sup> Department of Biochemistry, State University of New York at 
Geneseo, Geneseo, New York 14454<br>
<sup>2</sup> Department of Biology, State University of New York at Geneseo, 
Geneseo, New York 14454


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

<center>
<img align=center width=500 src=./documentation/images/network_plot.png />
<div style="text-align:justify;width:450px;font-size:0.8em;"><b>Figure 1: Plot 
of the modified network exhibits degree stratification..</b> The light blue 
vertices represent susceptible airports. Simulations start with 10 airports with 
infectious individuals. The edges these individuals can travel on are colored 
green. The separation of the tight clusters in the center from the rest of the
graph suggests that the network is highly stratified in the types of connections
between airports.</div>
</center>

Examination of the properties of our network yielded 39,467 edges, 3,308 
vertices and has a diameter of three. We also examined the degree distribution
for the airports in the network. As shown in Figure 2, the distribution follows
a power law, where the vast majority of airports worldwide are smaller regional
airports with 2 to 20 inbound and outbound flights. Additionally, there are few
airports that have between 300 and 400 inbound and outbound flights.

<center>
<img align=center width=500 src=./documentation/images/degree_disributions/network-dist.png />
<div style="text-align:justify;width:450px;font-size:0.8em;"><b>Figure 2: Degree
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
<img align=center width=500 src=./documentation/images/degree_disributions/atl-dist.png />
<img align=center width=500 src=./documentation/images/degree_disributions/roc-dist.png />
<div style="text-align:justify;width:450px;font-size:0.8em;"><b>Figure 3: Degree
distribution for the neighbors of ATL and ROC highlights differences in network
connectivity.</b> As shown in the degree distribution for ATL's neighbors, ATL 
is connected to a wide variety of different size airports. Although ATL is 
primarily connected to airports with a degree between 12 and 54, it is also 
connected to airports with a degree of two as well as airports with a degree 
over 400. This is contrasted to ROC, which is connected primarily to airports 
with between 200 and 400 connections. This would suggest that high-degree hubs 
like ATL serve as connecting stops for flights to and from smaller airports. 
From this, we can assume that there will be more traffic traveling from small 
airports to international hubs, than there will be from an international hub to 
a particular regional airport.</div>
</center>

To examine the propagation of disease through the network, we implemented an
example disease based upon influenza A. As such, there were four states,
susceptible, exposed, infectious, and recovered. The infection took three days
to incubate and 7 days to recover from, and as such was observed to have a basic 
reproductive rate of 2.32. Using the aforementioned example disease, 
the model propagated the disease by following five rules:

  1. The disease starts randomly from 10 airports, with a higher probability of
        starting in high-degree airports.
  2. One plane leaves on every outbound route every tick and arrives at its 
        destination.
  3. Airports are used as proxies for individuals harbored within the airport.
        This assumption is held witht he understanding that transient passengers
        are able to infect perminant employees of the airport. These employees
        are then able to pass the infection on to other individuals in the 
        airport.
  4. Routes are directional, and each edge is weighted to represent the 
        probability of carrying infectious individuals based upon the degree of 
        the source and destination airports.
  5. Edge weights are recalculated after each cancellation strategy to properly
        model the flow of individuals around our cancellations. 

To test our ability to limit the number of infections in the network, we
implemented a strategy of canceling flights based on edge betweenness 
centrality, clustering coefficient summation, and a random implementation. For
this, we calculated the betweenness centrality from each edge *i* to each other
edge *j*, where σ<sub>ij</sub> is the number of shortest paths from *i* to *j*,
such that

<center>
<img src="./documentation/equations/betweenness_centrality.png" width="146" height="53" />
</center>

The clustering coefficient summation is assigned to edge<sub>ij</sub> as the 
sum of the clustering coefficients of vertices *i* and *j*. The random strategy 
randomly closes an edge in the network.


## Results

Propagation of the infection through the network can be seen in Figure 4. 
Starting with 10 randomly chosen airports of high degree, infectious 
individuals, as depicted in green, propagate along the green-colored edges. As
time continues, the majority of airports infected are the highly connected 
airports in the center of the graph. This can be compared to the less-connected
regional or public airports seen on the periphery of the graph. The simulation
seen in Figure 4 suggests that infections are most likely to occur in the larger
international airports, making them valid targets for cancellations and 
quarantine.

<center style="page-break-before:always;">
<img align=center width=500 src=./documentation/images/infection.gif />
<div style="text-align:justify;width:450px;font-size:0.8em;"><b>Figure 4: 
Propagation of the disease favors highly connected airports near the graph
center..</b> The infection starts in the 10 most connected airports in the 
network. As a result of how network weights are applied, infections are more 
likely to be spread to other large airports. This trend is easily visualized as 
the airports in the center of the graph are infected more often. Smaller 
regional airports on the periphery of the graph are sometimes infected, but the 
rate of infection is much lower when distant from the central, highly connected 
networks. </div>
</center>

To collect a large enough sample size, 343 simulations were performed for each
strategy at 20 cancellation efforts at a 5% interval. As such a total of 20,580
simulations were performed, and the number of airports harboring infected 
individuals was used to the determination of the optimal strategy. We chose to
examine five cancellation efforts, namely 0%, 10%, 30%, 50%, and 80%. All
strategies within the five cancellation efforts were determined to be normally 
distributed using a Shaprio-Wilk test (a = 0.05). An analysis of variance was
performed, and we found that strategy (D = 4, F = 52730.4, p < 0.001 ), effort 
(D = 2, F = 640.8, p < 0.001) and their interaction 
(D = 8, F = 2037.3, p < 0.001) were all significant factors.

A TukeyHSD was performed, as shown in Figure 5. We found that in efforts under
30%, the application of a random cancellation strategy will increase the number
for airports containing infectious individuals. We also found that in high
cancellation efforts the clustering coefficient strategy resulted in the 
greatest decrease in the number of airports with infectious individuals. In
a moderate effort, such as 30%, betweenness centrality is the optimal strategy,
and is statistically different that the 30% clustering coefficient-based 
strategy. 

<center>
<img align=center width=500 src=./documentation/images/strategy_comparison.png />
<div style="text-align:justify;width:450px;font-size:0.8em;"><b>Figure 5: 
Betweenness centrality and clustering coefficient cancellation strategies prove
effective at moderate cancellation efforts.</b> </div>
</center>


## Discussion

## References
   1.  He D, Dushoff J, Eftimie R, Earn DJD. 2013 *Patterns of sperad of 
influenza A in Canada*. *Proc. R. Soc. B.* 280, 20131174; 
DOI: 10.1098/rspb.2013.1174
   2. Brockmann D, Helbing D. 2013 *The Hidden Geometry of Complex, 
Network-Driven Contagion Phenomena*. *Science* 342, 1337; 
DOI: 10.1126/science.1245200
 

