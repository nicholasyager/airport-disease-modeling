<style>

img {
    width:100%;
}

</style>

Airport Flu Modeling
=====



Methods
-----

A network consisting of airports as nodes and airline routes as edges was 
constructed with data from [openflights.org][openflights] as shown in Figure 1.
Nodes without inbound edges were removed.

![world-network](documentation/images/world_crop.png "Figure 1.")

Two models were created for the investigation of flu-dynamics in air travel,
namely a "simple model" and a complex model. The simple model has 4 rules:

1. One plane leaves on every outbound route every tick and arrives at its 
destination.
2. One sick plane will completely infect a destination airport.
3. Routes are directional.
4. Airports take 5 - 7 ticks to recover.


An example SIR plot of this model looks as follows

![simple-model](documentation/images/simple.png "Figure 2.")

The complex model follows 6 rules:

1. One plane leaves on each route per hour.
2. All routes are directional.
3. Each airport has 10 people per route that can be either susceptible, 
infectious, or recovered.
4. There is a parameter &beta; for the number of people an infected person can
infect per hour, and a parameter &gamma; for the number of hours required to 
recover from being infectious.
5. Each plane has 10 passengers, and will travel along its route at a speed of
&mu; over the distance between airports.
6. An infected individual is only capable of infecting other passengers of a 
plane during the flight.


This model's SIR plot is as follows

![complex-model](documentation/images/dummy.png "Figure 3.")

[openflights]: http://openflights.org/data
