 This repo contains code to solve **maximum flow problem** with Dijkstra's algorithm (for finding the shortest paths between nodes in a weighted graph) in a bidirectional Manhattan street mesh topology for different combinations of 
 - number of nodes (N) and
 - number of transmitters/receivers per node (Δ)
 
 (for both cases when flow splitting is allowed and not allowed)
 
 ```
 e.g. scenario 1: N = 16;      Δ = 4
 
      scenario 2: N = 8,10,12; Δ = 4
 
      scenario 3: N = 8;       Δ = 1,2,3,4,5 
      
      and so on
 ```
 A uniform traffic matrix is considered, in which the traffic sent from any source to any destination is a uniform random variable in the range [1,10],
 i.e., traffic sent from node s to node d can be expressed as tsd = Uniform[1,10].   
 
 Maximum flow is plotted for different scenarios obtained with different seeds to generate the traffic matrix.
