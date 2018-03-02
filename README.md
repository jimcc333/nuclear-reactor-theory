# Nuclear-Reactor-Theory
Code used in my Nuclear Reactor Theory class to simulate neutron flux under various conditions and assumptions. All files can be used with the default parameters to quickly see their funcitonality.

This code is used to understand how flux and fission rate density behave under various nuclear reactor assumptions/conditions. The functions in this repository increase in complexity for the underlying neutron behavior in the following order, where the first is the simplest:
- lecture9.m
- Cylinder3.m
- Cylinder3Crit.m
- Cylinder3MG.m

### lecture9.m
Performs a simple matrix solution to a 1D neutron flux problem. The solution is simple enough to have a mathematical function that describes it is a hyperbolic cosine.

### Cylinder3.m
This program solves the reactor eqution for a 3-region infinite cylinder where the two inner zones undergo fission and the outer region does not. The provided material properties for the inner two regions are for uranium (with user-defined enrichment) and the outer zone is a lead reflector. The dimensions of the zones are also specified by the user. The program returns the flux profile and multiplication factor as well as plotting the final neutron flux.

### Cylinder3Crit.m
Cylinder3Crit.m uses most of the workflow of Cylinder3.m to find an enrichment that makes the setup critical. It does this by solving the reactor eqution for a 3-region infinite cylinder where the two inner zones undergo fission and the outer region does not. The provided material properties for the inner two regions are for uranium (with user-defined enrichment) and the outer zone is a lead reflector. The dimensions of the zones are also specified by the user. The program returns the needed enrichment for criticality, as well as the flux profile. 

### Cylinder3MG.m
Cylinder3MG.m builds on Cylinder3.m by adding 4 energy groups. It solves the 4-group reactor eqution for a 3-region infinite cylinder where the two inner zones undergo fission and the outer regiondoes not. The provided material properties for the inner two regions are for uranium (with user-defined enrichment and water fraction) and the outer zone is a lead reflector. The dimensions of the zones are also specified by the user. The program returns the flux profile and multiplication factor.
