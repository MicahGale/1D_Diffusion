# 1D_Diffusion
A 1-dimensional diffusion neutron transport solver. Solves the multi-group diffusion equation for neutron transport for a 1-D slab problem. 
## Solution Method
1. user specificies an array of material cells with the required material properties. 
1. The input is automatically expanded into a series of slabs with the desired mesh spacing.
1. The non-fission matrix(H) and an initial flux guess are generated.
1. The fission matrix(F) is generated based on the flux guess
1. The problem HΦ=FΦ is inverted to solve for Φ, the flux. 
1. This process is iterated until the flux converges.
