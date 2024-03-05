## Laminar vortex shedding across a 2D cylinder

Based on: https://github.com/meinarsve/CFDwOpenFoam/tree/master/LaminarVortexShedding

Adapted for comparison of parallel runs for different mesh sizes.

The example case corresponds to Re=100 and mesh size m=0.5.
- Re can be changed by changing the value of `nu` in `constant/transportProperties`.
- Mesh setup can be changed in `cylinder.m0.5.geo`.

The test can be run in two-ways:
1. Using the pre-supplied snapshots in the time range t = [150, 170].
   User needs to only run the scripts for POD calculation and reconstruction:

        sbatch runscript.pod
        sbatch runscript.pod.reconstruct

2. Starting from data generation by running the OpenFOAM simulation.
   Launch the simulation:

        sbatch runscript.solve

   and then run the two POD scripts as mentioned in (1).

**Comparison between simulated (left) and POD reconstructed (right) fields:**

![Comparison between fields](plot/Umag.field.recon.gif)