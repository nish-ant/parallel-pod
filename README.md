**MPI-based parallel proper orthogonal decomposition (POD) in C++**

Adapted from: [marrov/parallel-pod](https://github.com/marrov/parallel-pod)

## Build
After cloning the repository, run the following from the terminal:

        cd path/to/rom4wt/POD
        mkdir cmake-build-release
        cd cmake-build-release
        cmake -DCMAKE_BUILD_TYPE=Release ../
        make -j 6 

## Test
Refer an OpenFOAM test case in `test/example.laminarVortexShedding` for an exmple of POD calculation using snapshot data.