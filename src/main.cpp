#include <iostream>
#include <omp.h>
#include "state.hpp"
#include "kdtree.hpp"
#include <cstdlib>
#include <vector>
#include <integrator.hpp>
#include <iterator.hpp>

int main(int argc, char** argv){

    int n = 1000;

    srand(0);

    const std::size_t NDIM = 3;
    double dt = 0.001 * 0.03;

    State<NDIM> state(n);
    for (std::size_t i=0;i<state.particles.size();i++){
        for (std::size_t j=0;j<NDIM;j++){
            state.particles[i].pos[j] = ((double) rand())/RAND_MAX;
        }
        state.particles[i].mass = 1.0;
    }

    AdaptiveTimestepVelocityVerletIntegrator<NDIM> integrator(state, dt);
    Iterator<NDIM> iterator(state, integrator, 10000);

    iterator.start();

    state.save_state("state.dat");

    return 0;

}