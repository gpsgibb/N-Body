#include <iostream>
#include <cstdlib>
#include "state.hpp"
#include "kdtree.hpp"
#include "integrator.hpp"
#include "iterator.hpp"
#include "forces.hpp"
#include "io.hpp"

int main(int argc, char** argv){

    int n = 1000;

    srand(0);

    const std::size_t NDIM = 3;
    double dt = 0.001;

    State<NDIM> state(n);
    for (std::size_t i=0;i<state.particles.size();i++){
        for (std::size_t j=0;j<NDIM;j++){
            state.particles[i].pos[j] = ((double) rand())/RAND_MAX;
        }
        state.particles[i].mass = 1.0;
    }

    KDTree<NDIM> tree(state.particles);

    KDTreeGravForce<NDIM> force(state.particles, tree);
    BlockTimestepIntegrator<NDIM> integrator(state, force, dt);
    Iterator<NDIM> iterator(state, integrator, 1000);

    iterator.start();

    return 0;

}