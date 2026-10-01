#include <iostream>
#include <omp.h>
#include "state.hpp"
#include "kdtree.hpp"
#include <cstdlib>
#include <vector>
#include <integrator.hpp>

int main(int argc, char** argv){

    int n = 1000;

    srand(0);

    const std::size_t NDIM = 2;
    double dt = 0.001 * 0.03;

    State<NDIM> s2(n);
    for (std::size_t i=0;i<s2.particles.size();i++){
        for (std::size_t j=0;j<NDIM;j++){
            s2.particles[i].pos[j] = ((double) rand())/RAND_MAX;
        }
        s2.particles[i].mass = 1.0;
    }

    double t0, t1;
    t0 = omp_get_wtime();
    KDTree<NDIM> tree(s2.particles);
    t1 = omp_get_wtime();
    std::cout << "Construction time taken " << t1-t0 << "s\n";

    // get forces
    t0 = omp_get_wtime();
    #pragma omp parallel for
    for (std::size_t i=0; i<s2.particles.size(); i++){
            s2.particles[i].acc = tree.evaluate_force(i);
        }

    t1 = omp_get_wtime();
    std::cout << "Force time taken " << t1-t0 << "s\n";

    tree.save_tree("tree.dat");
    s2.save_state("state.dat");

    EulerIntegrator<2> integrator(s2, dt);
    for (int i=0; i<10000;i++){
        std::cout << i << "\n";
        integrator.step();
    }
    s2.save_state("state.dat");

    return 0;

}