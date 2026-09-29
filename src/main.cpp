#include <iostream>
#include <omp.h>
#include "state.hpp"
#include "kdtree.hpp"
#include <cstdlib>

int main(int argc, char** argv){

    int n = 10000;

    srand(0);

    State<3> s2(n);
    for (std::size_t i=0;i<s2.particles.size();i++){
        s2.particles[i].pos[0] = ((double) rand())/RAND_MAX;
        s2.particles[i].pos[1] = ((double) rand())/RAND_MAX;
        s2.particles[i].pos[2] = ((double) rand())/RAND_MAX;
        s2.particles[i].mass = 1.0;
    }

    double t0, t1;
    t0 = omp_get_wtime();
    KDTree<3> tree(s2.particles);
    t1 = omp_get_wtime();
    std::cout << "Construction time taken " << t1-t0 << "s\n";

    // get forces
    t0 = omp_get_wtime();
    #pragma omp parallel
    {
        std::array<double, 3> force;
        #pragma omp for
        for (std::size_t i=0; i<s2.particles.size(); i++){
            force = tree.evaluate_force(i);
        }
    }
    t1 = omp_get_wtime();
    std::cout << "Force time taken " << t1-t0 << "s\n";

    return 0;

}