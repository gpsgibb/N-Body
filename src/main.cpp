#include <iostream>
#include <omp.h>
#include "state.hpp"
#include "kdtree.hpp"


int main(int argc, char** argv){

    int n = 1000000;

    State<3> s2(n);
    for (int i=0;i<s2.particles.size();i++){
        s2.particles[i].pos[0] = i;
        s2.particles[i].pos[1] = n-i;
    }

    double t0, t1;
    t0 = omp_get_wtime();
    KDTree<3> tree(s2.particles);
    t1 = omp_get_wtime();
    std::cout << "Time taken " << t1-t0 << "s\n";

    return 0;

}