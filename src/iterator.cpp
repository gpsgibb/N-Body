#include <iostream>
#include "iterator.hpp"
#include "omp.h"

template<std::size_t NDIM>
Iterator<NDIM>::Iterator(State<NDIM> &s, Integrator<NDIM> &intgrtr, int mxits): state(s), integrator(intgrtr) {
    maxits = mxits;
}


template<std::size_t NDIM>
void Iterator<NDIM>::start(){
    double t0, t1;
    t0 = omp_get_wtime();
    while (state.it < maxits){
        std::cout << "Iteration " << state.it << " of " << maxits << std::endl;
        integrator.step();
        state.it += 1;
    }
    t1 = omp_get_wtime();
    std::cout << "Total runtime = " << t1-t0 << "s" << std::endl;
}


template class Iterator<2>;
template class Iterator<3>;