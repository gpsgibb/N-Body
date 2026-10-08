#include <iostream>
#include "iterator.hpp"
#include "omp.h"

template<std::size_t NDIM>
Iterator<NDIM>::Iterator(State<NDIM> &s, Integrator<NDIM> &intgrtr, int mxits): state(s), integrator(intgrtr), writer(Writer<NDIM>(s)) {
    maxits = mxits;
}


template<std::size_t NDIM>
void Iterator<NDIM>::start(){
    double t0, t1;
    t0 = omp_get_wtime();
    while (state.it < maxits){
        state.it += 1;
        std::cout << "Iteration " << state.it << " of " << maxits << std::endl;
        integrator.step();
        writer.write();
    }
    t1 = omp_get_wtime();
    std::cout << "Total runtime = " << t1-t0 << "s" << std::endl;
}


template class Iterator<2>;
template class Iterator<3>;