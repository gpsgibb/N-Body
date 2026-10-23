#include <iostream>
#include "iterator.hpp"
#include "omp.h"

template<std::size_t NDIM>
Iterator<NDIM>::Iterator(State<NDIM> &s, Integrator<NDIM> &intgrtr, int mxits): state(s), integrator(intgrtr), writer(Writer<NDIM>(s)) {
    maxits = mxits;
}


template<std::size_t NDIM>
void Iterator<NDIM>::start(){
    double t0, t1, t00, t11;
    t0 = omp_get_wtime();
    while (state.it < maxits){
        t00 = omp_get_wtime();
        state.it += 1;
        std::cout << "Iteration " << state.it << " of " << maxits << std::endl;
        integrator.step();
        writer.write();
        t11 = omp_get_wtime();
        std::cout << "Time taken for iteration: " << t11-t00 << " s" << std::endl;
    }
    t1 = omp_get_wtime();
    std::cout << "Total runtime = " << t1-t0 << "s" << std::endl;
}


template class Iterator<2>;
template class Iterator<3>;