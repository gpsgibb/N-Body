#include <iostream>
#include "iterator.hpp"


template<std::size_t NDIM>
Iterator<NDIM>::Iterator(State<NDIM> &s, Integrator<NDIM> &intgrtr, int mxits): state(s), integrator(intgrtr) {
    maxits = mxits;
}


template<std::size_t NDIM>
void Iterator<NDIM>::start(){
    while (state.it < maxits){
        std::cout << "Step " << state.it << " of " << maxits << std::endl;
        integrator.step();
        state.it += 1;
    }
}


template class Iterator<2>;
template class Iterator<3>;