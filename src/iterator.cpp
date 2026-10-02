#include <iostream>
#include "iterator.hpp"


template<std::size_t NDIM>
Iterator<NDIM>::Iterator(State<NDIM> &s, Integrator<NDIM> &intgrtr, int mxits): state(s), integrator(intgrtr) {
    it = s.it;
    maxits = mxits;
}


template<std::size_t NDIM>
void Iterator<NDIM>::start(){
    while (it < maxits){
        std::cout << "Step " << it << " of " << maxits << std::endl;
        integrator.step();
        it += 1;
    }
}


template class Iterator<2>;
template class Iterator<3>;