#ifndef ITERATOR_H
#define ITERATOR_H

#include <iostream>
#include "integrator.hpp"


template<std::size_t NDIM>
class Iterator{
    private:
        int maxits;
        State<NDIM> &state;
        Integrator<NDIM> &integrator;
    
    public:
        Iterator(State<NDIM> &s, Integrator<NDIM> &intgrtr, int mxits);
        void start();
};

#endif