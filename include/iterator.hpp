#ifndef ITERATOR_H
#define ITERATOR_H

#include <iostream>
#include "integrator.hpp"
#include "io.hpp"
#include "config.hpp"


template<std::size_t NDIM>
class Iterator{
    private:
        int maxits;
        State<NDIM> &state;
        Integrator<NDIM> &integrator;
        Writer<NDIM> writer;

    public:
        Iterator(State<NDIM> &s, Integrator<NDIM> &intgrtr, Config &config);
        void start();
};

#endif