#ifndef STATE_H

#define STATE_H

#include <vector>
#include "particle.hpp"
#include <string>
#include <fstream>
#include <iostream>

template<std::size_t NDIM>
class State{
    public:
        double t;
        int it;
        std::vector<Particle<NDIM>> particles; 
    public: 
        State(int n){
            t = 0.0;
            it = 0;
            particles.resize(n);
        }
};

#endif