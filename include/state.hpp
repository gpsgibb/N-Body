#ifndef STATE_H

#define STATE_H

#include <vector>
#include "particle.hpp"

template<std::size_t NDIM>
class State{
    public:
        double time;
        std::vector<Particle<NDIM>> particles; 
    public: 
        State(int n){
            time = 0.0;
            particles.resize(n);
        }
};

#endif