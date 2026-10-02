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
            it = 0.0;
            particles.resize(n);
        }

        void save_state(std::string filename){
            std::ofstream out(filename, std::ios::binary);
            if (!out) {
                std::cerr << "Failed to open file\n";
                throw 1;
            }

            std::size_t ndim = NDIM, n=particles.size(), size = sizeof(Particle<NDIM>);

            out.write((char*) &ndim, sizeof(ndim));
            out.write((char*) &n, sizeof(std::size_t));
            out.write((char*) particles.data(), n * size);
        }
};

#endif