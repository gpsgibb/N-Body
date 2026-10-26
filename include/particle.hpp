#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>

template <std::size_t NDIM>
struct Particle{
    double mass;
    std::array<double, NDIM> pos;
    std::array<double, NDIM> vel;
    std::array<double, NDIM> acc;
    double potential_energy;
    double kinetic_energy;
    double density;
    double pressure;
    double smoothing_length;
    int tag;
};

#endif
