#ifndef PARTICLE_H
#define PARTICLE_H

#include <array>

template <std::size_t NDIM>
struct Particle{
    double mass;
    std::array<double, NDIM> pos;
    std::array<double, NDIM> vel;
};

#endif
