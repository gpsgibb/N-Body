#ifndef FORCES_H
#define FORCES_H

#include "particle.hpp"
#include "kdtree.hpp"
#include <cmath>
#include <array>

template<std::size_t NDIM>
class Force{
    public:
        Force() {};
        virtual std::array<double, NDIM> get_force(std::size_t i);
        virtual void full_update(){};
        virtual void partial_update(){};
        virtual double get_interaction_lengthscale(std::size_t) {return INFINITY;}
};


template<std::size_t NDIM>
class KDTreeGravForce: public Force<NDIM>{
    private:
        std::vector<Particle<NDIM>> &particles;
        KDTree<NDIM> &tree;
    public:
        KDTreeGravForce(std::vector<Particle<NDIM>> &particles, KDTree<NDIM> &tree);
        std::array<double, NDIM> get_force(std::size_t i) override;
        void full_update() override;
        void partial_update() override;
        double get_interaction_lengthscale(std::size_t i) override;
};

#endif