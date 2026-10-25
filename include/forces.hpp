#ifndef FORCES_H
#define FORCES_H

#include "particle.hpp"
#include "kdtree.hpp"
#include <cmath>
#include <array>

template<std::size_t NDIM>
class Force{
    protected:
        KDTree<NDIM> &kdtree;
    public:
        Force(KDTree<NDIM> &kdtree): kdtree(kdtree) {};
        virtual std::array<double, NDIM> get_force(std::size_t i);
        virtual void full_update(){};
        virtual void partial_update(){};
        virtual double get_interaction_lengthscale(std::size_t) {return INFINITY;}
        virtual double get_energy(std::size_t) = 0;
};


template<std::size_t NDIM>
class ForceContainer: public Force<NDIM>{
    public:
        std::vector<Force<NDIM>*> forces;

    public:
        ForceContainer(KDTree<NDIM> &kdtree);
        void register_force(Force<NDIM>*);
        std::array<double, NDIM> get_force(std::size_t i) override;
        void full_update() override;
        void partial_update() override;
        double get_interaction_lengthscale(std::size_t) override;
        double get_energy(std::size_t) override;
};


template<std::size_t NDIM>
class KDTreeGravForce: public Force<NDIM>{
    public:
        KDTreeGravForce(KDTree<NDIM> &kdtree);
        std::array<double, NDIM> get_force(std::size_t i) override;
        double get_energy(std::size_t i) override;
};

#endif