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
        std::vector<Particle<NDIM>> &particles;
    public:
        Force(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles): kdtree(kdtree), particles(particles){};
        virtual std::array<double, NDIM> get_force(std::size_t i) = 0;
        virtual void full_update(){};
        virtual void partial_update(){};
        virtual double get_interaction_lengthscale(std::size_t) {return INFINITY;}
        virtual double get_energy(std::size_t) = 0;
        virtual double compute_interaction_timescale(std::size_t){return INFINITY;};
};


// Contains an arbitrary number of force objects and loops over them to compute total forces/energies etc
template<std::size_t NDIM>
class ForceContainer: public Force<NDIM>{
    public:
        std::vector<Force<NDIM>*> forces;

    public:
        ForceContainer(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles);
        void register_force(Force<NDIM>*);
        std::array<double, NDIM> get_force(std::size_t i) override;
        void full_update() override;
        void partial_update() override;
        double get_energy(std::size_t) override;
        double compute_interaction_timescale(std::size_t i) override;
};


// Gravitational force as computed by a KDTree
template<std::size_t NDIM>
class KDTreeGravForce: public Force<NDIM>{
    public:
        KDTreeGravForce(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles);
        std::array<double, NDIM> get_force(std::size_t i) override;
        double get_energy(std::size_t i) override;
        double compute_interaction_timescale(std::size_t) override;
        double get_interaction_lengthscale(std::size_t) override;
};


// Constant gravitational force, e.g. f = (cx, cy, cz)
template<std::size_t NDIM>
class ConstantGravForce: public Force<NDIM>{
    private:
        std::array<double, NDIM> f;
    public:
        ConstantGravForce(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles, std::array<double, NDIM> f);
        std::array<double, NDIM> get_force(std::size_t i) override;
        double get_energy(std::size_t i) override;
};

#endif