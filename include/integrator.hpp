#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <state.hpp>
#include <kdtree.hpp>


template<std::size_t NDIM>
class Integrator{
    protected:
        double dt;
        State<NDIM> &state;
        KDTree<NDIM> tree;

    public:
        Integrator(State<NDIM> &s, double timestep);
        virtual void step() = 0;
};


template<std::size_t NDIM>
class EulerIntegrator: public Integrator<NDIM>{
    public:
        EulerIntegrator(State<NDIM> &s, double dt);
        void step() override;
};


template<std::size_t NDIM>
class VelocityVerletIntegrator: public Integrator<NDIM>{
    private:
        std::vector<Particle<NDIM>> old_particles;
    public:
        VelocityVerletIntegrator(State<NDIM> &s, double dt);
        void step() override;
};

#endif