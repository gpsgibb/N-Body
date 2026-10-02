#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <state.hpp>
#include <kdtree.hpp>

const double TIMESTEP_ETA = 0.05;

template<std::size_t NDIM>
class Integrator{
    protected:
        double dt;
        State<NDIM> &state;
        KDTree<NDIM> tree;

    public:
        Integrator(State<NDIM> &s, double timestep);
        virtual void step();

    protected:
        virtual void _step(double dt) = 0;
        void update_forces_for_all_particles();
        double timestep_for_particle(std::size_t i);
        double get_min_timestep();
};


template<std::size_t NDIM>
class EulerIntegrator: public Integrator<NDIM>{
    public:
        EulerIntegrator(State<NDIM> &s, double dt);
    private:
        void _step(double dt) override;
};


template<std::size_t NDIM>
class VelocityVerletIntegrator: public Integrator<NDIM>{
    private:
        std::vector<Particle<NDIM>> old_particles;
    public:
        VelocityVerletIntegrator(State<NDIM> &s, double dt);
    private:
        void _step(double dt) override;
};

#endif