#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <state.hpp>
#include "kdtree.hpp"
#include "forces.hpp"

// Parameter to control the minimum timestep. Smaller = shorter timestep
const double TIMESTEP_ETA = 0.05;

// Minimum timestep for an integrator =  dt/2**MAXDEPTH
const int MAXDEPTH = 10;

// ratio of Integrator's timestep and the minimum allowed timestep
const int SMALLEST_TIMESTEP_RATIO = pow(2, MAXDEPTH);

template<std::size_t NDIM>
class Integrator{
    protected:
        double dt;
        State<NDIM> &state;
        Force<NDIM> &force;

    public:
        Integrator(State<NDIM> &s, Force<NDIM> &force, double timestep);
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
        EulerIntegrator(State<NDIM> &s, Force<NDIM> &force, double dt);
    protected:
        void _step(double dt) override;
};


template<std::size_t NDIM>
class AdaptiveTimestepEulerIntegrator: public EulerIntegrator<NDIM>{
    public:
        AdaptiveTimestepEulerIntegrator(State<NDIM> &s, Force<NDIM> &force, double dt);
        void step() override;
};


template<std::size_t NDIM>
class VelocityVerletIntegrator: public Integrator<NDIM>{
    private:
        std::vector<Particle<NDIM>> old_particles;
    public:
        VelocityVerletIntegrator(State<NDIM> &s, Force<NDIM> &force, double dt);
    protected:
        void _step(double dt) override;
};


template<std::size_t NDIM>
class AdaptiveTimestepVelocityVerletIntegrator: public VelocityVerletIntegrator<NDIM>{
    public:
        AdaptiveTimestepVelocityVerletIntegrator(State<NDIM> &s, Force<NDIM> &force, double dt);
        void step() override;
};


template<std::size_t NDIM>
class BlockTimestepIntegrator: public Integrator<NDIM>{
    private:
        std::vector<Particle<NDIM>> old_pars;
        std::vector<double> lastt;
        std::vector<int> depths;

    public:
        BlockTimestepIntegrator(State<NDIM> &s, Force<NDIM> &force, double dt);

    private:
        void _step(double dt) override;
        double assign_timestep_to_particle(std::size_t i);
};

#endif