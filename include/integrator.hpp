#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "state.hpp"
#include "kdtree.hpp"
#include "forces.hpp"
#include "config.hpp"

template<std::size_t NDIM>
class Integrator{
    protected:
        double dt;
        double eta;
        // NOTE: the following two are only needed in variable timestep subclasses but are included here for simplicity
        int max_allowed_depth;
        int smallest_timestep_ratio;
        State<NDIM> &state;
        Force<NDIM> &force;
        double inactivate_dist;

    public:
        Integrator(State<NDIM> &s, Force<NDIM> &force, Config &config);
        virtual void step();
        void compute_energy();
        void inactivate_particles();

    protected:
        virtual void _step(double dt) = 0;
        void update_forces_for_all_particles();
        double timestep_for_particle(std::size_t i);
        double get_min_timestep();
};


template<std::size_t NDIM>
class EulerIntegrator: public Integrator<NDIM>{
    public:
        EulerIntegrator(State<NDIM> &s, Force<NDIM> &force, Config &config);
    protected:
        void _step(double dt) override;
};


template<std::size_t NDIM>
class AdaptiveTimestepEulerIntegrator: public EulerIntegrator<NDIM>{
    public:
        AdaptiveTimestepEulerIntegrator(State<NDIM> &s, Force<NDIM> &force, Config &config);
        void step() override;
};


template<std::size_t NDIM>
class VelocityVerletIntegrator: public Integrator<NDIM>{
    private:
        std::vector<Particle<NDIM>> old_particles;
    public:
        VelocityVerletIntegrator(State<NDIM> &s, Force<NDIM> &force, Config &config);
    protected:
        void _step(double dt) override;
};


template<std::size_t NDIM>
class AdaptiveTimestepVelocityVerletIntegrator: public VelocityVerletIntegrator<NDIM>{
    public:
        AdaptiveTimestepVelocityVerletIntegrator(State<NDIM> &s, Force<NDIM> &force, Config &config);
        void step() override;
};


template<std::size_t NDIM>
class BlockTimestepIntegrator: public Integrator<NDIM>{
    private:
        std::vector<Particle<NDIM>> old_pars;
        std::vector<double> lastt;
        std::vector<int> depths;

    public:
        BlockTimestepIntegrator(State<NDIM> &s, Force<NDIM> &force, Config &config);

    private:
        void _step(double dt) override;
        double assign_timestep_to_particle(std::size_t i);
};

#endif