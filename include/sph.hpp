#ifndef SPH_H
#define SPH_H

#include "forces.hpp"
#include "config.hpp"

// SPH kernels

template<std::size_t NDIM>
class SPHKernel{
    protected:
        double norm;
    public:
        SPHKernel(){};
        double evaluate(double r, double h);
        std::array<double, NDIM> evaluate_gradiant(std::array<double, NDIM> r, double h);
    protected:
        virtual double kernel_function(double q) = 0;
        virtual double kernel_function_deriv(double q) = 0;
};


template<std::size_t NDIM>
class CubicKernel: public SPHKernel<NDIM>{
    public:
        CubicKernel();
    private:
        double kernel_function(double q) override;
        double kernel_function_deriv(double q) override;
};


template<std::size_t NDIM>
class WendlandKernel: public SPHKernel<NDIM>{
    public:
        WendlandKernel();
    private:
        double kernel_function(double q) override;
        double kernel_function_deriv(double q) override;
};


// Equations of state

class SPHEquationOfState{
    public:
        virtual double get_pressure(double rho) = 0;
        virtual double get_soundspeed(double rho) = 0;
        virtual double get_energy(double rho) =0;
};


class TaitEOS: public SPHEquationOfState{
    public:
        TaitEOS(double rho0, double gamma, double cs);
        double get_pressure(double rho) override;
        double get_soundspeed(double rho) override;
        double get_energy(double rho) override;
    private:
        const double rho0;
        const double gamma;
        const double cs;
};


class IsothermalEOS: public SPHEquationOfState{
    public:
        IsothermalEOS(double rho0, double cs);
        double get_pressure(double rho) override;
        double get_soundspeed(double rho) override;
        double get_energy(double rho) override;
    private:
        const double rho0;
        const double cs;
};


class IdealGasEOS: public SPHEquationOfState{
    public:
        IdealGasEOS(double gamma, double k=1.0);
        double get_pressure(double rho) override;
        double get_soundspeed(double rho) override;
        double get_energy(double rho) override;
    private:
        const double gamma;
        const double k;
};


template<std::size_t NDIM>
class SPHForce: public Force<NDIM>{
    protected:
        SPHKernel<NDIM> &kernel;
        SPHEquationOfState &eos;
        double smoothing_length;
        int n_neighbours;
    public:
        SPHForce(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles, SPHKernel<NDIM> &kernel, SPHEquationOfState &eos, Config &config);
        std::array<double, NDIM> get_force(std::size_t i) override;
        void full_update() override;
        void partial_update(std::vector<std::size_t> &inds) override;
        double get_energy(std::size_t) override;
        double compute_interaction_timescale(std::size_t) override;
};


#endif