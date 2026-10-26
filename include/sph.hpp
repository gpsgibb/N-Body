#ifndef SPH_H
#define SPH_H

#include <array>

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
};


class TaitEOS: public SPHEquationOfState{
    public:
        TaitEOS(double rho0, double gamma, double cs);
        double get_pressure(double rho) override;
        double get_soundspeed(double rho) override;
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
    private:
        const double rho0;
        const double cs;
};


class IdealGasEOS: public SPHEquationOfState{
    public:
        IdealGasEOS(double rho0, double gamma, double k);
        double get_pressure(double rho) override;
        double get_soundspeed(double rho) override;
    private:
        const double rho0;
        const double gamma;
        const double k;
};


#endif