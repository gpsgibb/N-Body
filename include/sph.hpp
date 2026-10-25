#ifndef SPH_H
#define SPH_H

#include <array>


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




#endif