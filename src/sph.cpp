#include "sph.hpp"
#include <cmath>


template<std::size_t NDIM>
double SPHKernel<NDIM>::evaluate(double r, double h){
    double q = r / h;
    return norm / pow(h, NDIM) * kernel_function(q);
}

template<std::size_t NDIM>
std::array<double, NDIM> SPHKernel<NDIM>::evaluate_gradiant(std::array<double, NDIM> r, double h){
    double absr = 0.0;
    for (std::size_t k=0;k<NDIM;k++){
        absr += r[k]*r[k];
    }
    absr = sqrt(absr);

    double q = absr/h;

    double fact = norm / pow(h, NDIM+1) * kernel_function_deriv(q) / absr;

    std::array<double, NDIM> deriv;
    for (std::size_t k=0;k<NDIM;k++){
        deriv[k] = fact * r[k];
    }
    return deriv;

}


template<std::size_t NDIM>
CubicKernel<NDIM>::CubicKernel(){
    if (NDIM == 2){
        this->norm = 10.0 / 7.0 / M_PI;
    } else if (NDIM == 3) {
        this -> norm = 1.0 / M_PI;
    }
}

template<std::size_t NDIM>
double CubicKernel<NDIM>::kernel_function(double q){
    if (q > 2){
        return 0.0;
    } else if (q >= 1){
        return 0.25 * pow(2-q, 3);
    } else {
        return 1 - 1.5*q*q + 0.75*q*q*q;
    }
}


template<std::size_t NDIM>
double CubicKernel<NDIM>::kernel_function_deriv(double q){
    if (q > 2){
        return 0.0;
    } else if (q >= 1){
        return -0.75 * pow(2-q, 2);
    } else {
        return -3*q + 2.25*q*q;
    }
}


template class CubicKernel<2>;
template class CubicKernel<3>;


template<std::size_t NDIM>
WendlandKernel<NDIM>::WendlandKernel(){
    if (NDIM == 2){
        this->norm = 7.0 / 4.0 /  M_PI;
    } else if (NDIM == 3) {
        this -> norm = 21.0 / 16.0 /  M_PI;
    }
}


template<std::size_t NDIM>
double WendlandKernel<NDIM>::kernel_function(double q){
    return pow(1 - 0.5*q, 4) * (2*q + 1);
}


template<std::size_t NDIM>
double WendlandKernel<NDIM>::kernel_function_deriv(double q){
    return 4 * pow(1 - 0.5*q, 3) * (2*q + 1) + 2 * pow(1 - 0.5*q, 4);
}


template class WendlandKernel<2>;
template class WendlandKernel<3>;