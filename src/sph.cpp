#include "sph.hpp"
#include <cmath>

// kernels

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
    if (q > 2) return 0.0;
    return pow(1 - 0.5*q, 4) * (2*q + 1);
}


template<std::size_t NDIM>
double WendlandKernel<NDIM>::kernel_function_deriv(double q){
    if (q > 2) return 0.0;
    return -2 * pow(1 - 0.5*q, 3) * (2*q + 1) + 2 * pow(1 - 0.5*q, 4);
}


template class WendlandKernel<2>;
template class WendlandKernel<3>;


// Equations of state

TaitEOS::TaitEOS(double rho0, double gamma, double cs): rho0(rho0), gamma(gamma), cs(cs){}
double TaitEOS::get_pressure(double rho){
    return rho0*cs*cs/gamma * (pow(rho/rho0, gamma) - 1);
}
double TaitEOS::get_soundspeed(double rho){
    return cs*pow(rho/rho0, (gamma-1)/2);
}
double TaitEOS::get_energy(double rho){
    return rho0*cs*cs/gamma/(gamma-1) * (pow(rho/rho0, gamma-1) -1);
}


IsothermalEOS::IsothermalEOS(double rho0, double cs): rho0(rho0), cs(cs){}
double IsothermalEOS::get_pressure(double rho){
    return cs*cs*(rho-rho0);
}
double IsothermalEOS::get_soundspeed(double rho){
    (void) rho;
    return cs;
}
double IsothermalEOS::get_energy(double rho){
    // system is isothermal so no energy
    (void) rho;
    return 0;
}


IdealGasEOS::IdealGasEOS(double gamma, double k): gamma(gamma), k(k){}
double IdealGasEOS::get_pressure(double rho){
    return k * pow(rho, gamma);
}
double IdealGasEOS::get_soundspeed(double rho){
    return sqrt(gamma * pow(rho, gamma-1));
}
double IdealGasEOS::get_energy(double rho){
    double pressure = get_pressure(rho);
    return k*pressure / (1-gamma)/rho;
}


// SPH force

template<std::size_t NDIM>
SPHForce<NDIM>::SPHForce(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles, SPHKernel<NDIM> &kernel, SPHEquationOfState &eos, Config &config): Force<NDIM>(kdtree, particles), kernel(kernel), eos(eos){
    smoothing_length = config.get<double>("sph", "smoothing_length");
    n_neighbours = config.get<int>("sph", "n_neighbours");

    // fix smoothing length
    if (smoothing_length){
        for (Particle<NDIM> &par: this->particles){
            par.smoothing_length = smoothing_length;
        }
    }

    // recompute smoothing length, densities and pressures
    full_update();

}


template<std::size_t NDIM>
std::array<double, NDIM> SPHForce<NDIM>::get_force(std::size_t i){
    std::array<double, NDIM> force{0};
    Particle<NDIM> &par = this->particles[i];

    // TODO: Should we increase the search radius here to minimise mismatches of pairs with different smoothing lengths?
    std::vector<Neighbour> neighbours = this->kdtree.ball_search(i, 2*par.smoothing_length, n_neighbours);
    for (Neighbour neighbour: neighbours){
        Particle<NDIM> &other = this->particles[neighbour.ind];

        std::array<double, NDIM> rvec{};
        for (std::size_t k=0;k<NDIM;k++){
            rvec[k] = par.pos[k] - other.pos[k];
        }

        // compute gradient with smoothing length set to the mean of the smoothing length between the pair
        std::array<double, NDIM> grad{};
        grad = kernel.evaluate_gradiant(rvec, (par.smoothing_length + other.smoothing_length)/2);

        for (std::size_t k=0;k<NDIM;k++){
            grad[k] *= other.mass;
            grad[k] *= par.pressure / pow(par.density, 2) + other.pressure / pow(other.density, 2);
            // f = -grad(p) so we need a minus sign here
            force[k] -= grad[k];
        }
    }

    return force;
}


template<std::size_t NDIM>
void SPHForce<NDIM>::full_update(){
    this->kdtree.construct_tree();

    // update smoothing length
    if (!smoothing_length){
        #pragma omp parallel for
        for (std::size_t i=0;i<this->particles.size();i++){
            std::vector<Neighbour> neighbours = this->kdtree.k_nearest_neighbours(i, n_neighbours);
            // note: neighbours is sorted so the last element will be the furthest away
            this->particles[i].smoothing_length = neighbours.back().r / 2;
        }
    }

    // compute densities
    #pragma omp parallel for
    for (std::size_t i=0;i<this->particles.size();i++){
        Particle<NDIM> &par = this->particles[i];
        // remember to include self in the density computation
        par.density = par.mass * kernel.evaluate(0, par.smoothing_length);
        std::vector<Neighbour> neighbours = this->kdtree.ball_search(i, 2*par.smoothing_length, n_neighbours);
        for (Neighbour neighbour: neighbours){
            Particle<NDIM> &other = this->particles[neighbour.ind];
            par.density += other.mass * kernel.evaluate(neighbour.r, par.smoothing_length);
        }
    }

    //compute pressure
    #pragma omp parallel for
    for (Particle<NDIM> &par: this->particles){
        par.pressure = eos.get_pressure(par.density);
    }

}


template<std::size_t NDIM>
void SPHForce<NDIM>::partial_update(std::vector<long unsigned int> &inds){

    std::vector<double> old_densities;

    if (!smoothing_length){
        old_densities.reserve(inds.size());
        for (std::size_t i: inds){
            old_densities.push_back(this->particles[i].density);
        }
    }

    // compute densities
    #pragma omp parallel for
    for (std::size_t i: inds){
        Particle<NDIM> &par = this->particles[i];
        // remember to include self in the density computation
        par.density = par.mass * kernel.evaluate(0, par.smoothing_length);
        std::vector<Neighbour> neighbours = this->kdtree.ball_search(i, 2*par.smoothing_length, n_neighbours);
        for (Neighbour neighbour: neighbours){
            Particle<NDIM> &other = this->particles[neighbour.ind];
            par.density += other.mass * kernel.evaluate(neighbour.r, par.smoothing_length);
        }
    }

    // compute update to smoothing_length
    if (!smoothing_length){
        for (std::size_t i=0;i<inds.size();i++){
            this->particles[inds[i]].density *= pow(this->particles[inds[i]].density / old_densities[i], -1.0/NDIM);
        }
    }

    //compute pressure
    #pragma omp parallel for
    for (std::size_t i: inds){
        Particle<NDIM> &par = this->particles[i];
        par.pressure = eos.get_pressure(par.density);
    }


}


template<std::size_t NDIM>
double SPHForce<NDIM>::get_energy(std::size_t i){
    Particle<NDIM> &par = this->particles[i];
    double u = eos.get_energy(par.density);

    // This is the energy per unit mass so multiply by the particle mass
    return u * par.mass;
};


template<std::size_t NDIM>
double SPHForce<NDIM>::compute_interaction_timescale(std::size_t i){
    Particle<NDIM> &par = this->particles[i];
    return par.smoothing_length / eos.get_soundspeed(par.density);
}


template class SPHForce<2>;
template class SPHForce<3>;