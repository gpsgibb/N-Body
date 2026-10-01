#include <integrator.hpp>
#include <iostream>
#include <omp.h>


template<std::size_t NDIM>
EulerIntegrator<NDIM>::EulerIntegrator(State<NDIM> &s, double dt): Integrator<NDIM>(s, dt){
    std::cout << "Initialising EulerIntegrator\n";
}

template<std::size_t NDIM>
void EulerIntegrator<NDIM>::step(){

    std::vector<Particle<NDIM>> &particles = this->state.particles;

    this->tree.construct_tree();

    #pragma omp parallel for
    for (std::size_t i=0; i<particles.size();i++){
        particles[i].acc = this->tree.evaluate_force(i);
    }

    // NOTE: OpenMP threading actually makes this loop slower
    for(Particle<NDIM> &par: particles){
        for(std::size_t k=0;k<NDIM;k++){
            par.vel[k] += par.acc[k]*this->dt;
            par.pos[k] += par.vel[k]*this->dt;
        }
    }

    this->state.t += this->dt;

}


template class EulerIntegrator<2>;
template class EulerIntegrator<3>;