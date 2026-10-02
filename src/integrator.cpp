#include <integrator.hpp>
#include <iostream>
#include <omp.h>
#include <cmath>


template<std::size_t NDIM>
Integrator<NDIM>::Integrator(State<NDIM> &s, double timestep): state(s), tree(s.particles) {
    dt = timestep;

    //Compute force on all particles
    // This is because all integrators assume the initial state is pos(t), vel(t), acc(t)
    #pragma omp parallel for
    for(std::size_t i=0;i<state.particles.size();i++){
        state.particles[i].acc = tree.evaluate_force(i);
    }
};


template<std::size_t NDIM>
EulerIntegrator<NDIM>::EulerIntegrator(State<NDIM> &s, double dt): Integrator<NDIM>(s, dt){
    std::cout << "Initialising EulerIntegrator\n";
}

template<std::size_t NDIM>
void EulerIntegrator<NDIM>::step(){

    std::vector<Particle<NDIM>> &particles = this->state.particles;

    // v(t+dt) = v(t) + a(t)*dt
    // x(t+dt) = x(t) + v(t+dt)*dt
    for(Particle<NDIM> &par: particles){
        for(std::size_t k=0;k<NDIM;k++){
            par.vel[k] += par.acc[k]*this->dt;
            par.pos[k] += par.vel[k]*this->dt;
        }
    }

    // get a(t+dt)
    this->tree.construct_tree();
    #pragma omp parallel for
    for (std::size_t i=0; i<particles.size();i++){
        particles[i].acc = this->tree.evaluate_force(i);
    }

    this->state.t += this->dt;

}


template class EulerIntegrator<2>;
template class EulerIntegrator<3>;


template<std::size_t NDIM>
VelocityVerletIntegrator<NDIM>::VelocityVerletIntegrator(State<NDIM> &s, double dt): Integrator<NDIM>(s, dt){
    std::cout << "Initialising VelocityVerletIntegrator\n";
}


template<std::size_t NDIM>
void VelocityVerletIntegrator<NDIM>::step(){

    std::vector<Particle<NDIM>> &particles = this->state.particles;
    // copy old particles as we need this info when computing the updated velocities
    old_particles = particles;

    // x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2 
    for(Particle<NDIM> &par: particles){
        for(std::size_t k=0;k<NDIM;k++){
            par.pos[k] += par.vel[k]*this->dt + 0.5*par.acc[k]*pow(this->dt, 2);
        }
    }

    // update forces (e.g. get a(t+dt))
    this->tree.construct_tree();
    #pragma omp parallel for
    for (std::size_t i=0; i<particles.size();i++){
        particles[i].acc = this->tree.evaluate_force(i);
    }

    // v(t+dt) = v(t) + 0.5(a(t) + a(t+dt))
    for (std::size_t i=0;i<particles.size();i++){
        for(std::size_t k=0;k<NDIM;k++){
            particles[i].vel[k] += 0.5*(particles[i].acc[k] + old_particles[i].acc[k])*this->dt;
        }
    }

    this->state.t += this->dt;

}

template class VelocityVerletIntegrator<2>;
template class VelocityVerletIntegrator<3>;