#include <integrator.hpp>
#include <iostream>
#include <cmath>


template<std::size_t NDIM>
Integrator<NDIM>::Integrator(State<NDIM> &s, double timestep): state(s), tree(s.particles) {
    dt = timestep;

    // Compute force on all particles
    // This is because all integrators assume the initial state is pos(t), vel(t), acc(t)
    update_forces_for_all_particles();
};


// Step the system of particles by one timestep
// This can be overridden by child classes that use e.g. variable timesteps
template<std::size_t NDIM>
void Integrator<NDIM>::step(){
    _step(dt);
}


template<std::size_t NDIM>
void Integrator<NDIM>::update_forces_for_all_particles(){
    #pragma omp parallel for
    for (std::size_t i=0;i<state.particles.size();i++){
        state.particles[i].acc = tree.evaluate_force(i);
    }
}


// get the timestep for particle i
template<std::size_t NDIM>
double Integrator<NDIM>::timestep_for_particle(std::size_t i){
    Particle<NDIM> &par = state.particles[i];

    // compute |acceleration| on particle
    // (initialise to small non-zero value to avoid division by 0)
    double acc = 1E-5;
    for (std::size_t k=0; k<NDIM;k++){
        acc += pow(par.acc[k], 2);
    }
    acc = sqrt(acc);

    double rmin = tree.nearest_neighbour(i);

    return TIMESTEP_ETA * sqrt(rmin / acc);
}


// get the minimum timestep across all particles
template<std::size_t NDIM>
double Integrator<NDIM>::get_min_timestep(){
    double min_dt = INFINITY;

    #pragma omp parallel for reduction(min: min_dt)
    for (std::size_t i=0;i<state.particles.size();i++){
        double par_dt = timestep_for_particle(i);
        if (par_dt < min_dt) min_dt = par_dt; 
    }

    return min_dt;
}


template<std::size_t NDIM>
EulerIntegrator<NDIM>::EulerIntegrator(State<NDIM> &s, double dt): Integrator<NDIM>(s, dt){
    std::cout << "Initialising EulerIntegrator\n";
}

template<std::size_t NDIM>
void EulerIntegrator<NDIM>::_step(double dt){

    std::vector<Particle<NDIM>> &particles = this->state.particles;

    // v(t+dt) = v(t) + a(t)*dt
    // x(t+dt) = x(t) + v(t+dt)*dt
    for(Particle<NDIM> &par: particles){
        for(std::size_t k=0;k<NDIM;k++){
            par.vel[k] += par.acc[k]*dt;
            par.pos[k] += par.vel[k]*dt;
        }
    }

    // get a(t+dt)
    this->tree.construct_tree();
    this->update_forces_for_all_particles();

    this->state.t += dt;

}


template class EulerIntegrator<2>;
template class EulerIntegrator<3>;


template<std::size_t NDIM>
VelocityVerletIntegrator<NDIM>::VelocityVerletIntegrator(State<NDIM> &s, double dt): Integrator<NDIM>(s, dt){
    std::cout << "Initialising VelocityVerletIntegrator\n";
}


template<std::size_t NDIM>
void VelocityVerletIntegrator<NDIM>::_step(double dt){

    std::vector<Particle<NDIM>> &particles = this->state.particles;
    // copy old particles as we need this info when computing the updated velocities
    old_particles = particles;

    // x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*dt^2 
    for(Particle<NDIM> &par: particles){
        for(std::size_t k=0;k<NDIM;k++){
            par.pos[k] += par.vel[k]*dt + 0.5*par.acc[k]*pow(dt, 2);
        }
    }

    // update forces (e.g. get a(t+dt))
    this->tree.construct_tree();
    this->update_forces_for_all_particles();

    // v(t+dt) = v(t) + 0.5(a(t) + a(t+dt))
    for (std::size_t i=0;i<particles.size();i++){
        for(std::size_t k=0;k<NDIM;k++){
            particles[i].vel[k] += 0.5*(particles[i].acc[k] + old_particles[i].acc[k])*dt;
        }
    }

    this->state.t += dt;

}

template class VelocityVerletIntegrator<2>;
template class VelocityVerletIntegrator<3>;


template<std::size_t NDIM>
AdaptiveTimestepEulerIntegrator<NDIM>::AdaptiveTimestepEulerIntegrator(
    State<NDIM> &s, double dt
): EulerIntegrator<NDIM>(s, dt)
{
    std::cout << "Initialising AdaptiveTimestepEulerIntegrator" << std::endl;
}


template<std::size_t NDIM>
void AdaptiveTimestepEulerIntegrator<NDIM>::step(){
    double stoptime = this->state.t + this->dt;
    double dt_local;
    int n_steps = 0;
    while (this->state.t < stoptime){
        dt_local = this->get_min_timestep();

        if (dt_local < this->dt/SMALLEST_TIMESTEP_RATIO){
            dt_local = this->dt/SMALLEST_TIMESTEP_RATIO;
            std::cout << "WARNING: minimum timestep reached!" << std::endl;
        }

        if (dt_local > stoptime - this->state.t){
            dt_local = stoptime - this->state.t;
        }

        // advance time forward by dt_local
        this->_step(dt_local);

        n_steps += 1;

    }
    std::cout << "Number of sub-timesteps = " << n_steps << std::endl;
}


template class AdaptiveTimestepEulerIntegrator<2>;
template class AdaptiveTimestepEulerIntegrator<3>;




template<std::size_t NDIM>
AdaptiveTimestepVelocityVerletIntegrator<NDIM>::AdaptiveTimestepVelocityVerletIntegrator(
    State<NDIM> &s, double dt
): VelocityVerletIntegrator<NDIM>(s, dt)
{
    std::cout << "Initialising AdaptiveTimestepVelocityVerletIntegrator" << std::endl;
}


template<std::size_t NDIM>
void AdaptiveTimestepVelocityVerletIntegrator<NDIM>::step(){
    double stoptime = this->state.t + this->dt;
    double dt_local;
    int n_steps = 0;
    while (this->state.t < stoptime){
        dt_local = this->get_min_timestep();

        if (dt_local < this->dt/SMALLEST_TIMESTEP_RATIO){
            dt_local = this->dt/SMALLEST_TIMESTEP_RATIO;
            std::cout << "WARNING: minimum timestep reached!" << std::endl;
        }

        if (dt_local > stoptime - this->state.t){
            dt_local = stoptime - this->state.t;
        }

        // advance time forward by dt_local
        this->_step(dt_local);

        n_steps += 1;

    }
    std::cout << "Number of sub-timesteps = " << n_steps << std::endl;
}


template class AdaptiveTimestepVelocityVerletIntegrator<2>;
template class AdaptiveTimestepVelocityVerletIntegrator<3>;