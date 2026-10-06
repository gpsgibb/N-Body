#include <integrator.hpp>
#include <iostream>
#include <cmath>
#include <map>

template<std::size_t NDIM>
Integrator<NDIM>::Integrator(State<NDIM> &s, Force<NDIM> &force, double timestep): state(s), force(force) {
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
        state.particles[i].acc = force.get_force(i);
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

    double rmin = force.get_interaction_lengthscale(i);

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
EulerIntegrator<NDIM>::EulerIntegrator(State<NDIM> &s, Force<NDIM> &f, double dt): Integrator<NDIM>(s, f, dt){
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
    this->force.full_update();
    this->update_forces_for_all_particles();

    this->state.t += dt;

}


template class EulerIntegrator<2>;
template class EulerIntegrator<3>;


template<std::size_t NDIM>
VelocityVerletIntegrator<NDIM>::VelocityVerletIntegrator(State<NDIM> &s, Force<NDIM> &f, double dt): Integrator<NDIM>(s, f, dt){
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
    this->force.full_update();
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
    State<NDIM> &s, Force<NDIM> &f, double dt
): EulerIntegrator<NDIM>(s, f, dt)
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
    State<NDIM> &s, Force<NDIM> &f, double dt
): VelocityVerletIntegrator<NDIM>(s, f, dt)
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


template<std::size_t NDIM>
BlockTimestepIntegrator<NDIM>::BlockTimestepIntegrator(State<NDIM> &s, Force<NDIM> &f, double dt): Integrator<NDIM>(s, f, dt){

    // copy of the particles. This is needed as we need to store the old positions/velocities/forces of the particles
    old_pars = this->state.particles;

    // the last time each particle was updated
    lastt = std::vector<double>(old_pars.size(), this->state.t);

    /// for storing the depths of each particle
    depths.resize(old_pars.size());
}


// assign timestep block to particle i
template<std::size_t NDIM>
double BlockTimestepIntegrator<NDIM>::assign_timestep_to_particle(std::size_t i){
    int n;
    double dt_i = this->timestep_for_particle(i);
    n = std::max(0, int(ceil(log2(this->dt/dt_i))));
    if (n > MAXDEPTH){
        //std::cerr << "Warning timestep for particle smaller than the minimum allowed timestep" << std::endl;
        n = MAXDEPTH;
    }
    depths[i] = n;
    return dt_i;
}


template<std::size_t NDIM>
void BlockTimestepIntegrator<NDIM>::_step(double dt){
    int maxdepth=0;
    // copy of the particles. This is needed as we need to store the old positions/velocities/forces of the particles
    old_pars = this->state.particles;

    // the last time each particle was updated
    lastt = std::vector<double>(old_pars.size(), this->state.t);

    std::size_t n = old_pars.size();

    // assign timesteps for all particles
    double dt_i;
    #pragma omp parallel for reduction(max: maxdepth)
    for(std::size_t i=0;i<n;i++){
        dt_i = assign_timestep_to_particle(i);
        if (depths[i] > maxdepth) maxdepth=depths[i];
    }

    int step, nsteps;
    nsteps = pow(2, maxdepth);

    double t_now = this->state.t;
    double step_dt = dt / nsteps;

    step = 1;

    while (step <= nsteps){
        t_now += step_dt;

        // predict positions of all particles to t_now
        for (std::size_t i=0;i<n;i++){
            // time since particle was last updated
            dt_i = t_now - lastt[i];
            // x(t+dt) = x(t) + v(t)*dt + 0.5*a(t)*t**2
            for (std::size_t k=0;k<NDIM;k++){
                this->state.particles[i].pos[k] = old_pars[i].pos[k] + old_pars[i].vel[k]*dt_i + 0.5*old_pars[i].acc[k]*pow(dt_i, 2);
            }
        }

        // Do we want this?
        this -> force.partial_update();

        // get indices of active particles
        std::vector<std::size_t> inds;
        inds.reserve(n);
        for (std::size_t i=0;i<n;i++){
            bool active = (step % int(pow(2, maxdepth-depths[i]))) == 0;
            if (active) inds.push_back(i);
        }

        // Now identify particles to be updated on this step
        int maxdepthnew = maxdepth;
        #pragma omp parallel for reduction(max:maxdepthnew) private(dt_i)
        for (std::size_t ind=0;ind<inds.size();ind++){
            int i = inds[ind];
            int depth = depths[i];

            // time since particle was last updated
            dt_i = t_now - lastt[i];

            // get force for particle
            this->state.particles[i].acc = this->force.get_force(i);

            // update velocity v(t+dt) = v(t) + 0.5*(a(t) + a(t+dt))*t**2
            for (std::size_t k=0;k<NDIM;k++){
                this->state.particles[i].vel[k] = old_pars[i].vel[k] + 0.5*(old_pars[i].acc[k] + this->state.particles[i].acc[k])*dt_i;
            }

            // copy positon, velocity and acceleration back to old_pars
            old_pars[i] = this->state.particles[i];
            lastt[i] = t_now;

            // Update timestep of particle
            assign_timestep_to_particle(i);
            // We don't want to decrease the depth of the particle
            if (depth < depths[i]) depths[i] = depth;
            if (maxdepth < depths[i]) maxdepthnew = depths[i];

        }

        // update depth, step etc...
        if (maxdepthnew > maxdepth){
            nsteps *= pow(2, maxdepthnew - maxdepth);
            step_dt /= pow(2, maxdepthnew - maxdepth);
            step *= pow(2, maxdepthnew - maxdepth);
            maxdepth = maxdepthnew;
        }

        step += 1;
    }

    std::cout << "Number of minor timesteps " << nsteps << std::endl;

    std::map<int, int> stepdict;
    for (int d=0;d<=maxdepth;d++){
        stepdict[d] = 0;
    }
    for (std::size_t i=0;i<n;i++){
        stepdict[depths[i]] += 1;
    }
    for (int d=0;d<=maxdepth;d++){
        std::cout << "Depth " << d << " Count " << stepdict[d] << std::endl;
    }
   
    this -> state.t += dt;
    this -> force.full_update();
 
}


template class BlockTimestepIntegrator<2>;
template class BlockTimestepIntegrator<3>;