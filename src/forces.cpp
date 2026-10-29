#include <algorithm>
#include "forces.hpp"


template<std::size_t NDIM>
ForceContainer<NDIM>::ForceContainer(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles): Force<NDIM>(kdtree, particles){};


template<std::size_t NDIM>
void ForceContainer<NDIM>::register_force(Force<NDIM> *force){
    forces.push_back(force);
}


template<std::size_t NDIM>
std::array<double, NDIM> ForceContainer<NDIM>::get_force(std::size_t i){
    std::array<double, NDIM> f{}, fp;
    for (Force<NDIM> *force: forces){
        fp = force->get_force(i);
        for (std::size_t k=0;k<NDIM;k++){
            f[k] += fp[k];
        }
    }
    return f;
}


template<std::size_t NDIM>
void ForceContainer<NDIM>::full_update(){
    this->kdtree.construct_tree();
    for (Force<NDIM> *force: forces){
        force->full_update();
    }
}


template<std::size_t NDIM>
void ForceContainer<NDIM>::partial_update(std::vector<std::size_t> &inds){
    this->kdtree.update();
    for (Force<NDIM> *force: forces){
        force->partial_update(inds);
    }
}


template<std::size_t NDIM>
double ForceContainer<NDIM>::get_energy(std::size_t i){
    double e = 0.0;
    for (Force<NDIM> *force: forces){
        e += force->get_energy(i);
    }
    return e;
}


template<std::size_t NDIM>
double ForceContainer<NDIM>::compute_interaction_timescale(std::size_t i){
    std::vector<double> timesteps;

    for (Force<NDIM> *force: forces){
        timesteps.push_back(force->compute_interaction_timescale(i));
    }

    auto result = std::min_element(timesteps.begin(), timesteps.end());

    return *result;
}


template class ForceContainer<2>;
template class ForceContainer<3>;


template<std::size_t NDIM>
KDTreeGravForce<NDIM>::KDTreeGravForce(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles): Force<NDIM>(kdtree, particles){}


template<std::size_t NDIM>
std::array<double, NDIM> KDTreeGravForce<NDIM>::get_force(std::size_t i){
    return this->kdtree.evaluate_force(i);
}


template<std::size_t NDIM>
double KDTreeGravForce<NDIM>::get_energy(std::size_t i){
    return this->kdtree.potential_energy(i);
}


template<std::size_t NDIM>
double KDTreeGravForce<NDIM>::get_interaction_lengthscale(std::size_t i){
    return this->kdtree.nearest_neighbour(i);
}


template<std::size_t NDIM>
double KDTreeGravForce<NDIM>::compute_interaction_timescale(std::size_t i){
    // r_interact = |a| * t^2 -> t = sqrt(r_interact / |a|)
    double r_interact = this->kdtree.nearest_neighbour(i);

    double acc = 1E-10;
    for (std::size_t k=0;k<NDIM;k++){
        acc += pow(this->particles[i].acc[k], 2);
    }

    acc = sqrt(acc);

    if (acc == 0) return INFINITY;

    return sqrt(r_interact / acc);

}


template class KDTreeGravForce<2>;
template class KDTreeGravForce<3>;


template<std::size_t NDIM>
ConstantGravForce<NDIM>::ConstantGravForce(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles, std::array<double, NDIM> f): Force<NDIM>(kdtree, particles), f(f){}


template<std::size_t NDIM>
std::array<double, NDIM> ConstantGravForce<NDIM>::get_force(std::size_t i){
    (void) i;
    return f;
}


template<std::size_t NDIM>
double ConstantGravForce<NDIM>::get_energy(std::size_t i){
    double e = 0.0;
    for (std::size_t k=0;k<NDIM;k++){
        e += this->particles[i].pos[k] * f[k];
    }
    return e;
}


template class ConstantGravForce<2>;
template class ConstantGravForce<3>;
