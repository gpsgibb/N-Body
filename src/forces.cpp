#include "forces.hpp"


template<std::size_t NDIM>
ForceContainer<NDIM>::ForceContainer(KDTree<NDIM> &kdtree): Force<NDIM>(kdtree){};


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
void ForceContainer<NDIM>::partial_update(){
    this->kdtree.update();
    for (Force<NDIM> *force: forces){
        force->partial_update();
    }
}


template<std::size_t NDIM>
double ForceContainer<NDIM>::get_interaction_lengthscale(std::size_t i){
    return this->kdtree.nearest_neighbour(i);
}


template<std::size_t NDIM>
double ForceContainer<NDIM>::get_energy(std::size_t i){
    double e = 0.0;
    for (Force<NDIM> *force: forces){
        e += force->get_energy(i);
    }
    return e;
}


template class ForceContainer<2>;
template class ForceContainer<3>;


template<std::size_t NDIM>
KDTreeGravForce<NDIM>::KDTreeGravForce(KDTree<NDIM> &kdtree): Force<NDIM>(kdtree){}


template<std::size_t NDIM>
std::array<double, NDIM> KDTreeGravForce<NDIM>::get_force(std::size_t i){
    return this->kdtree.evaluate_force(i);
}


template<std::size_t NDIM>
double KDTreeGravForce<NDIM>::get_energy(std::size_t i){
    return this->kdtree.potential_energy(i);
}


template class KDTreeGravForce<2>;
template class KDTreeGravForce<3>;


template<std::size_t NDIM>
ConstantGravForce<NDIM>::ConstantGravForce(KDTree<NDIM> &kdtree, std::vector<Particle<NDIM>> &particles, std::array<double, NDIM> f): Force<NDIM>(kdtree), f(f), particles(particles) {}


template<std::size_t NDIM>
std::array<double, NDIM> ConstantGravForce<NDIM>::get_force(std::size_t i){
    (void) i;
    return f;
}


template<std::size_t NDIM>
double ConstantGravForce<NDIM>::get_energy(std::size_t i){
    double e = 0.0;
    for (std::size_t k=0;k<NDIM;k++){
        e += particles[i].pos[k] * f[k];
    }
    return e;
}


template class ConstantGravForce<2>;
template class ConstantGravForce<3>;
