#include "forces.hpp"


template<std::size_t NDIM>
KDTreeGravForce<NDIM>::KDTreeGravForce(std::vector<Particle<NDIM>> &particles, KDTree<NDIM> &tree): particles(particles), tree(tree){}


template<std::size_t NDIM>
std::array<double, NDIM> KDTreeGravForce<NDIM>::get_force(std::size_t i){
    return tree.evaluate_force(i);
}


template<std::size_t NDIM>
void KDTreeGravForce<NDIM>::partial_update(){
    tree.update();
}


template<std::size_t NDIM>
void KDTreeGravForce<NDIM>::full_update(){
    tree.construct_tree();
}


template<std::size_t NDIM>
double KDTreeGravForce<NDIM>::get_interaction_lengthscale(std::size_t i){
    return tree.nearest_neighbour(i);
}

template<std::size_t NDIM>
double KDTreeGravForce<NDIM>::get_energy(std::size_t i){
    return tree.potential_energy(i);
}

template class KDTreeGravForce<2>;
template class KDTreeGravForce<3>;