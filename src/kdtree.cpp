#include <numeric>
#include <iostream>
#include <algorithm>
#include "kdtree.hpp"

// sort indices (inplace) according to the positions of particles in the requested dimension
template<std::size_t NDIM>
void argsort(
    std::vector<Particle<NDIM>> &particles,
    std::vector<std::size_t> &indices,
    int dim,
    int start,
    int stop
){
    sort(
        indices.begin()+start, indices.begin()+stop,
        [dim, &particles](std::size_t i, std::size_t j){
            return particles[i].pos[dim] < particles[j].pos[dim];
        }
    );
}


template<std::size_t NDIM>
KDTree<NDIM>::KDTree(std::vector<Particle<NDIM>> &pars): particles(pars) {
    std::size_t n = particles.size();

    // initialise indices array
    indices.resize(n);
    std::iota(indices.begin(), indices.end(), 0);

    // initialise nodes vector
    nodes.resize(n);

    construct_tree();
}


template<std::size_t NDIM>
void KDTree<NDIM>::construct_tree(){
    int node_ptr = 0;

    create_node(node_ptr, 0, indices.size(), 0);

}


template<std::size_t NDIM>
int KDTree<NDIM>::create_node(int &idx, int start, int stop, int dim){

    int n = stop - start;
    int myidx = idx;
    int particle_index = start + n/2;
    int lstop, ustart;
    int newdim = (dim+1) % NDIM;

    // no particles in node
    if (n == 0){
        return -1;
    }

    argsort(particles, indices, dim, start, stop);

    nodes[myidx].dim = dim;
    nodes[myidx].idx = particle_index;
    // nodes[myidx].mass = ...
    // nodes[myidx].centre_of_mass = ...

    idx += 1;

    lstop = start + n/2;
    ustart = start + n/2+1;

    nodes[myidx].lower = create_node(idx, start, lstop, newdim);
    nodes[myidx].upper = create_node(idx, ustart, stop, newdim);

    return myidx;
}


template class KDTree<2>;
template class KDTree<3>;