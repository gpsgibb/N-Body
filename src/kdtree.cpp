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


// get min/max values of particle positions
template<std::size_t NDIM>
std::array<tuple, NDIM> minmax(std::vector<Particle<NDIM>> &particles){
    std::array<tuple, NDIM> bounds;
    
    // initialise min/max to first entry of particles
    for (std::size_t j=0;j<NDIM;j++){
        for (int i=0;i<2;i++){
            bounds[j][i] = particles[0].pos[j];
        }
    }

    // now loop through all particles and find minmax
    for (Particle par: particles){
        for (std::size_t j=0;j<NDIM;j++){
            if (par.pos[j] < bounds[j][0]) bounds[j][0] = par.pos[j];
            if (par.pos[j] > bounds[j][1]) bounds[j][1] = par.pos[j];
        }
    }

    return bounds;
}


// maximum size of bounding box
template<std::size_t NDIM>
double bounding_box_size(std::array<tuple, NDIM> box){
    double size = 0;
    double dsize;
    for (std::size_t j=0;j<NDIM;j++){
        dsize = box[j][1] - box[j][0];
        if (dsize > size) size = dsize;
    }

    return size;
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
    int node_idx = 0;

    std::array<tuple, NDIM> bounds = minmax(particles);

    // Create root node (which recursively creates all other nodes)
    create_node(node_idx, 0, indices.size(), 0, bounds);

    // Each node's centre of mass currently contains sum(pos*m), so divide sum(pos*m) by sum(m) to
    // get the centre of mass
    for (TreeNode<NDIM> &node: nodes){
        for (std::size_t i=0; i<NDIM;i++){
            node.centre_of_mass[i] /= node.mass;
        }
    }

    // sanity checking
    long n = 0;
    for (TreeNode node: nodes){
        n += node.idx;
    }
    std::cout << "Sum of and theoretical sum of indices ";
    std::cout << n << " " << nodes.size() * (nodes.size()-1) / 2 << "\n";

    std::cout << "Total mass " << nodes[0].mass << "\n";
    std::cout << "Centre of mass " << nodes[0].centre_of_mass[0] << " " << nodes[0].centre_of_mass[1] << " " <<nodes[0].centre_of_mass[2] <<"\n";
}


template<std::size_t NDIM>
int KDTree<NDIM>::create_node(int &idx, int start, int stop, int dim, std::array<tuple, NDIM> bounds){

    int n = stop - start;

    // no particles in node
    if (n == 0){
        return -1;
    }

    TreeNode<NDIM> &node = nodes[idx];

    // Sort the indices in the dimension of this node
    argsort(particles, indices, dim, start, stop);

    int myidx = idx;
    int particle_index = indices[start + n/2];

    node.dim = dim;
    node.idx = particle_index;

    // initialise mass and centre of mass to be the mass of the object and the mass * pos
    Particle mypar = particles[particle_index];
    node.mass = mypar.mass;
    for (std::size_t i=0; i<NDIM;i++){
        node.centre_of_mass[i] = mypar.pos[i] * mypar.mass;
    }

    int lstop, ustart;
    lstop = start + n/2;
    ustart = start + n/2+1;

    // increment the index of the "next" node and create the child nodes
    idx += 1;
    int newdim = (dim+1) % NDIM;

    std::array<tuple, NDIM> lowerbounds=bounds, upperbounds=bounds;
    lowerbounds = bounds;
    lowerbounds[dim][1] = particles[indices[lstop]].pos[dim];
    node.lower = create_node(idx, start, lstop, newdim, lowerbounds);
    
    upperbounds = bounds;
    upperbounds[dim][0] = particles[indices[ustart]].pos[dim];
    node.upper = create_node(idx, ustart, stop, newdim, upperbounds);

    // Add the masses and centre of masses * mass of the child nodes
    if (node.lower > 0){
        TreeNode lower = nodes[node.lower];
        node.mass += lower.mass;
        for (std::size_t i=0; i<NDIM;i++){
            node.centre_of_mass[i] += lower.centre_of_mass[i];
        }
    }
    if (node.upper > 0){
        TreeNode upper = nodes[node.upper];
        node.mass += upper.mass;
        for (std::size_t i=0; i<NDIM;i++){
            node.centre_of_mass[i] += upper.centre_of_mass[i];
        }
    }

    if (n==1){
        node.size = 0.0;
    } else {
        // get size of the node: maximum distance from centre of mass to a corner of its bounding box
        node.size = bounding_box_size<NDIM>(bounds);
    }

    return myidx;
}


template class KDTree<2>;
template class KDTree<3>;