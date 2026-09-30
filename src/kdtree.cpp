#include <numeric>
#include <iostream>
#include <algorithm>
#include "kdtree.hpp"
#include <cmath>
#include <string>

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


// checks if particle is contained by a bounding box
template<std::size_t NDIM>
bool contained_by_bounding_box(std::array<tuple, NDIM> box, std::array<double, NDIM> pos){
    bool is_contained = 1;
    for (std::size_t i=0;i<NDIM;i++){
        is_contained &= (pos[i] >= box[i][0]) & (pos[i] <= box[i][1]);
    }

    return is_contained;
}


// KDTree initialiser
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


// Construct the tree structure
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
    std::cout << "Centre of mass ";
    for (std::size_t i=0; i<NDIM;i++){
        std::cout << nodes[0].centre_of_mass[i] << " ";
    }
    std::cout << "\n";
}


// Create the node of a tree (recursively creating all child nodes if required)
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

    node.bounding_box = bounds;
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
    upperbounds[dim][0] = particles[indices[lstop]].pos[dim];
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

    node.size = bounding_box_size<NDIM>(bounds);

    return myidx;
}


// Calculate the force on the ith particle from all other particles using the KDTree
template<std::size_t NDIM>
std::array<double, NDIM> KDTree<NDIM>::evaluate_force(std::size_t i){
    return calc_force_from_node(0, i);
}


// Calculate the force on particle par_idx from node node_idx. The function is called recursively to evaluate the force
template<std::size_t NDIM>
std::array<double, NDIM> KDTree<NDIM>::calc_force_from_node(std::size_t node_idx, std::size_t par_idx){
    std::array<double, NDIM> force = {}, r;
    TreeNode node = nodes[node_idx];
    Particle par = particles[par_idx];
    double dist=0;
    bool mac, self, contained;

    contained = contained_by_bounding_box(node.bounding_box, par.pos);

    for (std::size_t i=0;i<NDIM;i++){
        r[i] = node.centre_of_mass[i] - par.pos[i];
        dist += r[i]*r[i];
    }
    dist = sqrt(dist);

    mac = (node.size / dist) > MAC_CONDITION;
    self = node.idx == par_idx;

    if (self || mac || contained){
        // go down tree
        std::array<double, NDIM> f1 = {}, f2 = {};
        if (node.lower >= 0) f1 = calc_force_from_node(node.lower, par_idx);
        if (node.upper >= 0) f2 = calc_force_from_node(node.upper, par_idx);

        for (std::size_t i=0;i<NDIM;i++){
            force[i] = f1[i] + f2[i];
        }

        // child nodes do not contain the particle represented by this node, so need to explicitly
        // include it in the force calculation
        if (!self){
            dist = 0;
            for (std::size_t i=0;i<NDIM;i++){
                r[i] = particles[node.idx].pos[i] - par.pos[i];
                dist += r[i]*r[i];
            }
            dist = sqrt(dist);
            double distcubed = pow(dist + SMOOTHING_SIZE, 3);
            for (std::size_t i=0;i<NDIM;i++){
                force[i] += particles[node.idx].mass* r[i] / distcubed;
            }
        }
    } else {
        // evaluate force from this node and exit
        // f = M r / (|r|+a)**3
        double distcubed = pow(dist + SMOOTHING_SIZE, 3);
        for (std::size_t i=0;i<NDIM;i++){
            force[i] = node.mass * r[i] / distcubed;
        }
    }

    return force;

}


// write the tree to file
template<std::size_t NDIM>
void KDTree<NDIM>::save_tree(std::string filename){
    std::ofstream out(filename, std::ios::binary);
    if (!out) {
        std::cerr << "Failed to open file\n";
        throw 1;
    }

    std::size_t ndim = NDIM, n=nodes.size(), size = sizeof(TreeNode<NDIM>);

    out.write((char*) &ndim, sizeof(ndim));
    out.write((char*) &n, sizeof(std::size_t));
    out.write((char*) &size, sizeof(std::size_t));
    out.write((char*) nodes.data(), nodes.size() * size);

}


template class KDTree<2>;
template class KDTree<3>;