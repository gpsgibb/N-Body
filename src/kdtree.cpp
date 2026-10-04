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
    std::size_t start,
    std::size_t stop
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

    // initialise nodes vector
    nodes.resize(n);

    construct_tree();

    // sanity checking
    n = 0;
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


// Construct the tree structure
template<std::size_t NDIM>
void KDTree<NDIM>::construct_tree(){
    std::size_t node_idx = 0;

    std::iota(indices.begin(), indices.end(), 0);

    std::array<tuple, NDIM> bounds = minmax(particles);

    // Create root node (which recursively creates all other nodes)
    create_node(node_idx, 0, indices.size(), bounds);

    // Each node's centre of mass currently contains sum(pos*m), so divide sum(pos*m) by sum(m) to
    // get the centre of mass
    for (TreeNode<NDIM> &node: nodes){
        for (std::size_t i=0; i<NDIM;i++){
            node.centre_of_mass[i] /= node.mass;
        }
    }

}


// Create the node of a tree (recursively creating all child nodes if required)
template<std::size_t NDIM>
int KDTree<NDIM>::create_node(std::size_t &idx, std::size_t start, std::size_t stop, std::array<tuple, NDIM> bounds){

    std::size_t n = stop - start;

    // no particles in node
    if (n == 0){
        return -1;
    }

    TreeNode<NDIM> &node = nodes.at(idx);

    // determine dimension to split (the one with the largest bounding box)
    int dim = 0;
    double boxl = bounds[0][1] - bounds[0][0];
    double l;
    for (std::size_t k=1;k<NDIM;k++){
        l = bounds[k][1] - bounds[k][0];
        if (boxl < l){
            boxl = l;
            dim = k;
        }
    }

    // Sort the indices in the dimension of this node
    argsort(particles, indices, dim, start, stop);

    std::size_t myidx = idx;
    std::size_t particle_index = indices[start + n/2];

    node.bounding_box = bounds;
    node.dim = dim;
    node.idx = particle_index;

    // initialise mass and centre of mass to be the mass of the object and the mass * pos
    Particle mypar = particles[particle_index];
    node.mass = mypar.mass;
    for (std::size_t i=0; i<NDIM;i++){
        node.centre_of_mass[i] = mypar.pos[i] * mypar.mass;
    }

    std::size_t lstop, ustart;
    lstop = start + n/2;
    ustart = start + n/2+1;

    // increment the index of the "next" node and create the child nodes
    idx += 1;

    std::array<tuple, NDIM> lowerbounds=bounds, upperbounds=bounds;
    lowerbounds = bounds;
    lowerbounds[dim][1] = particles[indices[lstop]].pos[dim];
    node.lower = create_node(idx, start, lstop, lowerbounds);
    
    upperbounds = bounds;
    upperbounds[dim][0] = particles[indices[ustart]].pos[dim];
    upperbounds[dim][0] = particles[indices[lstop]].pos[dim];
    node.upper = create_node(idx, ustart, stop, upperbounds);

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


// Returns the distance to particle i's nearest neighbour
template<std::size_t NDIM>
double KDTree<NDIM>::nearest_neighbour(std::size_t i){
    double mindist = INFINITY;

    // NOTE: This returns dist^2, so need to take the sqrt
    mindist = nearest_neighbour_from_node(0, i, mindist);
    mindist = sqrt(mindist);

    return mindist;
}


// returns the nearest neighbour distance (squared) to particle par_idx from node node_idx
template<std::size_t NDIM>
double KDTree<NDIM>::nearest_neighbour_from_node(std::size_t node_idx, std::size_t par_idx, double mindist){
    TreeNode node = nodes[node_idx];
    Particle par = particles[par_idx];
    Particle nodepar = particles[node.idx];
    int dim = node.dim;

    // distance between target particle and the particle associated with this node
    double dist = 0;
    for (std::size_t i=0;i<NDIM;i++){
        dist += pow(par.pos[i] - nodepar.pos[i], 2);
    }

    // if this node is closer than the current mindist, update mindist
    // UNLESS this node contains the target particle, in which case we ignore
    if (par_idx != node.idx){
        if (dist < mindist) mindist = dist;
    }

    // determine which of the child nodes the particle is in
    long node_containing, node_adjacent;
    if (par.pos[dim] <= nodepar.pos[dim]){
        node_containing = node.lower;
        node_adjacent = node.upper;
    } else {
        node_containing = node.upper;
        node_adjacent = node.lower;
    }

    // search in containing node
    if (node_containing > -1) mindist = nearest_neighbour_from_node(node_containing, par_idx, mindist);

    // search in adjacent node if it is closer to the particle than the current mindist
    double plane_dist = pow(par.pos[dim] - nodepar.pos[dim], 2);
    if ((node_adjacent > -1) & (plane_dist < mindist)){
        mindist = nearest_neighbour_from_node(node_adjacent, par_idx, mindist);
    }

    return mindist;
}


template<std::size_t NDIM>
void KDTree<NDIM>::reorder_particles(){

    std::size_t n = particles.size();

    // get indices of particles as they appear in the tree
    std::vector<std::size_t> tree_indices;
    tree_indices.reserve(n);
    for (TreeNode<NDIM> node: nodes){
        tree_indices.push_back(node.idx);
    }

    // reorder particles to be in that order, and update the nodes' particle indidces to match
    std::vector<Particle<NDIM>> particles_copy = particles;
    for (std::size_t i=0;i<n;i++){
        particles[i] = particles_copy[tree_indices[i]];
        nodes[i].idx = i;
    }

}


template<std::size_t NDIM>
void KDTree<NDIM>::update_centre_of_mass(){
    // this puts sum(m*x) into the centre of mass fields
    update_node_centre_of_mass(0);

    // divide by mass of node to get centre of mass
    for (TreeNode<NDIM> &node: nodes){
        for (std::size_t k=0;k<NDIM;k++){
            node.centre_of_mass[k] /= node.mass;
        }
    }
}


// return centre of mass (* total mass) of node idx
template<std::size_t NDIM>
std::array<double, NDIM> KDTree<NDIM>::update_node_centre_of_mass(std::size_t idx){
    TreeNode<NDIM> &node = nodes[idx];

    // compute own node particle's centre of mass
    Particle<NDIM> par = particles[node.idx];
    std::array<double, NDIM> mycom;
    for (std::size_t k=0;k<NDIM;k++){
        mycom[k] = par.pos[k] * par.mass;
    }

    //get centre of masses of children and sum
    std::array<double, NDIM> childcom;

    if (node.lower >= 0){
        childcom = update_node_centre_of_mass(node.lower);
        for (std::size_t k=0;k<NDIM;k++){
            mycom[k] += childcom[k];
        }
    }

    if (node.upper >= 0){
        childcom = update_node_centre_of_mass(node.upper);
        for (std::size_t k=0;k<NDIM;k++){
            mycom[k] += childcom[k];
        }
    }

    node.centre_of_mass = mycom;
    return mycom;
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