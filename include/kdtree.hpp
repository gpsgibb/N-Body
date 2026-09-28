#ifndef KDTREE_H

#define KDTREE_H

#include <vector>
#include "particle.hpp"
#include <algorithm>

// sort indices (inplace) according to the dimension
template<std::size_t NDIM>
void argsort(
    std::vector<Particle<NDIM>> &particles,
    std::vector<std::size_t> &indices,
    int dim
){
    sort(
        indices.begin(), indices.end(),
        [dim, &particles](std::size_t i, std::size_t j){
            return particles[i].pos[dim] < particles[j].pos[dim];
        }
    );
}

template<std::size_t NDIM>
class TreeNode{
    public:
        std::size_t idx;
        int dim;
        TreeNode<NDIM> *lower=NULL, *upper=NULL;
        double mass;

    public:
        TreeNode<NDIM>(
            std::vector<Particle<NDIM>> &particles,
            std::vector<std::size_t> &inds,
            int input_dim
        );
};


template<std::size_t NDIM>
class KDTree{
    private:
        std::vector<TreeNode<NDIM>> nodes;
    public:
        KDTree(std::vector<Particle<NDIM>> &particles);
        void construct(std::vector<Particle<NDIM>> &particles);
};

# endif