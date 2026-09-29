#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include <array>
#include "particle.hpp"

typedef std::array<double, 2> tuple;

template<std::size_t NDIM>
struct TreeNode{
    std::size_t idx;
    int dim;
    int lower, upper;
    double mass;
    std::array<double, NDIM> centre_of_mass;
    double size;
};


template<std::size_t NDIM>
class KDTree{
    private:
        std::vector<TreeNode<NDIM>> nodes;
        std::vector<Particle<NDIM>> &particles;
        std::vector<std::size_t> indices;

    public:
        KDTree(std::vector<Particle<NDIM>> &pars);
        void construct_tree();

    private:
        int create_node(int &index, int start, int stop, int dim, std::array<tuple, NDIM> bounds);
};

# endif