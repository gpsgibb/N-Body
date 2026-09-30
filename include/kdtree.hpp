#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include <array>
#include "particle.hpp"
#include <string>
#include <fstream>
#include <iostream>

typedef std::array<double, 2> tuple;

const double MAC_CONDITION = 0.5;
const double SMOOTHING_SIZE = 1E-3;

template<std::size_t NDIM>
struct TreeNode{
    std::size_t idx;
    int dim;
    int lower, upper;
    double mass;
    std::array<double, NDIM> centre_of_mass;
    double size;
    std::array<tuple, NDIM> bounding_box;
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
        std::array<double, NDIM> evaluate_force(std::size_t ind);
        void save_tree(std::string filename);

    private:
        int create_node(int &index, int start, int stop, int dim, std::array<tuple, NDIM> bounds);
        std::array<double, NDIM> calc_force_from_node(std::size_t node_idx, std::size_t par_idx);
};

# endif