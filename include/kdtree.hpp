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
    long lower, upper;
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
        double nearest_neighbour(std::size_t i);
        void reorder_particles();
        void update_centre_of_mass();

    private:
        int create_node(std::size_t &index, std::size_t start, std::size_t stop, std::array<tuple, NDIM> bounds);
        std::array<double, NDIM> calc_force_from_node(std::size_t node_idx, std::size_t par_idx);
        double nearest_neighbour_from_node(std::size_t node_idx, std::size_t par_idx, double mindist);
        std::array<double, NDIM> update_node_centre_of_mass(std::size_t node_idx);
};

# endif