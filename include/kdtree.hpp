#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include <array>
#include "particle.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include "config.hpp"

typedef std::array<double, 2> tuple;

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
        double mac_condition;
        double smoothing_size;

    public:
        KDTree(std::vector<Particle<NDIM>> &pars, Config &config);
        void construct_tree();
        std::array<double, NDIM> evaluate_force(std::size_t ind);
        void save_tree(std::string filename);
        double nearest_neighbour(std::size_t i);
        void reorder_particles();
        void update();
        double potential_energy(std::size_t i);
        std::vector<std::size_t> ball_search(std::size_t i, double r, std::size_t nneighbour_guess=50);

    private:
        long create_node(std::size_t index, std::size_t start, std::size_t stop, std::array<tuple, NDIM> bounds, int depth);
        std::array<double, NDIM> calc_force_from_node(std::size_t node_idx, std::size_t par_idx);
        double nearest_neighbour_from_node(std::size_t node_idx, std::size_t par_idx, double mindist);
        void update_node(std::size_t node_idx);
        double potential_energy_from_node(std::size_t node_idx, std::size_t par_idx);
        void ball_search_node(std::size_t node_idx, std::size_t i, double r, std::vector<std::size_t> &inds);
};

# endif