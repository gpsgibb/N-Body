#include <iostream>
#include "particle.hpp"
#include "state.hpp"
#include "kdtree.hpp"


int main(int argc, char** argv){

    // Particle<2> p2;
    // Particle<3> p3;

    State<2> s2(10);

    std::vector<std::size_t> indices(10);

    std::cout << indices.size() << "\n";

    for (int i=0;i<s2.particles.size();i++){
        s2.particles[i].pos[0] = i;
        s2.particles[i].pos[1] = 10-i;
        indices[i] = i;
    }

    for (int i=0;i<s2.particles.size();i++){
        std::cout << s2.particles[i].pos[0] << " " << s2.particles[i].pos[1] << "\n";
    }
    int dim =0;
    argsort(s2.particles, indices, dim);
    for (int i=0; i<indices.size(); i++){
        std::cout << indices[i] << " ";
    }
    std::cout << "\n";

    dim =1;
    argsort(s2.particles, indices, dim);
    for (int i=0; i<indices.size(); i++){
        std::cout << indices[i] << " ";
    }
    std::cout << "\n";

    // State<2> s2(100);
    // State<3> s3(100);

    // std::cout << s2.particles.size() << "\n";
    // std::cout << s3.particles.size() << "\n";

    // std::cout << p2.pos.size() << "\n";
    // std::cout << p3.pos.size() << "\n";

    // std::cout << sizeof(p2.pos) << "\n";
    // std::cout << sizeof(p3.pos) << "\n";


    // std::cout << p2.pos[0] << p2.pos[1] << "\n";
    // std::cout << p3.pos[0] << p3.pos[2] << p3.pos[2] << "\n";
    return 0;

}