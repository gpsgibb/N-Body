#include <state.hpp>
#include <kdtree.hpp>


template<std::size_t NDIM>
class Integrator{
    protected:
        double dt;
        State<NDIM> &state;
        KDTree<NDIM> tree;

    public:
        Integrator(State<NDIM> &s, double timestep): state(s), tree(s.particles) {
            dt = timestep;
        };
        virtual void step() = 0;
};


template<std::size_t NDIM>
class EulerIntegrator: public Integrator<NDIM>{
    public:
        EulerIntegrator(State<NDIM> &s, double dt);
        void step() override;
};
