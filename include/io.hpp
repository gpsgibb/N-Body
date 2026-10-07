
#include <H5Cpp.h>
#include <string>
#include "state.hpp"

template<std::size_t NDIM>
class Writer{
    private:
        State<NDIM> &state;
        std::size_t n;
        std::vector<std::array<double, NDIM>> vec;
        std::vector<double> scalar;
    public:
        Writer(State<NDIM> &state);
        void write(std::string filename);
    private:
        void create_vector_dataset(std::string dataset_name, H5::H5File file);
        void create_scalar_dataset(std::string dataset_name, H5::H5File file);
        void create_attribute(std::string key, double value,  H5::H5File file);
        void create_attribute(std::string key, hsize_t value,  H5::H5File file);
};