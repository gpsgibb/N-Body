#ifndef IO_H
#define IO_H

#include <H5Cpp.h>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include "state.hpp"

template<std::size_t NDIM>
class Writer{
    private:
        State<NDIM> &state;
        std::size_t n;
    public:
        Writer(State<NDIM> &state);
        void write();
    private:
        void create_dataset(std::string dataset_name, std::vector<double> &buf,  H5::H5File file);
        void create_dataset(std::string dataset_name, std::vector<std::array<double, NDIM>> &buf, H5::H5File file);
        void create_attribute(std::string key, double value,  H5::H5File file);
        void create_attribute(std::string key, hsize_t value,  H5::H5File file);
        std::string get_filename(std::string root="state", int w=5, std::string ext=".h5");
};


template<std::size_t NDIM>
class Reader{
    public:
        State<NDIM> read(std::string filename);
    private:
        void read_attribute(std::string key, int &value, H5::H5File file);
        void read_attribute(std::string key, std::size_t &value, H5::H5File file);
        void read_attribute(std::string key, double &value, H5::H5File file);
        void read_dataset(std::string name, std::vector<double> &buf, H5::H5File file);
        void read_dataset(std::string name, std::vector<std::array<double, NDIM>> &buf, H5::H5File file);
};

#endif