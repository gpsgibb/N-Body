#include "io.hpp"

template<std::size_t NDIM>
Writer<NDIM>::Writer(State<NDIM> &s): state(s) {
    n = state.particles.size();
};


template<std::size_t NDIM>
void Writer<NDIM>::write() {

    std::string filename = get_filename();

    H5::H5File file(filename, H5F_ACC_TRUNC);

    create_attribute("NDIM", hsize_t(NDIM), file);
    create_attribute("time", state.t, file);
    create_attribute("iteration", hsize_t(state.it), file);
    create_attribute("n", hsize_t(n), file);

    std::vector<double> buf(n);
    std::vector<std::array<double, NDIM>> vec_buf(n);

    // dataset for positions
    for (std::size_t i=0;i<n;i++){
        vec_buf[i] = state.particles[i].pos;
    }
    create_dataset("pos", vec_buf, file);

    // dataset for velocities
    for (std::size_t i=0;i<n;i++){
        vec_buf[i] = state.particles[i].vel;
    }
    create_dataset("vel", vec_buf, file);

    // dataset for accelerations
    for (std::size_t i=0;i<n;i++){
        vec_buf[i] = state.particles[i].acc;
    }
    create_dataset("acc", vec_buf, file);

    // dataset for potential energy
    for (std::size_t i=0;i<n;i++){
        buf[i] = state.particles[i].mass;
    }
    create_dataset("mass", buf, file);

    // dataset for potential energy
    for (std::size_t i=0;i<n;i++){
        buf[i] = state.particles[i].potential_energy;
    }
    create_dataset("potential_energy", buf, file);

    // dataset for kinetic energy
    for (std::size_t i=0;i<n;i++){
        buf[i] = state.particles[i].kinetic_energy;
    }
    create_dataset("kinetic_energy", buf, file);

    file.close();

    std::cout << "Written state to " << filename << std::endl;

};


template<std::size_t NDIM>
void Writer<NDIM>::create_dataset(std::string dataset_name, std::vector<std::array<double, NDIM>> &buf, H5::H5File file) {
    hsize_t dims[2] = {n, 3};
    H5::DataSpace space(2, dims);
    H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, space);
    dataset.write(buf.data(), H5::PredType::NATIVE_DOUBLE);
}

template<std::size_t NDIM>
void Writer<NDIM>::create_dataset(std::string dataset_name, std::vector<double> &buf, H5::H5File file) {
    hsize_t dims[1] = {n};
    H5::DataSpace space(1, dims);
    H5::DataSet dataset = file.createDataSet(dataset_name, H5::PredType::NATIVE_DOUBLE, space);
    dataset.write(buf.data(), H5::PredType::NATIVE_DOUBLE);
}


template<std::size_t NDIM>
void Writer<NDIM>::create_attribute(std::string key, double value,  H5::H5File file){
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute attribute = file.createAttribute(key, H5::PredType::NATIVE_DOUBLE, scalar);
    attribute.write(H5::PredType::NATIVE_DOUBLE, &value);
}


template<std::size_t NDIM>
void Writer<NDIM>::create_attribute(std::string key, hsize_t value,  H5::H5File file){
    H5::DataSpace scalar(H5S_SCALAR);
    H5::Attribute attribute = file.createAttribute(key, H5::PredType::NATIVE_HSIZE, scalar);
    attribute.write(H5::PredType::NATIVE_HSIZE, &value);
}


template<std::size_t NDIM>
std::string Writer<NDIM>::get_filename(std::string root, int w, std::string ext){
    int i = state.it;
    std::ostringstream oss;
    oss << root + "_" << std::setw(w) << std::setfill('0') << i << ext;
    return oss.str();
}


template class Writer<3>;
template class Writer<2>;


template<std::size_t NDIM>
State<NDIM> Reader<NDIM>::read(std::string filename){
    H5::H5File file(filename, H5F_ACC_RDONLY);

    std::size_t ndim, n;
    int it;
    double t;

    read_attribute("NDIM", ndim, file);
    if (ndim != NDIM) {
        std::cerr << "File has NDIM of " << ndim << " but this object is set up for NDIM = " << NDIM << std::endl;
        std::exit(1);
    }

    read_attribute("time", t, file);
    read_attribute("iteration", it, file);
    read_attribute("n", n, file);

    State<NDIM> state(n);

    state.it = it;
    state.t = t;

    std::vector<double> scalarbuff(n);
    std::vector<std::array<double, NDIM>> vecbuff(n);

    // read positions
    read_dataset("pos", vecbuff, file);
    for (std::size_t i=0;i<n;i++){
        state.particles[i].pos = vecbuff[i];
    }

    // read velocities
    read_dataset("vel", vecbuff, file);
    for (std::size_t i=0;i<n;i++){
        state.particles[i].vel = vecbuff[i];
    }

    // read accelerations
    read_dataset("acc", vecbuff, file);
    for (std::size_t i=0;i<n;i++){
        state.particles[i].acc = vecbuff[i];
    }

    // read mass
    read_dataset("mass", scalarbuff, file);
    for (std::size_t i=0;i<n;i++){
        state.particles[i].mass = scalarbuff[i];
    }

    // read pe
    read_dataset("potential_energy", scalarbuff, file);
    for (std::size_t i=0;i<n;i++){
        state.particles[i].potential_energy = scalarbuff[i];
    }

    // read ke
    read_dataset("kinetic_energy", scalarbuff, file);
    for (std::size_t i=0;i<n;i++){
        state.particles[i].kinetic_energy = scalarbuff[i];
    }

    file.close();

    std::cout << "Read state from file " << filename << std::endl;

    return state;
}


template<std::size_t NDIM>
void Reader<NDIM>::read_attribute(std::string key, int &value, H5::H5File file){
    H5::Attribute attr = file.openAttribute(key);
    hsize_t val;
    attr.read(H5::PredType::NATIVE_HSIZE, &val);
    value = int(val);
}


template<std::size_t NDIM>
void Reader<NDIM>::read_attribute(std::string key, std::size_t &value, H5::H5File file){
    H5::Attribute attr = file.openAttribute(key);
    hsize_t val;
    attr.read(H5::PredType::NATIVE_HSIZE, &val);
    value = std::size_t(val);
}


template<std::size_t NDIM>
void Reader<NDIM>::read_attribute(std::string key, double &value, H5::H5File file){
    H5::Attribute attr = file.openAttribute(key);
    attr.read(H5::PredType::NATIVE_DOUBLE, &value);
}


template<std::size_t NDIM>
void Reader<NDIM>::read_dataset(std::string name, std::vector<double> &buf, H5::H5File file){
    H5::DataSet dataset = file.openDataSet(name);
    dataset.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
}


template<std::size_t NDIM>
void Reader<NDIM>::read_dataset(std::string name, std::vector<std::array<double, NDIM>> &buf, H5::H5File file){
    H5::DataSet dataset = file.openDataSet(name);
    dataset.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
}

template class Reader<3>;
template class Reader<2>;
