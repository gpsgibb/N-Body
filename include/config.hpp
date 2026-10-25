#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <string>
#include <memory>
#include <sstream>

class Config{
    private:
        std::map<std::string, std::unique_ptr<std::map<std::string, std::string>>> config_dict;
    
    public:
        Config(std::string filename);
        template<typename T>
        T get(std::string section, std::string key);
        template<typename T>
        T get(std::string section, std::string key, T default_value);

    private:
        void check_section_exists(std::string section);
};



// get a value of specified type from the config (no default)
template<typename T>
T Config::get(std::string section, std::string key){
    check_section_exists(section);
    if (!(*config_dict[section]).count(key)){
        std::cerr << "Requested key " << key << " does not exist in section " << section << std::endl;
        std::exit(1);
    }
    std::string value = (*config_dict[section]).at(key);
    std::istringstream iss(value);
    T output;
    iss >> output;

    return output;
}


// get a value of specified type from the config (with provided default)
template<typename T>
T Config::get(std::string section, std::string key, T default_value){
    check_section_exists(section);
    if ((*config_dict[section]).count(key)) return get<T>(section, key);

    return default_value;
}


#endif