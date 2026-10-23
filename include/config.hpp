#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <string>
#include <memory>

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


#endif