#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <sstream>

#include "config.hpp"


// load file, returning a list of lines
std::list<std::string> load_file(std::string filename){
    std::ifstream file(filename);
    std::string line;
    std::list<std::string> lines;

    while (std::getline(file, line)){
        lines.push_back(line);
    }

    file.close();

    return lines;
}


void strip_leading_and_trailing_whitespace(std::string &line){
    std::size_t start = line.find_first_not_of(' ');
    std::size_t stop = line.find_last_not_of(' ')+1;

    if (start >= stop){
        start=line.size();
    }

    line.erase(0, start);
    line.erase(stop-start);

}


void strip_whitespace(std::list<std::string> &lines){
    for (std::string &line: lines){
        strip_leading_and_trailing_whitespace(line);
    }
}


void strip_comments(std::list<std::string> &lines){
    for (std::string &line: lines){
        std::size_t pos = line.find("#");
        if (pos != std::string::npos) line.erase(pos);
    }
}


void strip_blanklines(std::list<std::string> &lines){
    lines.erase(std::remove_if(lines.begin(), lines.end(), [](std::string s){return s.size() == 0;}), lines.end());
}


bool is_section_header(std::string s){
    if (s.front() == '[' && s.back() == ']') return true;
    return false;
}


bool is_key_value_pair(std::string s){
    if (s.find('=') != s.npos) return true;
    return false;
}


// read a config file, returning a map of section names to a map of key-value pairs
std::map<std::string, std::unique_ptr<std::map<std::string, std::string>>> read_config(std::string filename){
    std::list<std::string> lines = load_file(filename);

    strip_comments(lines);
    strip_whitespace(lines);
    strip_blanklines(lines);

    std::map<std::string, std::unique_ptr<std::map<std::string, std::string>>> config_dict;

    if (!is_section_header(lines.front())){
        std::cerr << "First entry in file is not a section header: " << lines.front() << std::endl;
        std::exit(1);
    }

    std::string section;
    std::string key, value;
    for (std::string line: lines){
        if (is_section_header(line)){
            section = line.substr(1, line.size()-2);
            config_dict.insert({section, std::make_unique<std::map<std::string, std::string>>()});

        } else if (is_key_value_pair(line)) {
            std::size_t loc = line.find('=');
            key = line.substr(0, loc);
            value = line.substr(loc+1);

            // handle trailing whitespace around the "="
            strip_leading_and_trailing_whitespace(key);
            strip_leading_and_trailing_whitespace(value);

            (*config_dict[section])[key] = value;
            
        } else {
            std::cerr << "Invalid line: " << line << std::endl;
            std::exit(1);
        }
    }

    return config_dict;
}


Config::Config(std::string filename){
    config_dict = read_config(filename);
}


void Config::check_section_exists(std::string section){
    if (config_dict.count(section)) return;
    std::cerr << "Requested section " << section << " is not in the config" << std::endl;
    std::exit(1);
}


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
