#ifndef INI_PARSER_H
#define INI_PARSER_H

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cctype>

class IniParser {
private:
    std::map<std::string, std::map<std::string, std::string>> sections;
    
    // Helper function to trim whitespace
    std::string trim(const std::string& str) {
        size_t start = str.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        size_t end = str.find_last_not_of(" \t\r\n");
        return str.substr(start, end - start + 1);
    }
    
    // Helper function to remove comments
    std::string removeComments(const std::string& line) {
        size_t commentPos = line.find(';');
        if (commentPos != std::string::npos) {
            return line.substr(0, commentPos);
        }
        return line;
    }

public:
    bool loadFile(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "ini_parser.h:" << __LINE__ << ": Error: Could not open file " << filename << std::endl;
            return false;
        }
        
        std::string currentSection = "";
        std::string line;
        
        while (std::getline(file, line)) {
            line = trim(line);
            if (line.empty() || line[0] == ';') continue;
            
            // Check if it's a section header
            if (line[0] == '[') {
                size_t endBracket = line.find(']');
                if (endBracket != std::string::npos) {
                    currentSection = trim(line.substr(1, endBracket - 1));
                    continue;
                }
            }
            
            // Parse key-value pairs
            size_t equalPos = line.find('=');
            if (equalPos != std::string::npos) {
                std::string key = trim(line.substr(0, equalPos));
                std::string value = trim(removeComments(line.substr(equalPos + 1)));
                
                if (!key.empty() && !value.empty()) {
                    sections[currentSection][key] = value;
                }
            }
        }
        
        file.close();
        return true;
    }
    
    // Get string value
    std::string getString(const std::string& section, const std::string& key) {
        auto sectionIt = sections.find(section);
        if (sectionIt == sections.end()) {
            throw std::runtime_error("Section '" + section + "' not found in INI file");
        }
        
        auto keyIt = sectionIt->second.find(key);
        if (keyIt == sectionIt->second.end()) {
            throw std::runtime_error("Key '" + key + "' not found in section '" + section + "'");
        }
        
        return keyIt->second;
    }
    
    // Get integer value
    int getInt(const std::string& section, const std::string& key) {
        std::string value = getString(section, key);
        
        try {
            return std::stoi(value);
        } catch (const std::exception& e) {
            throw std::runtime_error("Could not parse '" + value + "' as integer for key '" + key + "' in section '" + section + "'");
        }
    }
    
    // Get double value
    double getDouble(const std::string& section, const std::string& key) {
        std::string value = getString(section, key);
        
        try {
            return std::stod(value);
        } catch (const std::exception& e) {
            throw std::runtime_error("Could not parse '" + value + "' as double for key '" + key + "' in section '" + section + "'");
        }
    }
    
    // Get vector of doubles (for comma-separated values)
    std::vector<double> getDoubleVector(const std::string& section, const std::string& key) {
        std::string value = getString(section, key);
        std::vector<double> result;
        
        if (value.empty()) {
            throw std::runtime_error("Empty value for key '" + key + "' in section '" + section + "'");
        }
        
        std::stringstream ss(value);
        std::string item;
        
        while (std::getline(ss, item, ',')) {
            item = trim(item);
            if (!item.empty()) {
                try {
                    result.push_back(std::stod(item));
                } catch (const std::exception& e) {
                    throw std::runtime_error("Could not parse '" + item + "' as double in vector for key '" + key + "' in section '" + section + "'");
                }
            }
        }
        
        if (result.empty()) {
            throw std::runtime_error("No valid values found for key '" + key + "' in section '" + section + "'");
        }
        
        return result;
    }
    
    // Check if section exists
    bool hasSection(const std::string& section) {
        return sections.find(section) != sections.end();
    }
    
    // Check if key exists in section
    bool hasKey(const std::string& section, const std::string& key) {
        auto sectionIt = sections.find(section);
        if (sectionIt == sections.end()) return false;
        return sectionIt->second.find(key) != sectionIt->second.end();
    }
    
    // Print all sections and keys (for debugging)
    void printAll() {
        for (const auto& section : sections) {
            std::cout << "[" << section.first << "]" << std::endl;
            for (const auto& keyValue : section.second) {
                std::cout << "  " << keyValue.first << " = " << keyValue.second << std::endl;
            }
        }
    }
};

#endif // INI_PARSER_H 