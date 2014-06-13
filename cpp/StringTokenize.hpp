// RealtimeCppAV - Developed by Sigurd Storve (sigurd.storve@ntnu.no)
// This file contains convenience code to split ut a string. 

#include <iostream>
#include <string>
#include <vector>
#include <boost/tokenizer.hpp>

inline std::vector<std::string> StringTokenize(const std::string& str, const std::string& separators) {
    boost::char_separator<char> sep(separators.c_str());
    boost::tokenizer<boost::char_separator<char> > tokens(str, sep);
    std::vector<std::string> res;
    for (const auto& t : tokens) {
        res.push_back(t);
    }    
    return res;
}
