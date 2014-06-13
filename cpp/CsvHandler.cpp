#include <fstream>
#include <boost/lexical_cast.hpp>
#include "StringTokenize.hpp"
#include "CsvHandler.hpp"

// TODO: Refactor. loadTimeSignal should not duplicate
// functionality in loadFloatVector!

void stringVectorToFloatVector(const std::vector<std::string>& strVec,
                               std::vector<float>& /*out*/ res) {
    res.clear();
    for (auto it = strVec.begin(); it != strVec.end(); ++it) {
        res.push_back(boost::lexical_cast<float>(*it));
    }
}

// TODO: Should ideally test that the file exists rather than
// test length of token vector...
std::vector<float> loadFloatVector(const std::string& filename) {
    std::vector<float> res;
                      
    std::ifstream infile(filename);
    std::string line;
    std::vector<std::string> stringVector;
    
    // Read the times
    std::getline(infile, line);
    stringVector = StringTokenize(line, ", ");
    stringVectorToFloatVector(stringVector, res);
    if (res.size() == 0) {
        throw std::runtime_error("loadFloatVector: empty file: " + filename);
    }
    return res;
}

template <typename T>
void csvWriteSingleVector(std::ofstream& outStream, const std::vector<T>& v) {
    size_t n = v.size();
    for (size_t i = 0; i < n-1; i++) {
        outStream << v[i] << " ";
    }
    outStream << v[n-1] << std::endl;
}

void saveFloatVector(const std::string& filename,
                     const std::vector<float>& v) {
    std::ofstream outStream(filename);
    csvWriteSingleVector(outStream, v);
    outStream.close();    
}                     

template <typename T>
void csvReadSingleVector(std::ifstream& inStream, std::vector<T>& v) {
    std::string line;
    std::getline(inStream, line);
    auto stringVector = StringTokenize(line, ", ");
    stringVectorToFloatVector(stringVector, v);
}

void loadTimeSignal(const std::string& filename,
                    std::vector<float>& /*out*/ times,
                    std::vector<float>& /*out*/ samples) {
    std::ifstream inStream(filename);
    csvReadSingleVector(inStream, times);
    csvReadSingleVector(inStream, samples);
    inStream.close();
}

void storeTimeSignal(const std::string& filename,
                     const std::vector<float>& times,
                     const std::vector<float>& samples) {

    std::ofstream outStream(filename);
    if (times.size() != samples.size()) {
        throw std::runtime_error("Vector length mismatch");
    }                    
    csvWriteSingleVector(outStream, times);
    csvWriteSingleVector(outStream, samples);    
    outStream.close();
}                     
