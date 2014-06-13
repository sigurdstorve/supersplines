// RealtimeCppAV - Developed by Sigurd Storve (sigurd.storve@ntnu.no)
// This file contains primitive facilities for CSV handling.

#pragma once
#include <string>
#include <vector>

// Convert std::vector of std::string to std::vector of float.
void stringVectorToFloatVector(const std::vector<std::string>& strVec,
                               std::vector<float>& /*out*/ res);

// Load single trace
// Throws if errors.
std::vector<float> loadFloatVector(const std::string& filename);

// Save a single trace.
// Throws if errors
void saveFloatVector(const std::string& filename,
                     const std::vector<float>& v);

// Load a trace (times, samples) from a CSV file.
// Throws if errors.
void loadTimeSignal(const std::string& filename,
                    std::vector<float>& /*out*/ times,
                    std::vector<float>& /*out*/ samples);

// Store a trace (times, samples) to a CSV file.
// Throws if errors.
void storeTimeSignal(const std::string& filename,
                     const std::vector<float>& times,
                     const std::vector<float>& samples); 
