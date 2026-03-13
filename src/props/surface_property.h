#pragma once
#include <any>
#include <string>
#include <unordered_map>

class Surface;

struct ISurfaceProperty {
    virtual ~ISurfaceProperty() = default;
    virtual std::string key() const = 0;
    virtual void compute(const Surface& s, std::unordered_map<std::string, std::any>& cache) = 0;
};