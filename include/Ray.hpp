#pragma once

#include "Vector.hpp"
#include <limits>

class Ray
{
public:
    Ray(Vector origin, Vector direction, double length = std::numeric_limits<double>::infinity());
    ~Ray() {}
    Vector origin;
    Vector direction;
    double length;
};