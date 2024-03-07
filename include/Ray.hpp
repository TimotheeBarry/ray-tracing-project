#pragma once

#include "Vector.hpp"
#include <limits>

class Ray
{
public:
    Ray(Vector origin, Vector direction, double length = std::numeric_limits<double>::infinity());
    Vector origin;
    Vector direction;
    double length;
};