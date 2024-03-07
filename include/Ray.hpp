#pragma once

#include "Vector.hpp"

class Ray
{
public:
    Ray(Vector origin, Vector direction);
    Vector origin;
    Vector direction;
};