#pragma once

#include "Vector.hpp"
#include "Ray.hpp"

class Object
{
public:
    virtual double intersect(Ray ray, Vector &N) = 0;
};