#include "../include/BoundingBox.hpp"
#include "../include/Functions.hpp"
#include <limits>
#include <iostream>
#include <cmath>

bool BoundingBox::intersect(Ray &ray) const
{
    double txMin = (min[0] - ray.origin[0]) / ray.direction[0];
    double txMax = (max[0] - ray.origin[0]) / ray.direction[0];
    if (txMin > txMax)
    {
        std::swap(txMin, txMax);
    }
    double tyMin = (min[1] - ray.origin[1]) / ray.direction[1];
    double tyMax = (max[1] - ray.origin[1]) / ray.direction[1];
    if (tyMin > tyMax)
    {
        std::swap(tyMin, tyMax);
    }
    if ((txMin > tyMax) || (tyMin > txMax))
    {
        return false;
    }
    double tzMin = (min[2] - ray.origin[2]) / ray.direction[2];
    double tzMax = (max[2] - ray.origin[2]) / ray.direction[2];
    if (tzMin > tzMax)
    {
        std::swap(tzMin, tzMax);
    }
    if ((txMin > tzMax) || (tzMin > txMax) || (tyMin > tzMax) || (tzMin > tyMax))
    {
        return false;
    }
    return true;
};

void BoundingBox::scale(double s, Vector &center)
{
    min = (min - center) * s + center;
    max = (max - center) * s + center;
};

void BoundingBox::translate(Vector &t)
{
    min = min + t;
    max = max + t;
};
