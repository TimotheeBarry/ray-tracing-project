#pragma once

#include "Vector.hpp"
#include "Ray.hpp"

class Sphere
{
public:
    Vector center;
    double radius;
    Vector albedo;
    double reflectance;
    double opacity;
    double refractiveIndex;

    // Constructor with default parameters
    Sphere(const Vector &center = Vector(),
           double radius = 0.0,
           const Vector &albedo = Vector(),
           double reflectance = 0.0,
           double opacity = 1.0,
           double refractiveIndex = 1.0);

    double intersectionDistance(Ray ray);

    Vector findOutPoint(Ray ray);
};