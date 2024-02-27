#pragma once

#include <tuple>
#include "Vector.hpp"
#include "Ray.hpp"
#include "Object.hpp"

class Sphere: public Object
{
public:
    Vector center;
    double radius;
    Vector albedo;
    double reflectance;
    double opacity;
    double refractiveIndex;
    double lightIntensity;

    // Constructor with default parameters
    Sphere(const Vector &center = Vector(),
           double radius = 0.0,
           const Vector &albedo = Vector(),
           double reflectance = 0.0,
           double opacity = 1.0,
           double refractiveIndex = 1.0,
           double lightIntensity = 0.0);

    double intersect(Ray &ray, Vector &P, Vector &N) const override;

};