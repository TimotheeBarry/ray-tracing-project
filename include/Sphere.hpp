#pragma once

#include <tuple>
#include "Vector.hpp"
#include "Ray.hpp"
#include "Object.hpp"
#include "BRDF.hpp"

class Sphere : public Object
{
public:
    Vector center;
    double radius;
    BRDF *brdf;

    // Constructor with default parameters
    Sphere(const Vector &center = Vector(),
           double radius = 0.0,
           BRDF *brdf = nullptr)
        : center(center), radius(radius), brdf(brdf)
    {
        if (brdf == nullptr)
        {
            this->brdf = new LambertianBRDF(Vector(1, 1, 1));
        }
    }

    double intersect(Ray &ray, Vector &P, Vector &N, Vector &albedo) const override;
    bool fastIntersect(Ray &ray) const override;
};