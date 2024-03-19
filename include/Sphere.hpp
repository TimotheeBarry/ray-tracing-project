#pragma once

#include <tuple>
#include "Vector.hpp"
#include "Ray.hpp"
#include "Object.hpp"
#include "Material.hpp"

class Sphere : public Object
{
public:
    Vector center;
    double radius;
    Vector albedo;
    Material *mat;

    // Constructor with default parameters
    Sphere(const Vector &center = Vector(),
           double radius = 0.0,
           const Vector &albedo = Vector(1, 1, 1),
           Material *mat = nullptr)
        : center(center), radius(radius), albedo(albedo), mat(mat)
    {
        if (mat == nullptr)
        {
            this->mat = new Diffuse();
        }
    }

    double intersect(Ray &ray, Vector &P, Vector &N, Vector &albedo) const override;
    bool fastIntersect(Ray &ray) const override;
};