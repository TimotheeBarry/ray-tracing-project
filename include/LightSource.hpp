#pragma once

#include <tuple>
#include "Vector.hpp"
#include "Ray.hpp"
#include "Object.hpp"

class LightSource: public Object
{
public:
    Vector center;
    double radius;
    double intensity;

    // Constructor with default parameters
    LightSource(const Vector &center = Vector(),
           double radius = 0.0,
           double intensity = 0.0);

    double intersect(Ray ray, Vector &P, Vector &N) const override;

    double realIntensity() const;

};