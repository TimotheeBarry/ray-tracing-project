#pragma once

#include "Vector.hpp"
#include "Ray.hpp"
#include "TriangleIndices.hpp"
#include <vector>

class TriangleIndices;

class BoundingBox
{
public:
    BoundingBox() : min(Vector()), max(Vector()){};
    ~BoundingBox() {}
    Vector min, max;

    bool intersect(Ray &ray) const;
    void scale(double s, Vector &center);
    void translate(Vector &t);
    // void computeDimensions(std::vector<Vector> &vertices);
};