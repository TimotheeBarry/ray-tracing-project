# pragma once

#include "Vector.hpp"
#include "Ray.hpp"

class BoundingBox
{
public:
    BoundingBox() : min(Vector()), max(Vector()) {}
    BoundingBox(Vector min, Vector max) : min(min), max(max) {}
    Vector min, max;
    bool intersect(Ray &ray) const;
    void scale(double s, Vector &center);
    void translate(Vector &t);
    void rotate(double angle, Vector &axis, Vector &center);
};