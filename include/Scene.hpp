#pragma once

#include "Sphere.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include <vector>

class Scene
{
public:
    std::vector<const Object*> objects;

    void addObject(const Object &obj);

    bool intersect(Ray &ray, Vector &P, Vector &N, int &objectIndex, double &t);

    Vector getColor(Ray &ray, int depth, bool isIndirect = false);

};


