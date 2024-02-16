#pragma once

#include "Sphere.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include <vector>

class Scene
{
public:
    std::vector<Sphere> spheres;

    void addSphere(Sphere s);

    bool intersect(Ray &ray, Vector &P, Vector &N, int &objectIndex, double &t);

    // double lightVisibility(Vector &P, Vector &L);

    Vector getColor(Ray &ray, int depth, bool isIndirect = false);

};


