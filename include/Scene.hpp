#pragma once

#include "Sphere.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include <vector>

class Scene
{
public:
    std::vector<Sphere> spheres;
    double intensity;
    Vector lightSource;

    void add(Sphere s);

    bool intersect(Ray &ray, Vector &P, Vector &N, int &objectIndex);

    double lightVisibility(Vector &P);

    Vector getColor(Ray &ray, int depth);

};


