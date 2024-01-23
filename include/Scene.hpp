#pragma once

#include "Sphere.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include <vector>

class Scene
{
public:
    std::vector<Sphere> spheres;

    void add(Sphere s);

    bool intersect(Ray ray, Vector &P, Vector &N, int &objectIndex);

    double lightVisibility(Vector P, Vector L);

    Vector getColor(Vector color, Vector lightSource, Ray ray, double intensity, int depth);

};