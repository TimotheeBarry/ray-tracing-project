#pragma once

#include "Sphere.hpp"
#include "Vector.hpp"
#include "Ray.hpp"
#include <vector>

class Scene
{
public:
    std::vector<const Object *> objects;
    void addObject(const Object &obj);
    bool intersect(Ray &ray, Vector &P, Vector &N, int &objectIndex, double &t, Vector &albedo);
    bool intersectObjectOnly(Ray &ray);
    Vector getColor(Ray &ray, int depth, bool isIndirect = false);

private:
    Vector getIndirectColor(Vector &P, Vector &N, int depth, Vector &albedo);
    Vector getTransparentSphereColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, double refractiveIndex);
    Vector getDiffusedColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, Vector &albedo);
    Vector getReflectionColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect);
};
