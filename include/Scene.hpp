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
    Vector getIndirectColor(Ray &ray, Vector &P, Vector &N, int depth, Vector &albedo, Material *mat);
    Vector getTransparentColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, double refractiveIndex);
    Vector getDiffusedColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, Vector &albedo, Material *mat);
    Vector getReflectionColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect);
    Vector getMirrorColor(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, double reflectance, Vector &albedo, Material *mat);
    Vector computeColorFromBRDF(Ray &ray, Vector &P, Vector &N, int depth, bool isIndirect, Vector &albedo, Material *mat);
};
