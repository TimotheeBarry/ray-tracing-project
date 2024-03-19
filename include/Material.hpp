#pragma once

#include "Vector.hpp"
#include "Constants.hpp"

class Material
{
public:
    virtual Vector brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const = 0;
    virtual ~Material() {}
};

class Diffuse : public Material
{
public:
    Diffuse() {}
    Vector brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const override;
};

class BlinnPhong : public Material
{
public:
    double alpha;
    double shininess;
    BlinnPhong(double alpha, double shininess) : alpha(alpha), shininess(shininess) {}
    Vector brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const override;
};

class Transparent : public Material
{
public:
    double refractiveIndex;
    Transparent(double refractiveIndex) : refractiveIndex(refractiveIndex) {}
    Vector brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const override;
};

class Mirror : public Material
{
public:
    double reflectance;

    Mirror(double reflectance) : reflectance(reflectance) {}
    Vector brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const override;
};
