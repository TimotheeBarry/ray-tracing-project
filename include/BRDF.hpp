#pragma once

#include "Vector.hpp"
#include "Constants.hpp"

class BRDF
{
public:
    virtual Vector f(Vector &wi, Vector &wo, Vector &N) const = 0;
    virtual ~BRDF() {}
    virtual Vector getAlbedo() const = 0;
};

class LambertianBRDF : public BRDF
{
public:
    Vector albedo;
    LambertianBRDF(Vector albedo) : albedo(albedo) {}
    Vector f(Vector &wi, Vector &wo, Vector &N) const override;
    Vector getAlbedo() const override { return albedo; }


};

class BlinnPhongBRDF : public BRDF
{
public:
    double alpha;
    Vector albedo;
    BlinnPhongBRDF(double alpha, Vector albedo) : alpha(alpha), albedo(albedo) {}
    Vector f(Vector &wi, Vector &wo, Vector &N) const override;
    Vector getAlbedo() const override { return albedo; }
};

class Transparent : public BRDF
{
public:
    double refractiveIndex;
    Transparent(double refractiveIndex) : refractiveIndex(refractiveIndex) {}
    Vector f(Vector &wi, Vector &wo, Vector &N) const override;
    Vector getAlbedo() const override { return Vector(0, 0, 0); }
};

class Mirror : public BRDF
{
public:
    double reflectance;
    Vector albedo;
    Mirror(double reflectance, Vector albedo) : reflectance(reflectance), albedo(albedo) {}
    Vector f(Vector &wi, Vector &wo, Vector &N) const override;
    Vector getAlbedo() const override { return albedo; }
};
