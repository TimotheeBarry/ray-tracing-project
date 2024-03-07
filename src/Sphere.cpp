#include "../include/Vector.hpp"
#include "../include/Sphere.hpp"
#include "../include/Ray.hpp"
#include <cmath>

Sphere::Sphere(
    const Vector &center,
    double radius,
    const Vector &albedo,
    double reflectance,
    double isTransparent,
    double refractiveIndex,
    double lightIntensity) : center(center),
                             radius(radius),
                             albedo(albedo),
                             reflectance(reflectance),
                             opacity(isTransparent),
                             refractiveIndex(refractiveIndex),
                             lightIntensity(lightIntensity)
{
}

double Sphere::intersect(Ray &ray, Vector &P, Vector &N) const
{
    double a = 1;
    double b = 2 * dot(ray.direction, ray.origin - center);
    double c = (ray.origin - center).normSquared() - std::pow(radius, 2);
    double delta = std::pow(b, 2) - 4 * a * c;

    if (delta < 0)
    {
        return -1;
    }
    double t1 = (-b - std::sqrt(delta)) / (2 * a);
    double t2 = (-b + std::sqrt(delta)) / (2 * a);

    if (t1 < 0 && t2 < 0)
    {
        return -1;
    }

    double t = t1 < t2 ? t1 : t2;
    P = ray.origin + ray.direction * t;
    N = (P - center).normalized();
    return t;
}

bool Sphere::fastIntersect(Ray &ray) const
{
    double a = 1;
    double b = 2 * dot(ray.direction, ray.origin - center);
    double c = (ray.origin - center).normSquared() - std::pow(radius, 2);
    double delta = std::pow(b, 2) - 4 * a * c;

    if (delta < 0)
    {
        return false;
    }
    double t1 = (-b - std::sqrt(delta)) / (2 * a);
    double t2 = (-b + std::sqrt(delta)) / (2 * a);
    double t = t1 < t2 ? t1 : t2;
    if (t < 0 || t > ray.length)
    {
        return false;
    }
    return true;
}
