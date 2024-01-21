#include "../include/Vector.hpp"
#include "../include/Sphere.hpp"
#include "../include/Ray.hpp"
#include <cmath>

Sphere::Sphere(
    const Vector &center,
    double radius,
    const Vector &albedo,
    double reflectance,
    double opacity,
    double refractiveIndex) : center(center),
                              radius(radius),
                              albedo(albedo),
                              reflectance(reflectance),
                              opacity(opacity),
                              refractiveIndex(refractiveIndex) {}

double Sphere::intersectionDistance(Ray ray)
{
    double a = 1;
    double b = 2 * dot(ray.direction, ray.origin - center);
    double c = (ray.origin - center).norm2() - std::pow(radius, 2);
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
    return t;
}
