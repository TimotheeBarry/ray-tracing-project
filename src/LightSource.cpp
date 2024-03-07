#include "../include/LightSource.hpp"
#include "../include/Sphere.hpp"
#include "../include/Constants.hpp"

LightSource::LightSource(
    const Vector &center,
    double radius,
    double intensity) : center(center), radius(radius), intensity(intensity) {}

double LightSource::intersect(Ray &ray, Vector &P, Vector &N) const
{
    return Sphere(center, radius).intersect(ray, P, N);
}

bool LightSource::fastIntersect(Ray &ray) const
{
    return Sphere(center, radius).fastIntersect(ray);
}

double LightSource::realIntensity() const
{
    return intensity / (4 * PI * std::pow(radius, 2));
}
