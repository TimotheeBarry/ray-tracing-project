#include "../include/Functions.hpp"
#include "../include/Constants.hpp"
#include <cmath>

double sqr(double x)
{
    return x * x;
}

Vector gammaCorrection(Vector color)
{
    return pow(color / 255, 1 / GAMMA)*255;
}

Vector computeColor(Vector albedo, Vector lightVector, Vector N, double intensity)
{
    return albedo * intensity * (std::max(0., dot(lightVector.normalized(), N)) / (4 * sqr(PI) * lightVector.norm2()));
}
