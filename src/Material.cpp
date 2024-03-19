#include "../include/Material.hpp"

Vector Diffuse::brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const
{
    return albedo * (1 / PI);
}

Vector BlinnPhong::brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const
{
    Vector wr = (2 * dot(N, wo) * N - wo).normalized();
    double cosTheta = std::max(0.0, dot(wi, wr));
    return albedo * ((1 - shininess) * (1. / PI) + shininess * (alpha + 2) / (2. * PI) * pow(cosTheta, alpha));
}

Vector Transparent::brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const
{
    return Vector(0, 0, 0);
}

Vector Mirror::brdf(Vector &albedo, Vector &wi, Vector &wo, Vector &N) const
{
    return albedo * (1 / PI);
}