#include "../include/BRDF.hpp"

Vector LambertianBRDF::f(Vector &wi, Vector &wo, Vector &N) const
{
    return albedo * (1 / PI);
}

Vector BlinnPhongBRDF::f(Vector &wi, Vector &wo, Vector &N) const
{
    Vector H = (wi + wo).normalized(); // demi-vecteur
    double cosTheta = std::max(0.0, dot(N, H));
    return albedo * (alpha + 8) / (8 * PI) * pow(cosTheta, alpha);
}

Vector Transparent::f(Vector &wi, Vector &wo, Vector &N) const
{
    return Vector(0, 0, 0);
}

Vector Mirror::f(Vector &wi, Vector &wo, Vector &N) const
{
    return albedo * (1 / PI);
}