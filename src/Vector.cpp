#include "../include/Vector.hpp"
#include "../include/Functions.hpp"
#include "../include/Constants.hpp"
#include <cmath>
#include <sstream>

Vector::Vector(double x, double y, double z)
{
    coord[0] = x;
    coord[1] = y;
    coord[2] = z;
}
double &Vector::operator[](int i) { return coord[i]; }
double Vector::operator[](int i) const { return coord[i]; }

Vector &Vector::operator+=(const Vector &v)
{
    coord[0] += v[0];
    coord[1] += v[1];
    coord[2] += v[2];
    return *this;
}

Vector Vector::operator/=(double b)
{
    if (b == 0)
        throw "Division by zero";

    coord[0] /= b;
    coord[1] /= b;
    coord[2] /= b;
    return *this;
}

double Vector::normSquared() const
{
    return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
}

double Vector::norm() const
{
    return sqrt(normSquared());
}

Vector Vector::normalized()
{
    double norm = sqrt(normSquared());
    return Vector(coord[0] / norm, coord[1] / norm, coord[2] / norm);
}

void Vector::normalize()
{
    double norm = sqrt(normSquared());
    coord[0] /= norm;
    coord[1] /= norm;
    coord[2] /= norm;
}

void Vector::clip(double a, double b)
{
    coord[0] = std::max(a, std::min(b, coord[0]));
    coord[1] = std::max(a, std::min(b, coord[1]));
    coord[2] = std::max(a, std::min(b, coord[2]));
}

Vector Vector::generateRandomCosineVector()
{
    // génère un vecteur aléatoire suivant une loi cosinus selon le vecteur courant
    double r1 = uniform(gen);
    double r2 = uniform(gen);
    Vector randomVector(cos(2 * PI * r1) * sqrt(1 - r2), sin(2 * PI * r1) * sqrt(1 - r2), sqrt(r2));

    // on créé un repère local (u, v, N) à partir de N
    Vector u = cross(*this, generateRandomUniformVector()).normalized();
    Vector v = cross(*this, u).normalized();

    // on retourne le vecteur aléatoire dans le repère global
    return randomVector[0] * u + randomVector[1] * v + randomVector[2] * (*this);
}

std::string Vector::toString()
{
    return "(" + std::to_string(coord[0]) + ", " + std::to_string(coord[1]) + ", " + std::to_string(coord[2]) + ")";
}
Vector operator+(const Vector &a, const Vector &b)
{
    return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector &a, const Vector &b)
{
    return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector &a, double b)
{
    return Vector(a[0] * b, a[1] * b, a[2] * b);
}
Vector operator*(double a, const Vector &b)
{
    return Vector(a * b[0], a * b[1], a * b[2]);
}

Vector operator*(const Vector &a, const Vector &b)
{
    return Vector(a[0] * b[0], a[1] * b[1], a[2] * b[2]);
}

double dot(const Vector &a, const Vector &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector operator/(const Vector &a, double b)
{
    if (b == 0)
        throw "Division by zero";

    return Vector(a[0] / b, a[1] / b, a[2] / b);
}
Vector pow(const Vector &a, double b)
{
    return Vector(pow(a[0], b), pow(a[1], b), pow(a[2], b));
}

Vector cross(const Vector &a, const Vector &b)
{
    return Vector(a[1] * b[2] - a[2] * b[1],
                  a[2] * b[0] - a[0] * b[2],
                  a[0] * b[1] - a[1] * b[0]);
}

Vector generateRandomUniformVector()
{
    // génère un vecteur aléatoire dont chaque coordonnée est comprise entre -1 et 1 (uniformément)
    return Vector(uniform(gen) - .5, uniform(gen) - .5, uniform(gen) - .5);
}
