#include "../include/Vector.hpp"
#include "../include/Functions.hpp"

#include <cmath>

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

double Vector::norm2() const
{
    return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
}

Vector Vector::normalized()
{
    double norm = sqrt(norm2());
    return Vector(coord[0] / norm, coord[1] / norm, coord[2] / norm);
}

void Vector::normalize()
{
    double norm = sqrt(norm2());
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
