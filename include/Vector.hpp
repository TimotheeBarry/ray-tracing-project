#pragma once
#include <string>

class Vector
{
public:
    explicit Vector(double x = 0, double y = 0, double z = 0);
    double &operator[](int i);
    double operator[](int i) const;
    Vector &operator+=(const Vector &v);
    Vector operator/=(double b);
    double normSquared() const;
    double norm() const;
    Vector normalized();
    void normalize();
    void clip(double a, double b);

    std::string toString() const;
    Vector rotate(double angle, const Vector &axis);

    double coord[3];
};

Vector operator+(const Vector &a, const Vector &b);
Vector operator-(const Vector &a, const Vector &b);
Vector operator*(const Vector &a, double b);
Vector operator*(double a, const Vector &b);
Vector operator*(const Vector &a, const Vector &b);
Vector operator/(const Vector &a, double b);
Vector pow(const Vector &a, double b);
Vector cross(const Vector &a, const Vector &b);
Vector generateRandomUniformVector();
double dot(const Vector &a, const Vector &b);
Vector generateRandomCosineVector(Vector &N);
Vector generateRandomBlinnPhongVector(Vector &wr, double alpha);
