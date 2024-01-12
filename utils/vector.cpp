#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../stb/stb_image.h"

static inline double sqr(double x) { return x * x; }

class Vector
{
public:
	explicit Vector(double x = 0, double y = 0, double z = 0)
	{
		coord[0] = x;
		coord[1] = y;
		coord[2] = z;
	}
	double &operator[](int i) { return coord[i]; }
	double operator[](int i) const { return coord[i]; }

	Vector &operator+=(const Vector &v)
	{
		coord[0] += v[0];
		coord[1] += v[1];
		coord[2] += v[2];
		return *this;
	}

	double norm2() const
	{
		return sqr(coord[0]) + sqr(coord[1]) + sqr(coord[2]);
	}

	double coord[3];
};

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