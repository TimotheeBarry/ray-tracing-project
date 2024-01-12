#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

#define _USE_MATH_DEFINES
#include <cmath>

#include <iostream>

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

	Vector normalized()
	{
		double norm = sqrt(norm2());
		return Vector(coord[0] / norm, coord[1] / norm, coord[2] / norm);
	}

	void normalize()
	{
		double norm = sqrt(norm2());
		coord[0] /= norm;
		coord[1] /= norm;
		coord[2] /= norm;
	}

	void clamp(double a, double b)
	{
		coord[0] = std::max(a, std::min(b, coord[0]));
		coord[1] = std::max(a, std::min(b, coord[1]));
		coord[2] = std::max(a, std::min(b, coord[2]));
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

Vector gammaCorrection(Vector color)
{
	return Vector(std::pow(color[0] / 255, 1 / 2.2), std::pow(color[1] / 255, 1 / 2.2), std::pow(color[2] / 255, 1 / 2.2)) * 255;
}

class Ray
{
public:
	explicit Ray(Vector O, Vector u)
	{
		origin = O;
		direction = u;
	}

	Vector origin;
	Vector direction;
};

class Sphere
{
public:
	explicit Sphere(Vector C, double R, Vector a)
	{
		center = C;
		radius = R;
		albedo = a;
	}

	Vector center;
	double radius;
	Vector albedo;

	bool intersect(Ray ray, Vector &P, Vector &N)
	{
		double a = 1;
		double b = 2 * dot(ray.direction, ray.origin - center);
		double c = (ray.origin - center).norm2() - std::pow(radius, 2);
		double delta = std::pow(b, 2) - 4 * a * c;

		if (delta < 0)
		{
			return false;
		}
		double t1 = (-b - std::sqrt(delta)) / (2 * a);
		double t2 = (-b + std::sqrt(delta)) / (2 * a);

		if (t1 < 0 && t2 < 0)
		{
			return false;
		}

		double t = t1 < t2 ? t1 : t2;
		P = ray.origin + t * ray.direction;
		N = (P - center).normalized();

		return true;
	}

	double intersectionDistance(Ray ray)
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
};

class Scene
{
public:
	std::vector<Sphere> spheres;

	void add(Sphere s)
	{
		spheres.push_back(s);
	}

	bool intersect(Ray r, Vector &P, Vector &N, int &index)
	{
		double t = 1e10;
		bool intersect = false;

		for (int i = 0; i < spheres.size(); i++)
		{
			double tTemp = spheres[i].intersectionDistance(r);
			if (tTemp > 0 && tTemp < t)
			{
				t = tTemp;
				P = r.origin + t * r.direction;
				N = (P - spheres[i].center).normalized();
				index = i;
				intersect = true;
			}
		}
		return intersect;
	}
};

int main()
{
	int W = 512;
	int H = 512;
	double alpha = 60 * (M_PI) / 180;

	Vector center(0.2, 0.1, 0.);

	std::vector<unsigned char> image(W * H * 3, 0);
	const double intensity = 2e7;

	Vector O(0, 0, 55);	   // camera origin
	Vector L(-10, 20, 40); // light position (point light)
	Scene scene = Scene();
	Sphere sphere1(Vector(-12, -10, 10), 10, Vector(0.3, 0.4, 0.9));
	Sphere sphere2(Vector(12, 0, 0), 10, Vector(0.9, 0.2, 0.3));
	Sphere sphere3(Vector(-40, 20, -50), 30, Vector(0.9, 0.9, 0.1));
	Sphere floor(Vector(0, -10010, 0), 10000, Vector(0.9, 0.8, 0.8));
	Sphere wall1(Vector(0, 0, -10300), 10000, Vector(0.1, 0.8, 0.4));

	scene.add(sphere1);
	scene.add(sphere2);
	scene.add(sphere3);
	scene.add(floor);
	scene.add(wall1);

	Vector albedo = scene.spheres[0].albedo;

#pragma omp parallel for
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			Vector P, N;
			Vector vector(j - W / 2 + 0.5, -i + H / 2 - 0.5, -W / (2 * std::tan(alpha / 2)));
			Ray ray(O, vector.normalized());
			int intersectIndex = -1;
			bool intersect = scene.intersect(ray, P, N, intersectIndex);

			if (intersect)
			{

				Vector lightVector = L - P;
				Vector albedo = scene.spheres[intersectIndex].albedo;
				// std::cout << intersectIndex << ": "<<albedo[0] << " " << albedo[1] << " " << albedo[2] << std::endl;
				Vector color = albedo * intensity * (std::max(0., dot(lightVector.normalized(), N)) / (4 * std::pow(M_PI, 2) * lightVector.norm2()));
				color.clamp(0, 255);
				color = gammaCorrection(color);

				image[(i * W + j) * 3 + 0] = color[0]; // RED
				image[(i * W + j) * 3 + 1] = color[1]; // GREEN
				image[(i * W + j) * 3 + 2] = color[2]; // BLUE
			}
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}