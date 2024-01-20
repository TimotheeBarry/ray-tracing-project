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

const double epsilon = 1e-6;

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

	void clip(double a, double b)
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
	Vector center;
	double radius;
	Vector albedo;
	double reflectance;
	double opacity;
	double refractiveIndex;

	// Constructor with default parameters
	Sphere(const Vector &center = Vector(), double radius = 0.0,
		   const Vector &albedo = Vector(), double reflectance = 0.0,
		   double opacity = 1.0, double refractiveIndex = 1.0)
		: center(center), radius(radius), albedo(albedo),
		  reflectance(reflectance), opacity(opacity), refractiveIndex(refractiveIndex) {}

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

	Vector findOutPoint(Ray ray)
	{
		double t = this->intersectionDistance(ray);
		return ray.origin + t * ray.direction;
	}
};

Vector computeColor(Vector albedo, Vector lightVector, Vector N, double intensity)
{
	return albedo * intensity * (std::max(0., dot(lightVector.normalized(), N)) / (4 * std::pow(M_PI, 2) * lightVector.norm2()));
}

class Scene
{
public:
	std::vector<Sphere> spheres;

	void add(Sphere s)
	{
		spheres.push_back(s);
	}

	bool intersect(Ray ray, Vector &P, Vector &N, int &objectIndex)
	{
		double t = 1e10;
		bool intersect = false;

		for (int i = 0; i < spheres.size(); i++)
		{
			double tTemp = spheres[i].intersectionDistance(ray);
			if (tTemp > 0 && tTemp < t)
			{
				t = tTemp;
				P = ray.origin + t * ray.direction;
				N = (P - spheres[i].center).normalized();
				objectIndex = i;
				intersect = true;
			}
		}
		return intersect;
	}

	bool isLightVisible(Vector P, Vector L)
	{
		Vector lightVector = L - P;
		Ray lightRay(P + lightVector.normalized() * epsilon, lightVector.normalized());
		bool isLightVisible(true);
		for (int j = 0; j < spheres.size(); j++)
		{
			double tTemp = spheres[j].intersectionDistance(lightRay);
			if (tTemp > 0 && tTemp < std::sqrt(lightVector.norm2()))
			{
				// si il y a une intersection entre le point et la lumiere avec une autre sphere
				isLightVisible = false;
				break;
			}
		}
		return isLightVisible;
	}

	Vector getColor(Vector color, Vector L, Ray ray, double intensity, int depth)
	{
		const int maxDepth = 5;
		if (depth > maxDepth)
		{
			return color;
		}
		Vector P, N;
		int intersectIndex = -1;
		bool intersect = this->intersect(ray, P, N, intersectIndex);
		if (intersect)
		{
			Vector lightVector = L - P;

			bool isLightVisible = this->isLightVisible(P, L);
			if (!isLightVisible)
			{
				return color;
			}
			// si la sphere est reflechissante
			if (spheres[intersectIndex].reflectance > 0)
			{
				Vector reflexionVector = ray.direction - 2 * dot(ray.direction, N) * N;
				Ray reflexionRay(P, reflexionVector);
				Vector objectColor = computeColor(spheres[intersectIndex].albedo, lightVector, N, intensity);
				return (1-spheres[intersectIndex].reflectance)*objectColor + spheres[intersectIndex].reflectance*getColor(color, L, reflexionRay, intensity, depth + 1);
			}
			else if (spheres[intersectIndex].opacity < 1)
			{
				// entrée dans la sphère
				bool isEntering = dot(ray.direction, N) < 0;
				double n1, n2;
				if (isEntering)
				{
					n1 = 1.0; 
					n2 = spheres[intersectIndex].refractiveIndex;
				}
				else
				{
					n1 = spheres[intersectIndex].refractiveIndex;
					n2 = 1.0; 
				}
				double n = n1 / n2;
				double cosThetaI = dot(ray.direction, N);
				Vector tN = (std::sqrt(1 - sqr(n) * (1 - sqr(cosThetaI))) * N).normalized();
				if (isEntering)
				{
					tN = (-1)*tN;
				}
				Vector tT = (n * ray.direction).normalized();

				Ray refractionRay = Ray(P + tN * epsilon, tN + tT);
				
				return getColor(color, L, refractionRay, intensity, depth + 1);
				
			}
			else
			{
				Vector objectColor = computeColor(spheres[intersectIndex].albedo, lightVector, N, intensity);
				return color + objectColor;
			}
		}
		return color;
	}
};

int main()
{
	int W = 2000;
	int H = 2000;
	double alpha = 80 * (M_PI) / 180;

	Vector center(0.2, 0.1, 0.);

	std::vector<unsigned char> image(W * H * 3, 0);
	const double intensity = 3e7;

	Vector O(0, 0, 55);	   // camera origin
	Vector L(-10, 20, 40); // light position (point light)
	Scene scene = Scene();
	Sphere sphere1(Vector(0, 0, 0), 10, Vector(0, 0.6, 0), 0.0, 0.0, 2.417);
	Sphere floor(Vector(0, -1000, 0), 990, Vector(1, 1, 1), 0.5);
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(1, 0.1, 0.1), 0.5);
	Sphere wallBack(Vector(0, 0, -1000), 940, Vector(0.1, 0.1, 1),0.5);
	Sphere wallFront(Vector(0, 0, 1000), 940, Vector(0.1, 0.1, 1), 0.5);
	Sphere wallLeft(Vector(-1000, 0, 0), 940, Vector(0.5, 0.1, 1),0.5);
	Sphere wallRight(Vector(1000, 0, 0), 940, Vector(0.1, 0.1, 5),0.5);

	scene.add(sphere1);
	scene.add(floor);
	scene.add(ceiling);
	scene.add(wallLeft);
	scene.add(wallRight);
	scene.add(wallBack);
	scene.add(wallFront);

	Vector albedo = scene.spheres[0].albedo;

#pragma omp parallel for
	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			Vector P, N;
			Vector vector(j - W / 2 + 0.5, -i + H / 2 - 0.5, -W / (2 * std::tan(alpha / 2)));
			Ray ray(O, vector.normalized());

			Vector color = scene.getColor(albedo, L, ray, intensity, 0);
			color.clip(0, 255);
			color = gammaCorrection(color);

			image[(i * W + j) * 3 + 0] = color[0]; // RED
			image[(i * W + j) * 3 + 1] = color[1]; // GREEN
			image[(i * W + j) * 3 + 2] = color[2]; // BLUE
		}
	}
	stbi_write_png("image.png", W, H, 3, &image[0], 0);

	return 0;
}