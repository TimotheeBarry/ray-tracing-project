
#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "../stb/stb_image.h"

#include "../include/Vector.hpp"
#include "../include/Ray.hpp"
#include "../include/Sphere.hpp"
#include "../include/Scene.hpp"
#include "../include/Functions.hpp"
#include "../include/Constants.hpp"

#include <iostream>

int main()
{
    int W = 500;
	int H = 500;
	double alpha = 80 * (PI) / 180;

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
