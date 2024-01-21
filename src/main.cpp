
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

std::vector<unsigned char> subSampleImage(std::vector<unsigned char> image, int W, int H, int subSamplingFactor)
{
	if (subSamplingFactor == 1)
	{
		return image;
	}
	int W_bis = W / subSamplingFactor;
	int H_bis = H / subSamplingFactor;

	std::vector<unsigned char> subSampledImage(W_bis * H_bis * 3, 0);

	for (int i = 0; i < H; i += subSamplingFactor)
	{
		for (int j = 0; j < W; j += subSamplingFactor)
		{
			int index = (i * W + j) * 3;
			int subSampledIndex = (i / subSamplingFactor * W_bis + j / subSamplingFactor) * 3;
			for (int k = 0; k < 3; k++)
			{
				int sum = 0;
				for (int l = 0; l < subSamplingFactor; l++)
				{
					for (int m = 0; m < subSamplingFactor; m++)
					{
						sum += image[index + (l * W + m) * 3 + k];
					}
				}
				subSampledImage[subSampledIndex + k] = sum / (subSamplingFactor * subSamplingFactor);
			}
		}
	}

	return subSampledImage;
}

int main()
{
	int W = 2048;
	int H = 2048;
	int subSamplingFactor = 1;
	double alpha = 80 * (PI) / 180;

	Vector center(0.2, 0.1, 0.);

	std::vector<unsigned char> image(W * H * 3, 0);
	const double intensity = 3e7;

	Vector O(0, 0, 55);	   // camera origin
	Vector L(-10, 20, 40); // light position (point light)
	Scene scene = Scene();
	Sphere sphere1(Vector(0, 0, 0), 10, Vector(0, 0.6, 0), 0.0, 0.5, 2.4);
	Sphere sphere2(Vector(10, 10, -10), 5, Vector(0, 1, 1), 0.0, 0.5, 1.33);
	Sphere sphere3(Vector(-30, -2, -10), 8, Vector(1, 0, 1), 1.0);
	Sphere floor(Vector(0, -1000, 0), 990, Vector(1, 1, 1), 0.5);
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(1, 0.1, 0.1), 0.1);
	Sphere wallBack(Vector(0, 0, -1000), 940, Vector(0.1, 0.1, 1), 0.1);
	Sphere wallFront(Vector(0, 0, 1000), 940, Vector(0.1, 0.1, 1), 0.1);
	Sphere wallLeft(Vector(-1000, 0, 0), 940, Vector(0.5, 0.1, 1), 0.1);
	Sphere wallRight(Vector(1000, 0, 0), 940, Vector(0.1, 0.1, 5), 0.1);

	scene.add(sphere1);
	scene.add(sphere2);
	scene.add(sphere3);
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
			// if (i == 250)
			// {
			// 	std::cout << color[0] << color[1] << color[2] << std::endl;
			// }

			image[(i * W + j) * 3 + 0] = color[0]; // RED
			image[(i * W + j) * 3 + 1] = color[1]; // GREEN
			image[(i * W + j) * 3 + 2] = color[2]; // BLUE
		}
	}

	std::vector<unsigned char> subSampledImage = subSampleImage(image, W, H, subSamplingFactor);

	stbi_write_png("image.png", W / subSamplingFactor, H / subSamplingFactor, 3, &subSampledImage[0], 0);

	return 0;
}
