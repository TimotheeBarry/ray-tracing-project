
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
	bool showProgress = true;
	double n = 0.0; // pourcentage de complétion
	int s = 4096;
	int W = s;
	int H = s;
	int subSamplingFactor = 1;
	double alpha = 80 * (PI) / 180;

	Vector center(0.2, 0.1, 0.);

	std::vector<unsigned char> image(W * H * 3, 0);
	const double intensity = 3e7;

	Vector O(0, 0, 55);				 // camera origin
	Vector lightSource(-10, 20, 40); // light position (point light)
	Scene scene = Scene();
	Sphere sphere1(Vector(0, 0, 0), 10, Vector(0, 0.6, 0), 0.0, 0.5, 2.4);
	Sphere sphere2(Vector(10, 10, -10), 10, Vector(0, 1, 1), 0.0, 0.5, 1.33);
	Sphere sphere3(Vector(-30, -2, -10), 8, Vector(1, 0, 1), 1.0);
	Sphere floor(Vector(0, -1000, 0), 990, Vector(1, 1, 1), 0.5);
	Sphere ceiling(Vector(0, 1000, 0), 940, Vector(1, 0.1, 0.1), 0.1);
	Sphere wallBack(Vector(0, 0, 1000), 940, Vector(0.1, 0.1, 1), 0.1);
	Sphere wallFront(Vector(0, 0, -1000), 940, Vector(0.1, 0.1, 1), 1.0);
	Sphere wallLeft(Vector(-1000, 0, 0), 940, Vector(0.5, 0.1, 1), 0.1);
	Sphere wallRight(Vector(1000, 0, 0), 940, Vector(0.1, 0.1, 5), 0.1);

	// scene.add(sphere1);
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

			Vector color = scene.getColor(albedo, lightSource, ray, intensity, 0);

			color.clip(0, 255);
			color = gammaCorrection(color);

			image[(i * W + j) * 3 + 0] = color[0]; // RED
			image[(i * W + j) * 3 + 1] = color[1]; // GREEN
			image[(i * W + j) * 3 + 2] = color[2]; // BLUE

			if (showProgress)
			// afficher le pourcentage de complétion dans une barre de chargement
			{
				double progress = getPercentage(i, j, H, W);
				if (progress > (n / 100))
				{
					// std::cout << n << "% done" << std::endl;
					int barWidth = 70;

					std::cout << "[";
					int pos = barWidth * progress;
					for (int i = 0; i < barWidth; ++i)
					{
						if (i < pos)
							std::cout << "=";
						else if (i == pos)
							std::cout << ">";
						else
							std::cout << " ";
					}
					std::cout << "] " << int(progress * 100.0) << " %\r";
					std::cout.flush();
					n += 1;
				}
			}
		}
	}
	if (showProgress)
	{
		std::cout << std::endl;
	}

	std::vector<unsigned char> subSampledImage = subSampleImage(image, W, H, subSamplingFactor);

	stbi_write_png("image.png", W / subSamplingFactor, H / subSamplingFactor, 3, &subSampledImage[0], 0);

	return 0;
}
